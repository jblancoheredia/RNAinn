process KALLISTO_SE_GENE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-summarizedexperiment:1.32.0--r43hdfd78af_0' :
        'biocontainers/bioconductor-summarizedexperiment:1.32.0--r43hdfd78af_0' }"

    input:
    tuple val(meta) , path(counts_gene), path(tpm_gene), path(tx2gene)

    output:
    tuple val(meta), path("*.rds")              , emit: rds
    tuple val(meta), path("*.R_sessionInfo.log"), emit: log
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def matrix_files = "${counts_gene} ${tpm_gene}" 
    def coldata = ""
    def rowdata_path = tx2gene ? "${tx2gene}" : ''
    """
    #!/usr/bin/env Rscript

    # Written by Lorena Pantano and revised for flexibility in handling assays
    # Modified by Juan Blanco Heredia blancoj@mskcc.org for mskcc/RNAinn pipeline. 
    # Released under the MIT license.

    library(SummarizedExperiment)

    #' Flexibly read CSV or TSV files
    #'
    #' @param file Input file
    #' @param header Passed to read.delim()
    #' @param row.names Passed to read.delim()
    #'
    #' @return output Data frame

    read_delim_flexible <- function(file, header = TRUE, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE){
        # Safety check: ensure file exists and is not just a string
        if (!file.exists(file)) {
            stop(paste("File not found:", file))
        }

        ext <- tolower(tail(strsplit(basename(file), split = "\\\\.")[[1]], 1))

        if (ext == "tsv" || ext == "txt") {
            separator <- "\\t"
        } else if (ext == "csv") {
            separator <- ","
        } else {
            # Default to tab if extension is unclear
            warning(paste("Unknown file extension for", file, "- defaulting to tab separator"))
            separator <- "\\t"
        }

        read.delim(
            file,
            sep = separator,
            header = header,
            row.names = row.names,
            check.names = check.names,
            stringsAsFactors = stringsAsFactors
        )
    }

    #' Parse out options from a string without recourse to optparse
    #'
    #' @param x Long-form argument list like --opt1 val1 --opt2 val2
    #'
    #' @return named list of options and values similar to optparse

    parse_args <- function(x){
        if (x == "") {
            return(list())
        }
        
        args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
        if (length(args_list) == 0) {
            return(list())
        }
        
        args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

        # Ensure the option vectors are length 2 (key/ value) to catch empty ones
        args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

        parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
        parsed_args[! is.na(parsed_args)]
    }

    #' Find First Column Containing All Specified Entries
    #'
    #' This function searches through each column of a given data frame to find the
    #' first column that contains all of the specified entries in a vector. If such
    #' a column is found, the name of the column is returned. If no column matches,
    #' the function will try a more flexible matching approach before giving up.
    #'
    #' @param namesVector A character vector containing the names to be matched.
    #' @param df A data frame within which to search for the column containing all
    #'   names specified in `namesVector`.
    #'
    #' @return The name of the first column in `df` that contains all entries from
    #'   `namesVector`. If no such column exists, the function will try flexible matching.

    findColumnWithAllEntries <- function(namesVector, df) {
        # Try exact matching first
        for (colName in names(df)) {
            if (all(namesVector %in% df[[colName]])) {
                return(colName)
            }
        }
        
        # If exact matching fails, try matching without version numbers (for transcript IDs)
        # Extract base IDs without version numbers (e.g., ENST00000000233 from ENST00000000233.10)
        baseIds <- sub("\\.[0-9]+\$", "", namesVector)
        
        for (colName in names(df)) {
            # Extract base IDs from the column as well
            columnBaseIds <- sub("\\.[0-9]+\$", "", df[[colName]])
            
            # Check if all base IDs are in the column's base IDs
            if (all(baseIds %in% columnBaseIds)) {
                message("Using flexible matching (ignoring version numbers) for column: ", colName)
                return(colName)
            }
        }
        
        # If we get here, try to print some diagnostic info
        cat("First few entries in namesVector: ", paste(head(namesVector), collapse=", "), "\n", file=stderr())
        cat("Available columns in metadata: ", paste(names(df), collapse=", "), "\n", file=stderr())
        
        # If first column in df has values, show a sample
        if (ncol(df) > 0 && length(df[[1]]) > 0) {
            cat("Sample values from first column (", names(df)[1], "): ", 
                paste(head(df[[1]]), collapse=", "), "\n", file=stderr())
        }
        
        stop(paste("No column contains all vector entries (even with flexible matching). First few entries: ", 
                   paste(head(namesVector, 5), collapse=', ')))
    }

    #' Check Matrix Name Uniformity in List
    #'
    #' Verifies if all matrices in a list have identical row and column names.
    #' It returns TRUE if uniformity is found, otherwise FALSE.
    #'
    #' @param matrices List of matrices.
    #' @return Logical indicating uniformity of row and column names.
    #' @keywords matrix

    checkRowColNames <- function(matrices) {
        # Simplify the comparison process
        allEqual <- function(namesList) {
            all(sapply(namesList[-1], function(x) identical(x, namesList[[1]])))
        }

        rowNamesEqual <- allEqual(lapply(matrices, rownames))
        colNamesEqual <- allEqual(lapply(matrices, colnames))

        if ((! rowNamesEqual) || (! colNamesEqual)){
            stop("Rows or columns different among input matrices")
        }
    }

    #' Parse Metadata From File
    #'
    #' Reads metadata from a specified file and processes it to handle duplicate
    #' rows by aggregating them into a single row based on a unique identifier.
    #' The function dynamically identifies the appropriate ID column if not specified.
    #' It is designed to be flexible for processing either column (sample) or row (feature) metadata.
    #'
    #' @param metadata_path Character string specifying the path to the metadata file.
    #' @param ids Vector of identifiers (column names or row names) used to match against metadata columns.
    #' @param metadata_id_col Optional; character string specifying the column name in the metadata
    #'        to be used as the unique identifier. If NULL, the function attempts to
    #'        automatically find a suitable column based on `ids`.
    #'
    #' @return A data frame of processed metadata with duplicate rows aggregated, and row names set to the unique identifier.

    parse_metadata <- function(metadata_path, ids, metadata_id_col = NULL){
        # Check if file exists
        if (!file.exists(metadata_path)) {
            stop(paste("Metadata file not found:", metadata_path))
        }

        metadata <- read_delim_flexible(metadata_path, stringsAsFactors = FALSE, header = TRUE)
        if (is.null(metadata_id_col)){
            tryCatch({
                metadata_id_col <- findColumnWithAllEntries(ids, metadata)
            }, error = function(e) {
                cat("Error finding matching column:", conditionMessage(e), "\n", file=stderr())
                # Try using the first column as a fallback if it has a reasonable overlap with ids
                if (ncol(metadata) > 0) {
                    first_col <- names(metadata)[1]
                    overlap <- sum(ids %in% metadata[[first_col]])
                    overlap_pct <- overlap / length(ids) * 100
                    
                    if (overlap_pct > 50) {  # If more than 50% match, use it
                        cat("Using first column as fallback with", overlap_pct, "% overlap\n", file=stderr())
                        return(first_col)
                    }
                }
                # If we get here, re-throw the error
                stop(e)
            })
        }

        # Remove any all-NA columns
        metadata <- metadata[, colSums(is.na(metadata)) != nrow(metadata)]

        # Map IDs to the metadata IDs (handle version differences)
        baseMetadataIds <- sub("\\.[0-9]+\$", "", metadata[[metadata_id_col]])
        baseIds <- sub("\\.[0-9]+\$", "", ids)
        
        # Create a mapping from baseIds to original ids
        id_mapping <- setNames(ids, baseIds)
        
        # Find which metadata rows match our ids (either directly or by base id)
        matching_rows <- metadata[[metadata_id_col]] %in% ids | 
                         baseMetadataIds %in% baseIds
        
        if (sum(matching_rows) == 0) {
            warning("No matching rows found in metadata")
            # Return empty data frame with same columns
            empty_df <- data.frame(matrix(nrow=length(ids), ncol=ncol(metadata)))
            colnames(empty_df) <- colnames(metadata)
            rownames(empty_df) <- ids
            return(empty_df)
        }
        
        # Subset to only matching rows
        metadata_subset <- metadata[matching_rows,, drop=FALSE]
        
        # Allow for duplicate rows by the id column
        metadata_agg <- aggregate(
            . ~ metadata_subset[[metadata_id_col]],
            data = metadata_subset,
            FUN = function(x) paste(unique(x), collapse = ",")
        )[,-1]
        
        # Set rownames based on original IDs or mapped IDs
        original_ids <- metadata_subset[[metadata_id_col]]
        rownames(metadata_agg) <- original_ids
        
        # Create a result frame with all requested ids
        result <- data.frame(matrix(nrow=length(ids), ncol=ncol(metadata)))
        colnames(result) <- colnames(metadata)
        rownames(result) <- ids
        
        # Fill in the values for matching rows
        for (id in rownames(metadata_agg)) {
            if (id %in% ids) {
                # Direct match
                result[id,] <- metadata_agg[id,]
            } else {
                # Try matching by base id
                base_id <- sub("\\.[0-9]+\$", "", id)
                matching_base_ids <- names(id_mapping)[names(id_mapping) == base_id]
                
                if (length(matching_base_ids) > 0) {
                    for (matching_id in id_mapping[matching_base_ids]) {
                        result[matching_id,] <- metadata_agg[id,]
                    }
                }
            }
        }
        
        return(result)
    }

    ################################################
    ################################################
    ## Main script starts here                    ##
    ################################################
    ################################################

    # Matrices

    args_opt <- parse_args('${task.ext.args ?: ""}')
    matrix_files <- as.list(strsplit('${matrix_files}', ' ')[[1]])

    if ('assay_names' %in% names(args_opt)){
        names(matrix_files) <- unlist(strsplit(args_opt[['assay_names']], ',')[[1]])
    }else{
        names(matrix_files) <- unlist(lapply(matrix_files, tools::file_path_sans_ext))
    }

    # Build and verify the main assays list for the summarisedexperiment

    assay_list <- lapply(matrix_files, function(m){
        cat("Processing matrix file:", m, "\n")
        mat <- read_delim_flexible(m, row.names = 1, stringsAsFactors = FALSE)
        mat[,sapply(mat, is.numeric), drop = FALSE]
    })

    checkRowColNames(assay_list)

    # Construct SummarizedExperiment
    se <- SummarizedExperiment(
        assays = assay_list
    )

    # Add column (sample) metadata if provided

    coldata_file <- '${coldata}'
    if (coldata_file != '' && file.exists(coldata_file)) {
        cat("Processing column metadata file:", coldata_file, "\n")
        coldata <- parse_metadata(
            metadata_path = coldata_file,
            ids = colnames(assay_list[[1]]),
            metadata_id_col = args_opt\$coldata_id_col
        )

        colData(se) <- DataFrame(coldata)
    } else {
        cat("No column metadata file provided or file not found.\n")
    }

    # Add row (feature) metadata if provided

    rowdata_file <- '${rowdata_path}'
    if (rowdata_file != '' && file.exists(rowdata_file)) {
        cat("Processing row metadata file:", rowdata_file, "\n")
        
        # Try to parse row metadata with more debugging
        tryCatch({
            rowdata <- parse_metadata(
                metadata_path = rowdata_file,
                ids = rownames(assay_list[[1]]),
                metadata_id_col = args_opt\$rowdata_id_col
            )
            
            # Only set rowData if we actually got results
            if (nrow(rowdata) > 0) {
                rowData(se) <- DataFrame(rowdata)
                cat("Successfully added row metadata\n")
            } else {
                cat("Warning: Row metadata parsing returned empty results\n", file=stderr())
            }
        }, error = function(e) {
            cat("Warning: Error in row metadata processing - ", conditionMessage(e), "\n", file=stderr())
            cat("Continuing without row metadata\n", file=stderr())
        })
    } else {
        cat("No row metadata file provided or file not found.\n")
    }

    # Write outputs as RDS files

    prefix <- tools::file_path_sans_ext(matrix_files[1])
    if ('${task.ext.prefix ?: ""}' != '') {
        prefix = '${task.ext.prefix}'
    } else if ('${meta.id ?: ""}' != '') {
        prefix = '${meta.id}'
    }

    # Save the SummarizedExperiment object
    output_file <- paste0(prefix, ".SummarizedExperiment.rds")
    saveRDS(se, file = output_file)

    ################################################
    ################################################
    ## R SESSION INFO                             ##
    ################################################
    ################################################

    sink(paste(prefix, "R_sessionInfo.log", sep = '.'))
    citation("SummarizedExperiment")
    print(sessionInfo())
    sink()

    ################################################
    ################################################
    ## VERSIONS FILE                              ##
    ################################################
    ################################################

    r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
    summarizedexperiment.version <- as.character(packageVersion('SummarizedExperiment'))

    writeLines(
        c(
            '"${task.process}":',
            paste('    bioconductor-summarizedexperiment:', summarizedexperiment.version)
        ),
    'versions.yml')
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.SummarizedExperiment.rds
    touch ${prefix}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-summarizedexperiment: \$(Rscript -e "library(SummarizedExperiment); cat(as.character(packageVersion('SummarizedExperiment')))")
    END_VERSIONS
    """
}