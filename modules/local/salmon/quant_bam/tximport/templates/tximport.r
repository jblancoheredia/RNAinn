#!/usr/bin/env Rscript --vanilla

# Script for importing and processing transcript-level quantifications.
# Modified to accept direct files instead of directories
# All property access methods standardized to [[]] notation for consistency

# Loading required libraries
library(SummarizedExperiment)
library(tximport)

################################################
################################################
## Functions                                  ##
################################################
################################################

#' Build a table from a SummarizedExperiment object
#'
#' This function takes a SummarizedExperiment object and a specific slot name to extract
#' assay data. It then combines the first two columns of the rowData with the specified
#' assay data slot into a new data table.
#'
#' @param se.obj A SummarizedExperiment object from which to build the table.
#' @param slot The name of the slot in the assays list from which to extract data.
#'
#' @return A data frame combining the first two columns of the rowData with the assay data from the specified slot.

build_table <- function(se.obj, slot) {
    data.frame(cbind(rowData(se.obj)[,1:2], assays(se.obj)[[slot]]), check.names = FALSE)
}

#' Write a table to a file from a SummarizedExperiment object with given parameters
#'
#' This function generates a table from a SummarizedExperiment object using specified parameters
#' and writes the resulting table to a file. The file name is constructed using a prefix and a
#' suffix from the parameters, and the table is written with tab separation, without quoting text,
#' and without row names.
#'
#' @param params A list containing the parameters needed for file generation and table writing.
#' The list should include:
#' - `obj`: A SummarizedExperiment object from which to build the table.
#' - `slot`: The name of the slot in the assays list from which to extract data.
#' - `suffix`: Suffix to use for generating the file name.
#'
#' @return NULL The function is called for its side effect of writing a file and does not return anything.

write_se_table <- function(params, prefix) {
    file_name <- paste0(prefix, ".", params[["suffix"]])
    write.table(build_table(params[["obj"]], params[["slot"]]), file_name,
                sep="\t", quote=FALSE, row.names = FALSE)
}

#' Read Transcript Metadata from a Given Path
#'
#' This function reads transcript metadata from a specified file path. The file is expected to
#' be a tab-separated values file without headers, containing transcript information. The function
#' checks if the file is empty and stops execution with an error message if so.
#'
#' @param tinfo_path The file path to the transcript information file.
#'
#' @return A list containing three elements:
#' - `tx_info`: A data frame with transcript IDs, gene IDs, and gene names, indexed by transcript IDs.
#' - `gene`: A data frame with unique gene IDs and gene names.
#' - `tx2gene`: A data frame mapping transcript IDs to gene IDs.

read_transcript_info <- function(tinfo_path){
    info <- file.info(tinfo_path)
    if (info[["size"]] == 0) {
        stop("tx2gene file is empty")
    }

    transcript_info <- read.csv(tinfo_path, sep="\t", header = TRUE,
                                col.names = c("tx", "gene_id", "gene_name"))

    extra <- setdiff(rownames(txi[[1]]), as.character(transcript_info[["tx"]]))
    transcript_info <- rbind(transcript_info, data.frame(tx=extra, gene_id=extra, gene_name=extra))
    transcript_info <- transcript_info[match(rownames(txi[[1]]), transcript_info[["tx"]]), ]
    rownames(transcript_info) <- transcript_info[["tx"]]

    list(tx_info = transcript_info,
        gene = unique(transcript_info[,2:3]),
        tx2gene = transcript_info[,1:2])
}

#' Create a SummarizedExperiment Object
#'
#' Constructs a SummarizedExperiment object using provided matrices for counts, abundance, and length,
#' along with metadata for columns and rows.
#'
#' @param counts A matrix or DataFrame containing counts data
#' @param abundance A matrix or DataFrame containing abundance data
#' @param length A matrix or DataFrame containing feature lengths
#' @param col_data A DataFrame containing sample-level metadata
#' @param row_data A DataFrame containing feature-level metadata
#'
#' @return A SummarizedExperiment object containing the supplied data and metadata.

create_summarized_experiment <- function(counts, abundance, length, col_data, row_data) {
    SummarizedExperiment(assays = list(counts = counts, abundance = abundance, length = length),
        colData = col_data,
        rowData = row_data)
}

################################################
################################################
## Main script starts here                    ##
################################################
################################################

# Instead of searching for files, use the directly provided quant.sf file
quant_file <- '$quant'
sample_name <- basename(dirname(quant_file))
if (sample_name == ".") {
    # If the file is in the current directory, use the file basename without extension
    sample_name <- tools::file_path_sans_ext(basename(quant_file))
}

# Create a named list with just one quantification file
fns <- c(quant_file)
names(fns) <- sample_name
dropInfReps <- '$quant_type' == "kallisto"

# Import transcript-level quantifications
txi <- tximport(fns, type = '$quant_type', txOut = TRUE, dropInfReps = dropInfReps)

# Read transcript and sample data
transcript_info <- read_transcript_info('$tx2gene')

# Make coldata just to appease the summarizedexperiment
coldata <- data.frame(files = fns, names = names(fns))
rownames(coldata) <- coldata[["names"]]

# Create initial SummarizedExperiment object
se <- create_summarized_experiment(txi[["counts"]], txi[["abundance"]], txi[["length"]],
    DataFrame(coldata), transcript_info[["tx_info"]])

# Setting parameters for writing tables
params <- list(
    list(obj = se, slot = "abundance", suffix = "transcript_tpm.tsv"),
    list(obj = se, slot = "counts", suffix = "transcript_counts.tsv"),
    list(obj = se, slot = "length", suffix = "transcript_lengths.tsv")
)

# Process gene-level data if tx2gene mapping is available
if ("tx2gene" %in% names(transcript_info) && !is.null(transcript_info[["tx2gene"]])) {
    tx2gene <- transcript_info[["tx2gene"]]
    gi <- summarizeToGene(txi, tx2gene = tx2gene)
    gi.ls <- summarizeToGene(txi, tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
    gi.s <- summarizeToGene(txi, tx2gene = tx2gene, countsFromAbundance = "scaledTPM")

    gene_info <- transcript_info[["gene"]][match(rownames(gi[[1]]), transcript_info[["gene"]][["gene_id"]]),]
    rownames(gene_info) <- NULL
    col_data_frame <- DataFrame(coldata)

    # Create gene-level SummarizedExperiment objects
    gse <- create_summarized_experiment(gi[["counts"]], gi[["abundance"]], gi[["length"]],
        col_data_frame, gene_info)
    gse.ls <- create_summarized_experiment(gi.ls[["counts"]], gi.ls[["abundance"]], gi.ls[["length"]],
        col_data_frame, gene_info)
    gse.s <- create_summarized_experiment(gi.s[["counts"]], gi.s[["abundance"]], gi.s[["length"]],
        col_data_frame, gene_info)

    params <- c(params, list(
        list(obj = gse, slot = "length", suffix = "gene_lengths.tsv"),
        list(obj = gse, slot = "abundance", suffix = "gene_tpm.tsv"),
        list(obj = gse, slot = "counts", suffix = "gene_counts.tsv"),
        list(obj = gse.ls, slot = "counts", suffix = "gene_counts_length_scaled.tsv"),
        list(obj = gse.s, slot = "counts", suffix = "gene_counts_scaled.tsv")
    ))
}

# Writing tables for each set of parameters
prefix <- ''
if ('$task.ext.prefix' != 'null'){
    prefix = '$task.ext.prefix'
} else if ('$meta.id' != 'null'){
    prefix = '$meta.id'
}

done <- lapply(params, function(param) write_se_table(param, prefix))

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste(prefix, "R_sessionInfo.log", sep = '.'))
citation("tximeta")
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
tximeta.version <- as.character(packageVersion('tximeta'))

writeLines(
    c(
        '"${task.process}":',
        paste('    bioconductor-tximeta:', tximeta.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################