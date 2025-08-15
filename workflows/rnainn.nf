/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                      IMPORT PLUGINS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap                                                              } from 'plugin/nf-validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                      IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { FASTQC                                                                        } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                                                       } from '../modules/nf-core/multiqc/main'
include { CAT_FASTQ                                                                     } from '../modules/local/cat/fastq/main'
include { DOWNSAMPLING_COUNT                                                            } from '../modules/local/downsampling/count'
include { DOWNSAMPLING_SEQTK                                                            } from '../modules/local/downsampling/seqtk'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                   IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VARIANTDSCVRY                                                                 } from '../subworkflows/local/variantdscvry'
include { GENEXPRESSION                                                                 } from '../subworkflows/local/genexpression'
include { FUSION_SPLICE                                                                 } from '../subworkflows/local/fusion_splice'
include { UMIPROCESSING                                                                 } from '../subworkflows/local/umiprocessing'
include { DEDUPONLY4RNA                                                                 } from '../subworkflows/local/deduponly4rna'
include { paramsSummaryMultiqc                                                          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                                        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                                                        } from '../subworkflows/local/utils_nfcore_rnainn_pipeline'
include { validateInputSamplesheet                                                      } from '../subworkflows/local/utils_nfcore_rnainn_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                  CREATE CUSTOM CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_bed                                          = Channel.fromPath(params.bed).map                              { it -> [[id:it.Name], it] }.collect()
ch_fai                                          = Channel.fromPath(params.fai).map                              { it -> [[id:it.Name], it] }.collect()
ch_gtf                                          = Channel.fromPath(params.gtf).map                              { it -> [[id:it.Name], it] }.collect()
ch_bwa2                                         = Channel.fromPath(params.bwa2).map                             { it -> [[id:it.Name], it] }.collect()
ch_dict                                         = Channel.fromPath(params.dict).map                             { it -> [[id:it.Name], it] }.collect()
ch_dbsnp                                        = Channel.fromPath(params.dbsnp).map                            { it -> [[id:it.Name], it] }.collect()
ch_fasta                                        = Channel.fromPath(params.fasta).map                            { it -> [[id:it.Name], it] }.collect()
ch_chrgtf                                       = Channel.fromPath(params.chrgtf).map                           { it -> [[id:it.Name], it] }.collect()
ch_refflat                                      = Channel.fromPath(params.refflat).map                          { it -> [[id:it.Name], it] }.collect()
ch_targets                                      = Channel.fromPath(params.targets).map                          { it -> [[id:it.Name], it] }.collect()
ch_rsem_ref                                     = Channel.fromPath(params.rsem_ref).map                         { it -> [[id:it.Name], it] }.collect()
ch_hgnc_ref                                     = Channel.fromPath(params.hgnc_ref).map                         { it -> [[id:it.Name], it] }.collect()
ch_dbsnp_tbi                                    = Channel.fromPath(params.dbsnp_tbi).map                        { it -> [[id:it.Name], it] }.collect()
ch_intervals                                    = Channel.fromPath(params.intervals).map                        { it -> [[id:it.Name], it] }.collect()
ch_hgnc_date                                    = Channel.fromPath(params.hgnc_date).map                        { it -> [[id:it.Name], it] }.collect()
ch_star_index                                   = Channel.fromPath(params.star_index).map                       { it -> [[id:it.Name], it] }.collect()
ch_targets_bed                                  = Channel.fromPath(params.targets_bed).map                      { it -> [[id:it.Name], it] }.collect()
ch_known_indels                                 = Channel.fromPath(params.known_indels).map                     { it -> [[id:it.Name], it] }.collect()
ch_multiqc_logo                                 = Channel.fromPath(params.multiqc_logo).map                     { it -> [[id:it.Name], it] }.collect()
ch_salmon_index                                 = Channel.fromPath(params.salmon_index).map                     { it -> [[id:it.Name], it] }.collect()
ch_kallisto_index                               = Channel.fromPath(params.kallisto_index).map                   { it -> [[id:it.Name], it] }.collect()
ch_rrna_intervals                               = Channel.fromPath(params.rrna_intervals).map                   { it -> [[id:it.Name], it] }.collect()
ch_cnvkit_reference                             = Channel.fromPath(params.cnvkit_reference).map                 { it -> [[id:it.Name], it] }.collect()
ch_fusionreport_ref                             = Channel.fromPath(params.fusionreport_ref).map                 { it -> [[id:it.Name], it] }.collect()
ch_intervals_gunzip                             = Channel.fromPath(params.intervals_gunzip).map                 { it -> [[id:it.Name], it] }.collect()
ch_known_indels_tbi                             = Channel.fromPath(params.known_indels_tbi).map                 { it -> [[id:it.Name], it] }.collect()
ch_cnvkit_antitarget                            = Channel.fromPath(params.cnvkit_antitarget).map                { it -> [[id:it.Name], it] }.collect()
ch_gatk_interval_list                           = Channel.fromPath(params.gatk_interval_list).map               { it -> [[id:it.Name], it] }.collect()
ch_arriba_ref_blocklist                         = Channel.fromPath(params.arriba_ref_blocklist).map             { it -> [[id:it.Name], it] }.collect()
ch_arriba_ref_cytobands                         = Channel.fromPath(params.arriba_ref_cytobands).map             { it -> [[id:it.Name], it] }.collect()
ch_intervals_gunzip_index                       = Channel.fromPath(params.intervals_gunzip_index).map           { it -> [[id:it.Name], it] }.collect()
ch_arriba_ref_known_fusions                     = Channel.fromPath(params.arriba_ref_known_fusions).map         { it -> [[id:it.Name], it] }.collect()
ch_arriba_ref_protein_domains                   = Channel.fromPath(params.arriba_ref_protein_domains).map       { it -> [[id:it.Name], it] }.collect()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                      IMPORT ASSETS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_dummy_file                                   = file("$projectDir/workflows/assets/multiqc/deseq2_pca_header.txt"         , checkIfExists: true)
ch_pca_header_multiqc                           = file("$projectDir/workflows/assets/multiqc/deseq2_pca_header.txt"         , checkIfExists: true)
ch_biotypes_header_multiqc                      = file("$projectDir/workflows/assets/multiqc/biotypes_header.txt"           , checkIfExists: true)
sample_status_header_multiqc                    = file("$projectDir/workflows/assets/multiqc/sample_status_header.txt"      , checkIfExists: true)
ch_clustering_header_multiqc                    = file("$projectDir/workflows/assets/multiqc/deseq2_clustering_header.txt"  , checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNAINN {

    take:
    ch_samplesheet

    main:

    ch_reports                                  = Channel.empty()
    ch_versions                                 = Channel.empty()
    ch_multiqc_files                            = Channel.empty()

    //
    // Create channel from input file provided through params.input
    //
    Channel
        .fromSamplesheet("input")
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map {
            validateInputSamplesheet(it)
        }
        .branch {
            meta, fastqs ->
                single  : fastqs.size() == 1
                    return [ meta, fastqs.flatten() ]
                multiple: fastqs.size() > 1
                    return [ meta, fastqs.flatten() ]
        }
        .set { ch_fastq }
    ch_fastq_single     = ch_fastq.single.map   {meta, fastqs -> addReadgroupToMeta(meta, fastqs)}
    ch_fastq_multiple   = ch_fastq.multiple.map {meta, fastqs -> addReadgroupToMeta(meta, fastqs)}

    // PREPROCESSING starts here...

    //
    // MODULE: Run FastQC
    //
    FASTQC(ch_fastq.mix())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (ch_fastq_multiple)
    .reads
    .mix(ch_fastq_single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    // DOWNSAMPLINGS starts here...

    if (params.run_downsamplings) {

        if (params.ds_totalreads_aim) {

            //
            // MODULE: Run Downsampling with seqtk
            //
            DOWNSAMPLING_SEQTK(ch_cat_fastq, params.ds_totalreads_aim)
            ch_versions = ch_versions.mix(DOWNSAMPLING_SEQTK.out.versions)
            ch_downsampled_reads = DOWNSAMPLING_SEQTK.out.downsampled_reads

            ch_fastqs = ch_downsampled_reads

        } else {

            //
            // MODULE: Run in-house script for counting reads
            //
            DOWNSAMPLING_COUNT(ch_cat_fastq)
            ch_versions = ch_versions.mix(DOWNSAMPLING_COUNT.out.versions)
            ch_global_min_reads = DOWNSAMPLING_COUNT.out.total_reads
                .map { file -> 
                    def count = file.text.trim()
                    count.toInteger()
                }
                .collect()
                .map { counts -> 
                    def min_count = counts.min()
                    min_count
                }

            //
            // MODULE: Run Downsampling with seqtk
            //
            DOWNSAMPLING_SEQTK(ch_cat_fastq, ch_global_min_reads)
            ch_versions = ch_versions.mix(DOWNSAMPLING_SEQTK.out.versions)
            ch_downsampled_reads = DOWNSAMPLING_SEQTK.out.downsampled_reads

            ch_fastqs = ch_downsampled_reads

        }

    } else {

        ch_fastqs = ch_cat_fastq

    }

    // UMIPROCESSING starts here...

    // Subworkflow Channels
    ch_umiprocessing_output = Channel.empty()

    if (params.run_umiprocessing) {
        
        //
        // SUBWORKFLOW: UMI processing
        //
        UMIPROCESSING(
            ch_fai,
            ch_gtf,
            ch_bwa2,
            ch_dict,
            ch_fasta,
            ch_fastqs,
            ch_refflat,
            ch_star_index,
            ch_rrna_intervals
        )
        ch_ubam                                     = UMIPROCESSING.out.ubam
        ch_raw_bam					                = UMIPROCESSING.out.raw_bam
        ch_bam_grouped				                = UMIPROCESSING.out.group_bam
        ch_reads_finalized                          = UMIPROCESSING.out.reads_finalized
        ch_versions                                 = ch_versions.mix(UMIPROCESSING.out.versions)
        ch_multiqc_files                            = ch_multiqc_files.mix(UMIPROCESSING.out.multiqc_files)

    } else {

        //
        // SUBWORKFLOW: Deduplication & Recalibration
        //
        DEDUPONLY4RNA(
            ch_fai,
            ch_gtf,
            ch_dict,
            ch_fasta,
            ch_fastqs,
            ch_star_index
        )
        ch_versions					                = ch_versions.mix(DEDUPONLY4RNA.out.versions)
        ch_split_reads                              = Channel.empty()
        ch_multiqc_files			                = ch_multiqc_files.mix(DEDUPONLY4RNA.out.multiqc_files)
        ch_bam_finalized			                = DEDUPONLY4RNA.out.bam_final
        ch_reads_finalized			                = DEDUPONLY4RNA.out.reads_final

    }

    if (params.run_fusion_splice) {

        //
        // SUBWORKFLOW: Run Fusion/Junction/Splicing discovery
        //
        FUSION_SPLICE(
            ch_fai,
            ch_gtf,
            ch_fasta,
            ch_chrgtf,
            ch_refflat,
            ch_versions,
            ch_hgnc_ref,
            ch_hgnc_date,
            ch_reads_finalized,
            ch_star_index,
            ch_samplesheet,
            ch_multiqc_files,
            ch_rrna_intervals,
            ch_fusionreport_ref,
            ch_arriba_ref_blocklist,
            ch_arriba_ref_known_fusions,
            ch_arriba_ref_protein_domains
        )
        ch_versions                                 = ch_versions.mix(FUSION_SPLICE.out.versions)
        ch_multiqc_files                            = ch_multiqc_files.mix(FUSION_SPLICE.out.multiqc_files)
    } else {
        ch_versions                                 = ch_versions
        ch_multiqc_files                            = ch_multiqc_files
    }

    if (params.run_genexpression) {

        //
        // SUBWORKFLOW: Run Gene Expression / Counts
        //
        GENEXPRESSION(
            ch_gtf,
            ch_fasta,
            ch_rsem_ref,
            ch_versions,
            ch_star_index,
            ch_samplesheet,
            ch_salmon_index,
            ch_multiqc_files,
            ch_kallisto_index,
            ch_reads_finalized,
            ch_pca_header_multiqc,
            ch_biotypes_header_multiqc,
            ch_clustering_header_multiqc
        )
        ch_versions                                 = ch_versions.mix(GENEXPRESSION.out.versions)
        ch_multiqc_files                            = ch_multiqc_files.mix(GENEXPRESSION.out.multiqc_files)
    } else {
        ch_versions                                 = ch_versions
        ch_multiqc_files                            = ch_multiqc_files
    }

    if (params.run_variantdscvry) {

        //
        // SUBWORKFLOW: Run Variant Calling
        //
        VARIANTDSCVRY(
            ch_bed,
            ch_fai,
            ch_gtf,
            ch_dict,
            ch_dbsnp,
            ch_fasta,
            ch_dbsnp_tbi,
            ch_reads_finalized,
            ch_star_index,
            ch_known_indels,
            ch_rrna_intervals,
            ch_known_indels_tbi,
            ch_gatk_interval_list
        )
        ch_versions                                 = ch_versions.mix(VARIANTDSCVRY.out.versions)
        ch_multiqc_files                            = ch_multiqc_files.mix(VARIANTDSCVRY.out.multiqc_files)
    } else {
        ch_versions                                 = ch_versions
        ch_multiqc_files                            = ch_multiqc_files
    }

//    if (params.run_copynumberalt) {
//
//        //
//        // SUBWORKFLOW: Run Copy Number Alterations
//        //
//        COPUNUMBERALT()
//        ch_versions                                 = ch_versions.mix(COPUNUMBERALT.out.versions)
//        ch_multiqc_files                            = ch_multiqc_files.mix(COPUNUMBERALT.out.multiqc_files)
//    } else {
//        ch_versions                                 = ch_versions
//        ch_multiqc_files                            = ch_multiqc_files
//    }
//
//    if (params.run_immunoncology) {
//
//        //
//        // SUBWORKFLOW: Run HLAinn
//        //
//        IMMUNONCOLOGY()
//        ch_versions                                 = ch_versions.mix(IMMUNONCOLOGY.out.versions)
//        ch_multiqc_files                            = ch_multiqc_files.mix(IMMUNONCOLOGY.out.multiqc_files)
//    } else {
//        ch_versions                                 = ch_versions
//        ch_multiqc_files                            = ch_multiqc_files
//    }
//
//    if (params.run_telomerefeats) {
//
//        //
//        // SUBWORKFLOW: Run TEMPUS
//        //
//        TELOMEREFEATS()
//        ch_versions                                 = ch_versions.mix(TELOMEREFEATS.out.versions)
//        ch_multiqc_files                            = ch_multiqc_files.mix(TELOMEREFEATS.out.multiqc_files)
//    } else {
//        ch_versions                                 = ch_versions
//        ch_multiqc_files                            = ch_multiqc_files
//    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                           = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config                    = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                             = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                              = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                         = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description       = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) :
                                                  file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                      = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                            = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                            = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                            = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml',sort: true))

    ch_multiqc_files = ch_multiqc_files
        .flatten() 
        .map { it instanceof Map || it instanceof List ? it[-1] : it } 
        .mix(
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false),
            ch_collated_versions
        )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report                              = MULTIQC.out.report.toList()                                   // channel: /path/to/multiqc_report.html
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                        FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Import libraries
//
import groovy.json.JsonSlurper

//
// Create MultiQC tsv custom content from a list of values
//
def multiqcTsvFromList(tsv_data, header) {
    def tsv_string = ""
    if (tsv_data.size() > 0) {
        tsv_string += "${header.join('\t')}\n"
        tsv_string += tsv_data.join('\n')
    }
    return tsv_string
}

//
// Add readgroup to meta and remove lane
//
def addReadgroupToMeta(meta, files) {
    def flowcell = flowcellLaneFromFastq(files[0])
    def CN = params.seq_center ? "CN:${params.seq_center}\\t" : ''
    def read_group = "\"@RG\\tID:${meta.id}\\t${CN}PU:${flowcell}\\tSM:${meta.id}\\tLB:${meta.id}\\tPL:${params.seq_platform}\""
    meta  = meta + [read_group: read_group.toString()]

    return [ meta, files ]
}

//
// Parse first line of a FASTQ file, return the flowcell id and lane number.
//
def flowcellLaneFromFastq(path) {
    def line
    path.withInputStream {
        InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        BufferedReader buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(':')
    String fcid

    if (fields.size() >= 7) {
        // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
        // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
        fcid = fields[2]
    } else if (fields.size() == 5) {
        fcid = fields[0]
    }
    return fcid
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                         THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
