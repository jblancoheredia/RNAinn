/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                IMPORT PLUGINS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap                                                      } from 'plugin/nf-validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { FASTP                                                                 } from '../modules/nf-core/fastp/main'
include { ARRIBA                                                                } from '../modules/nf-core/arriba/main'
include { FASTQC                                                                } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                                               } from '../modules/nf-core/multiqc/main'
include { CAT_FASTQ                                                             } from '../modules/nf-core/cat/fastq/main'
include { STAR_ALIGN                                                            } from '../modules/nf-core/star/align/main'
include { PORTCULLIS                                                            } from '../modules/local/portcullis/main'
include { STARFUSION                                                            } from '../modules/local/starfusion/detect/main'
include { FGBIO_SORT                                                            } from '../modules/local/fgbio/sorts/main'
include { STAR_ARRIBA                                                           } from '../modules/local/arriba/star/main'
include { STAR_FUSION                                                           } from '../modules/local/starfusion/star/main'
include { STAR_FILBAM                                                           } from '../modules/local/align_fil_bam/star/main'
include { STAR_RAWBAM                                                           } from '../modules/local/align_raw_bam/star/main'
include { ARRIBA_INDEX                                                          } from '../modules/local/arriba/index/main'
include { FUSIONREPORT                                                          } from '../modules/local/fusionreport/detect/main'
include { FUSIONCATCHER                                                         } from '../modules/local/fusioncatcher/detect/main'
include { KALLISTO_QUANT                                                        } from '../modules/nf-core/kallisto/quant/main'
include { FGBIO_ZIPFILBAM                                                       } from '../modules/local/fgbio/zipfilbam/main'
include { FGBIO_ZIPRAWBAM                                                       } from '../modules/local/fgbio/ziprawbam/main'
include { STARFUSION_INDEX                                                      } from '../modules/local/starfusion/index/main'
include { FGBIO_FASTQTOBAM                                                      } from '../modules/nf-core/fgbio/fastqtobam/main'
include { SAMTOOLS_FASTQ_FIL                                                    } from '../modules/local/samtools/fastq_fil/main'
include { SAMTOOLS_FASTQ_RAW                                                    } from '../modules/local/samtools/fastq_raw/main'
include { ARRIBA_VISUALISATION                                                  } from '../modules/local/arriba/visualisation/main'
include { SAMTOOLS_COLLATEFASTQ                                                 } from '../modules/nf-core/samtools/collatefastq/main'
include { FGBIO_GROUPREADSBYUMI                                                 } from '../modules/nf-core/fgbio/groupreadsbyumi/main'
include { FGBIO_FILTERCONSENSUSREADS                                            } from '../modules/nf-core/fgbio/filterconsensusreads/main'
include { PICARD_COLLECTRNASEQMETRICS                                           } from '../modules/local/picard/collectrnaseqmetrics/main'
include { PICARD_COLLECTINSERTSIZEMETRICS                                       } from '../modules/nf-core/picard/collectinsertsizemetrics/main'
include { FGBIO_CALLMOLECULARCONSENSUSREADS                                     } from '../modules/nf-core/fgbio/callmolecularconsensusreads/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                             IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CALL_VARIANTS                                                         } from '../subworkflows/local/call_variants'
include { GENE_EXPRESSN                                                         } from '../subworkflows/local/gene_expressn'
include { FUSION_SPLICE                                                         } from '../subworkflows/local/fusion_splice'
include { PREPROCESSING                                                         } from '../subworkflows/local/preprocessing'
include { paramsSummaryMultiqc                                                  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                                                } from '../subworkflows/local/utils_nfcore_rnainn_pipeline'
include { validateInputSamplesheet                                              } from '../subworkflows/local/utils_nfcore_rnainn_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                            CREATE CUSTOM CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_bed = Channel.fromPath(params.bed).map { it -> [[id:it.Name], it] }.collect()
ch_fai = Channel.fromPath(params.fai).map { it -> [[id:it.Name], it] }.collect()
ch_gtf = Channel.fromPath(params.gtf).map { it -> [[id:it.Name], it] }.collect()
ch_dict = Channel.fromPath(params.dict).map { it -> [[id:it.Name], it] }.collect()
ch_dbsnp = Channel.fromPath(params.dbsnp).map { it -> [[id:it.Name], it] }.collect()
ch_fasta = Channel.fromPath(params.fasta).map { it -> [[id:it.Name], it] }.collect()
ch_chrgtf = Channel.fromPath(params.chrgtf).map { it -> [[id:it.Name], it] }.collect()
ch_rsem_ref = Channel.fromPath(params.rsem_ref).map { it -> [[id:it.Name], it] }.collect()
ch_hgnc_ref = Channel.fromPath(params.hgnc_ref).map { it -> [[id:it.Name], it] }.collect()
ch_dbsnp_tbi = Channel.fromPath(params.dbsnp_tbi).map { it -> [[id:it.Name], it] }.collect()
ch_hgnc_date = Channel.fromPath(params.hgnc_date).map { it -> [[id:it.Name], it] }.collect()
ch_star_index = Channel.fromPath(params.star_index).map { it -> [[id:it.Name], it] }.collect()
ch_known_indels = Channel.fromPath(params.known_indels).map { it -> [[id:it.Name], it] }.collect()
ch_multiqc_logo = Channel.fromPath(params.multiqc_logo).map { it -> [[id:it.Name], it] }.collect()
ch_salmon_index = Channel.fromPath(params.salmon_index).map { it -> [[id:it.Name], it] }.collect()
ch_kallisto_index = Channel.fromPath(params.kallisto_index).map { it -> [[id:it.Name], it] }.collect()
ch_rrna_intervals  = Channel.fromPath(params.rrna_intervals).map { it -> [[id:it.Name], it] }.collect()
ch_fusionreport_ref = Channel.fromPath(params.fusionreport_ref).map { it -> [[id:it.Name], it] }.collect()
ch_known_indels_tbi = Channel.fromPath(params.known_indels_tbi).map { it -> [[id:it.Name], it] }.collect()
ch_gatk_interval_list = Channel.fromPath(params.gatk_interval_list).map { it -> [[id:it.Name], it] }.collect()
ch_arriba_ref_blocklist = Channel.fromPath(params.arriba_ref_blocklist).map { it -> [[id:it.Name], it] }.collect()
ch_arriba_ref_cytobands = Channel.fromPath(params.arriba_ref_cytobands).map { it -> [[id:it.Name], it] }.collect()
ch_arriba_ref_known_fusions = Channel.fromPath(params.arriba_ref_known_fusions).map { it -> [[id:it.Name], it] }.collect()
ch_arriba_ref_protein_domains = Channel.fromPath(params.arriba_ref_protein_domains).map { it -> [[id:it.Name], it] }.collect()

// Header files for MultiQC
ch_dummy_file                   = file("$projectDir/workflows/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
ch_pca_header_multiqc           = file("$projectDir/workflows/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
ch_biotypes_header_multiqc      = file("$projectDir/workflows/assets/multiqc/biotypes_header.txt", checkIfExists: true)
sample_status_header_multiqc    = file("$projectDir/workflows/assets/multiqc/sample_status_header.txt", checkIfExists: true)
ch_clustering_header_multiqc    = file("$projectDir/workflows/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                               IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                  RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNAINN {

    take:
    ch_samplesheet

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

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
    
    ch_fastq_single     = ch_fastq.single
    ch_fastq_multiple   = ch_fastq.multiple

    //
    // SUBWORKFLOW: Run UMI handling
    //
    PREPROCESSING(
        ch_gtf,
        ch_fasta,
        ch_star_index,
        ch_samplesheet,
        ch_fastq_single,
        ch_fastq_multiple
        
    )

    ch_raw_bam					= PREPROCESSING.out.raw_bam
    ch_reads_all                = PREPROCESSING.out.all_reads
    ch_cat_fastq				= PREPROCESSING.out.raw_reads
    ch_bam_grouped				= PREPROCESSING.out.group_bam
    ch_ubam_sorted              = PREPROCESSING.out.ubam_sorted
    ch_bam_consensus			= PREPROCESSING.out.consensus_bam
    ch_bam_consensus_filtered	= PREPROCESSING.out.cons_filt_bam
    ch_reads_consensus_filtered = PREPROCESSING.out.reads_cons_fil

    ch_versions                 = ch_versions.mix(PREPROCESSING.out.versions)
    ch_multiqc_files            = ch_multiqc_files.mix(PREPROCESSING.out.multiqc_files)

    //
    // SUBWORKFLOW: Run Fusion/Junction/Splicing discovery
    //
    FUSION_SPLICE(
        ch_gtf,
        ch_fasta,
        ch_chrgtf,
        ch_versions,
        ch_hgnc_ref,
        ch_hgnc_date,
        ch_reads_all,
        ch_star_index,
        ch_samplesheet,
        ch_multiqc_files,
        ch_fusionreport_ref,
        ch_arriba_ref_blocklist,
        ch_arriba_ref_known_fusions,
        ch_arriba_ref_protein_domains
    )

    ch_versions                 = ch_versions.mix(FUSION_SPLICE.out.versions)
    ch_multiqc_files            = ch_multiqc_files.mix(FUSION_SPLICE.out.multiqc_files)

    //
    // SUBWORKFLOW: Run Gene Expression / Counts
    //
    GENE_EXPRESSN(
        ch_gtf,
        ch_fasta,
        ch_rsem_ref,
        ch_versions,
        ch_reads_all,
        ch_star_index,
        ch_samplesheet,
        ch_salmon_index,
        ch_multiqc_files,
        ch_kallisto_index,
        ch_pca_header_multiqc,
        ch_biotypes_header_multiqc,
        ch_clustering_header_multiqc
    )

    ch_versions                 = ch_versions.mix(GENE_EXPRESSN.out.versions)
    ch_multiqc_files            = ch_multiqc_files.mix(GENE_EXPRESSN.out.multiqc_files)
    
    //
    // SUBWORKFLOW: Run Variant Calling
    //
    CALL_VARIANTS(
        ch_bed,
        ch_fai,
        ch_gtf,
        ch_dict,
        ch_dbsnp,
        ch_fasta,
        ch_dbsnp_tbi,
        ch_reads_all,
        ch_star_index,
        ch_known_indels,
        ch_rrna_intervals,
        ch_known_indels_tbi,
        ch_gatk_interval_list

    )
    ch_versions                 = ch_versions.mix(CALL_VARIANTS.out.versions)
    ch_multiqc_files            = ch_multiqc_files.mix(CALL_VARIANTS.out.multiqc_files)

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
    ch_multiqc_config                       = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config                = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                         = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                          = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                     = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description   = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                  = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                        = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                        = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                        = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml',sort: true))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                       THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
