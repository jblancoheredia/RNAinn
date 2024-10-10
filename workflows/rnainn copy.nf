/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                       IMPORT PLUGINS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap                              } from 'plugin/nf-validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                       IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP                                         } from '../modules/nf-core/fastp/main'
include { ARRIBA                                        } from '../modules/nf-core/arriba/main'
include { FASTQC                                        } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                       } from '../modules/nf-core/multiqc/main'
include { CAT_FASTQ                                     } from '../modules/nf-core/cat/fastq/main'
include { STRINGTIE                                     } from '../modules/local/stringtie/main'
include { STAR_ALIGN                                    } from '../modules/nf-core/star/align/main'
include { PORTCULLIS                                    } from '../modules/local/portcullis/main'
include { STARFUSION                                    } from '../modules/local/starfusion/detect/main'
include { FGBIO_SORT                                    } from '../modules/local/fgbio/sorts/main'
include { STAR_ARRIBA                                   } from '../modules/local/arriba/star/main'
include { STAR_FUSION                                   } from '../modules/local/starfusion/star/main'
include { STAR_FILBAM                                   } from '../modules/local/align_fil_bam/star/main'
include { STAR_RAWBAM                                   } from '../modules/local/align_raw_bam/star/main'
include { ARRIBA_INDEX                                  } from '../modules/local/arriba/index/main'
include { FUSIONREPORT                                  } from '../modules/local/fusionreport/detect/main'
include { FUSIONCATCHER                                 } from '../modules/local/fusioncatcher/detect/main'
include { SALMON_TX2GENE                                } from '../modules/local/salmon/tx2gene/main'
include { KALLISTO_QUANT                                } from '../modules/nf-core/kallisto/quant/main'
include { SAMTOOLS_STATS                                } from '../modules/nf-core/samtools/stats/main'
include { STRINGTIE_MERGE                               } from '../modules/nf-core/stringtie/merge/main'
include { FGBIO_ZIPFILBAM                               } from '../modules/local/fgbio/zipfilbam/main'
include { FGBIO_ZIPRAWBAM                               } from '../modules/local/fgbio/ziprawbam/main'
include { SALMON_QUANT_FQS                              } from '../modules/local/salmon/quant_fqs/main'
include { STARFUSION_INDEX                              } from '../modules/local/starfusion/index/main'
include { FGBIO_FASTQTOBAM                              } from '../modules/nf-core/fgbio/fastqtobam/main'
include { SAMTOOLS_IDXSTATS                             } from '../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT                             } from '../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_FASTQ_FIL                            } from '../modules/local/samtools/fastq_fil/main'
include { SAMTOOLS_FASTQ_RAW                            } from '../modules/local/samtools/fastq_raw/main'
include { ARRIBA_VISUALISATION                          } from '../modules/local/arriba/visualisation/main'
include { SAMTOOLS_COLLATEFASTQ                         } from '../modules/nf-core/samtools/collatefastq/main'
include { FGBIO_GROUPREADSBYUMI                         } from '../modules/nf-core/fgbio/groupreadsbyumi/main'
include { FGBIO_FILTERCONSENSUSREADS                    } from '../modules/nf-core/fgbio/filterconsensusreads/main'
include { PICARD_COLLECTRNASEQMETRICS                   } from '../modules/local/picard/collectrnaseqmetrics/main'
include { PICARD_COLLECTINSERTSIZEMETRICS               } from '../modules/nf-core/picard/collectinsertsizemetrics/main'
include { FGBIO_CALLMOLECULARCONSENSUSREADS             } from '../modules/nf-core/fgbio/callmolecularconsensusreads/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                     IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNAINN_PREPROCESSING                          } from '../subworkflows/local/preprocessing'
include { paramsSummaryMultiqc                          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                        } from '../subworkflows/local/utils_nfcore_rnainn_pipeline'
include { validateInputSamplesheet                      } from '../subworkflows/local/utils_nfcore_rnainn_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                   CREATE CUSTOM CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_gtf = Channel.fromPath(params.gtf).map { it -> [[id:it.Name], it] }.collect()
ch_fasta = Channel.fromPath(params.fasta).map { it -> [[id:it.Name], it] }.collect()
ch_chrgtf = Channel.fromPath(params.chrgtf).map { it -> [[id:it.Name], it] }.collect()
ch_star_index = Channel.fromPath(params.star_index).map { it -> [[id:it.Name], it] }.collect()
ch_kallisto_index = Channel.fromPath(params.kallisto_index).map { it -> [[id:it.Name], it] }.collect()
ch_fusionreport_ref = Channel.fromPath(params.fusionreport_ref).map { it -> [[id:it.Name], it] }.collect()
ch_arriba_ref_blocklist = Channel.fromPath(params.arriba_ref_blocklist).map { it -> [[id:it.Name], it] }.collect()
ch_arriba_ref_cytobands = Channel.fromPath(params.arriba_ref_cytobands).map { it -> [[id:it.Name], it] }.collect()
ch_arriba_ref_known_fusions = Channel.fromPath(params.arriba_ref_known_fusions).map { it -> [[id:it.Name], it] }.collect()
ch_arriba_ref_protein_domains = Channel.fromPath(params.arriba_ref_protein_domains).map { it -> [[id:it.Name], it] }.collect()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                      RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    RNAINN_PREPROCESSING(
        ch_gtf,
        ch_fasta,
        ch_star_index,
        ch_samplesheet,
        ch_fastq_single,
        ch_fastq_multiple
        
    )

    ch_raw_bam					= RNAINN_PREPROCESSING.out.raw_bam
    ch_versions					= RNAINN_PREPROCESSING.out.versions
    ch_cat_fastq				= RNAINN_PREPROCESSING.out.raw_reads
    ch_bam_grouped				= RNAINN_PREPROCESSING.out.group_bam
    ch_ubam_sorted              = RNAINN_PREPROCESSING.out.ubam_sorted
    ch_multiqc_files			= RNAINN_PREPROCESSING.out.multiqc_files
    ch_bam_consensus			= RNAINN_PREPROCESSING.out.consensus_bam
    ch_bam_consensus_filtered	= RNAINN_PREPROCESSING.out.cons_filt_bam
    ch_reads_consensus_filtered = RNAINN_PREPROCESSING.out.reads_cons_fil

    //
    // MODULE: Run fastp
    //
    FASTP(ch_reads_consensus_filtered, params.adapter_fasta, false, false, false)
    ch_reads_all = FASTP.out.reads
    ch_fastp_html = FASTP.out.html
    ch_fastp_json = FASTP.out.json
    ch_versions = ch_versions.mix(FASTP.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(ch_fastp_html.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_fastp_json.collect{it[1]}.ifEmpty([]))

    //
    // MODULE: Run STAR for Arriba
    //
    STAR_ARRIBA(ch_reads_all, ch_star_index, ch_gtf, params.star_seq_platform, params.star_seq_center)
    ch_bam_star_arriba = STAR_ARRIBA.out.bam_sorted
    ch_versions = ch_versions.mix(STAR_ARRIBA.out.versions)
    ch_star_arriba_gene_count = STAR_ARRIBA.out.read_per_gene_tab
    ch_star_arriba_stats = STAR_ARRIBA.out.log_final
    ch_multiqc_files = ch_multiqc_files.mix(ch_star_arriba_stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_star_arriba_gene_count.collect{it[1]}.ifEmpty([]))

    //
    // MODULE: Run Arriba
    //
    ARRIBA(ch_bam_star_arriba, ch_fasta, ch_gtf, ch_arriba_ref_blocklist, ch_arriba_ref_known_fusions, [[],[]], [[],[]], ch_arriba_ref_protein_domains)
    ch_arriba_fusions = ARRIBA.out.fusions
    ch_arriba_fusion_fail = ARRIBA.out.fusions_fail.map{ meta, file -> return file}
    ch_versions = ch_versions.mix(ARRIBA.out.versions)

    //
    // MODULE: Index BAM file
    //
    ARRIBA_INDEX(ch_bam_star_arriba)
    ch_bai = ARRIBA_INDEX.out.bai
    ch_versions = ch_versions.mix(ARRIBA_INDEX.out.versions)
    ch_bam_star_arriba_indexed = STAR_ARRIBA.out.bam_sorted.join(ARRIBA_INDEX.out.bai)

    //
    // MODULE: Run Arriba visualisation tool
    //
    ARRIBA_VISUALISATION (ch_bam_star_arriba_indexed, ARRIBA.out.fusions, params.gtf, params.arriba_ref_protein_domains, params.arriba_ref_cytobands)
    ch_arriba_visualisation = ARRIBA_VISUALISATION.out.pdf
    ch_versions = ch_versions.mix(ARRIBA_VISUALISATION.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(ch_arriba_visualisation.collect{it[1]}.ifEmpty([]))

    //
    // MODULE: Run Portcullis
    //
    PORTCULLIS (ch_bam_star_arriba, params.bed, params.fasta)
    ch_versions = ch_versions.mix(PORTCULLIS.out.versions)

    //
    // WORKFLOW: Run STAR for STARfusion
    //
    STAR_FUSION(ch_reads_all, ch_star_index, ch_chrgtf, params.star_seq_platform, params.star_seq_center)
    ch_bam_star_fusion = STAR_FUSION.out.bam_sorted
    ch_star_fusion_stats = STAR_FUSION.out.log_final
    ch_versions = ch_versions.mix(STAR_FUSION.out.versions)
    ch_star_fusion_gene_count = STAR_FUSION.out.read_per_gene_tab
    ch_star_fusion_reads_junction = ch_reads_all.join(STAR_FUSION.out.junction)
    ch_multiqc_files = ch_multiqc_files.mix(ch_star_fusion_stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_star_fusion_gene_count.collect{it[1]}.ifEmpty([]))

    //
    // MODULE: Index BAM file
    //
    STARFUSION_INDEX(ch_bam_star_fusion)
    ch_versions = ch_versions.mix(STARFUSION_INDEX.out.versions)
    ch_bam_star_fusion_indexed = STAR_FUSION.out.bam_sorted.join(STARFUSION_INDEX.out.bai)

    //
    // MODULE: Run STARfusion
    //
    STARFUSION(ch_star_fusion_reads_junction, params.starfusion_ref)
    ch_versions = ch_versions.mix(STARFUSION.out.versions)
    ch_starfusion_fusions = STARFUSION.out.fusions

    //
    // WORKFLOW: Run FusionCatcher
    //
    FUSIONCATCHER(ch_reads_all, params.fusioncatcher_ref)
    ch_versions = ch_versions.mix(FUSIONCATCHER.out.versions)
    ch_fusioncatcher_fusions = FUSIONCATCHER.out.fusions

    //
    // WORKFLOW: Run StringTie
    //
    STRINGTIE(ch_bam_star_fusion_indexed, params.chrgtf)
    ch_versions = ch_versions.mix(STRINGTIE.out.versions)
    STRINGTIE
        .out
        .transcript_gtf
        .map { it -> it[1] }
        .set { stringtie_gtf }
    ch_versions = ch_versions.mix(STRINGTIE.out.versions)

    //
    // MODULE: Merge StringTie GTF files
    //
    STRINGTIE_MERGE (stringtie_gtf, params.chrgtf)
    ch_versions = ch_versions.mix(STRINGTIE_MERGE.out.versions)
    ch_stringtie_gtf = STRINGTIE_MERGE.out.gtf

    //
    // WORKFLOW: Run Fusion-Report
    //
    reads_fusions = ch_reads_all
        .join(ch_arriba_fusions, remainder: true)
        .join(ch_starfusion_fusions, remainder: true)
        .join(ch_fusioncatcher_fusions, remainder: true)

    FUSIONREPORT(reads_fusions, ch_fusionreport_ref, params.tools_cutoff)
    ch_fusion_list = FUSIONREPORT.out.fusion_list
    ch_fusion_list_filtered = FUSIONREPORT.out.fusion_list_filtered
    ch_versions = ch_versions.mix(FUSIONREPORT.out.versions)
    ch_report = FUSIONREPORT.out.report
    ch_csv = FUSIONREPORT.out.csv

    //
    // MODULE: Run Samtools Stat
    //
    SAMTOOLS_STATS(ch_bam_star_arriba_indexed, ch_fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)
    ch_samtools_stats = SAMTOOLS_STATS.out.stats

    //
    // MODULE: Run Samtools FlagStats
    //
    SAMTOOLS_FLAGSTAT(ch_bam_star_arriba_indexed)
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)
    ch_samtools_flagstat = SAMTOOLS_FLAGSTAT.out.flagstat

    //
    // MODULE: Run Samtools IdxStats
    //
    SAMTOOLS_IDXSTATS(ch_bam_star_arriba_indexed)
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions)
    ch_idxstats = SAMTOOLS_IDXSTATS.out.idxstats

    //
    // MODULE: Run STAR
    //
    STAR_ALIGN(ch_reads_all, ch_star_index, ch_gtf, params.star_seq_platform, params.star_seq_center)
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

    //
    // MODULE: Collect RNAseq metrics
    //
    PICARD_COLLECTRNASEQMETRICS(ch_bam_star_arriba_indexed, params.rrna_refflat, params.rrna_intervals)
    ch_versions = ch_versions.mix(PICARD_COLLECTRNASEQMETRICS.out.versions)
    ch_rnaseq_metrics = Channel.empty().mix(PICARD_COLLECTRNASEQMETRICS.out.metrics)

    //
    // MODULE: Collect Insert Size metrics
    //
    PICARD_COLLECTINSERTSIZEMETRICS(ch_bam_star_fusion_indexed)
    ch_versions = ch_versions.mix(PICARD_COLLECTINSERTSIZEMETRICS.out.versions)
    ch_insertsize_metrics = Channel.empty().mix(PICARD_COLLECTINSERTSIZEMETRICS.out.metrics)

    //
    // MODULE: Run Kallisto pseudoalignment
    //
    KALLISTO_QUANT (ch_reads_all, ch_kallisto_index, params.gtf, params.chromosomes, [], [])
    ch_pseudo_kallisto_results = KALLISTO_QUANT.out.results
    ch_pseudo_kallisto_multiqc = KALLISTO_QUANT.out.log
    ch_versions = ch_versions.mix(KALLISTO_QUANT.out.versions.first())

    //
    // MODULE: Run Salmon pseudoalignment
    //
    SALMON_QUANT_FQS(ch_reads_all, params.salmon_index, params.lib_type)
    ch_versions = ch_versions.mix(SALMON_QUANT_FQS.out.versions)
    ch_pseudo_salmon_results = SALMON_QUANT_FQS.out.results.collect{ it[1] }.map { [ [:], it ] }
    ch_pseudo_salmon_multiqc = ch_pseudo_salmon_results

    //
    // MODULE: Run Salmon TX2GENE
    //
    SALMON_TX2GENE(ch_gtf, ch_pseudo_salmon_results, 'salmon', 'gene_id', 'transcript_id')
    ch_versions = ch_versions.mix(SALMON_TX2GENE.out.versions)

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                      IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Import libraries
//
import groovy.json.JsonSlurper

//
// Function to determine library type by comparing type counts.
//
def calculateStrandedness(forwardFragments, reverseFragments, unstrandedFragments, stranded_threshold=0.8, unstranded_threshold=0.1) {
    def totalFragments = forwardFragments + reverseFragments + unstrandedFragments
    def totalStrandedFragments = forwardFragments + reverseFragments
    def library_strandedness = 'undetermined'
    if (totalStrandedFragments > 0) {
        def forwardProportion = forwardFragments / (totalStrandedFragments as double)
        def reverseProportion = reverseFragments / (totalStrandedFragments as double)
        def proportionDifference = Math.abs(forwardProportion - reverseProportion)
        if (forwardProportion >= stranded_threshold) {
            strandedness = 'forward'
        } else if (reverseProportion >= stranded_threshold) {
            strandedness = 'reverse'
        } else if (proportionDifference <= unstranded_threshold) {
            strandedness = 'unstranded'
        }
    }
    return [
        inferred_strandedness: strandedness,
        forwardFragments: (forwardFragments / (totalFragments as double)) * 100,
        reverseFragments: (reverseFragments / (totalFragments as double)) * 100,
        unstrandedFragments: (unstrandedFragments / (totalFragments as double)) * 100
    ]
}

//
// Function that parses Salmon quant 'lib_format_counts.json' output file to get inferred strandedness
//
def getSalmonInferredStrandedness(json_file, stranded_threshold = 0.8, unstranded_threshold = 0.1) {
    def libCounts = new JsonSlurper().parseText(json_file.text)
    def forwardKeys = ['SF', 'ISF', 'MSF', 'OSF']
    def reverseKeys = ['SR', 'ISR', 'MSR', 'OSR']
    def unstrandedKeys = ['IU', 'U', 'MU']
    def forwardFragments = forwardKeys.collect { libCounts[it] ?: 0 }.sum()
    def reverseFragments = reverseKeys.collect { libCounts[it] ?: 0 }.sum()
    def unstrandedFragments = unstrandedKeys.collect { libCounts[it] ?: 0 }.sum()
    return calculateStrandedness(forwardFragments, reverseFragments, unstrandedFragments, stranded_threshold, unstranded_threshold)
}

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                           THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
