/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                       UMIPROCESSING SUBWORKFLOW                                                    
*******************************************************************************************************************************************************************************************************

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                             IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMBLASTER                                                                                                                } from '../../modules/local/samblaster/main'
include { STAR_FILBAM                                                                                                               } from '../../modules/local/align_fil_bam/star/main'
include { STAR_RAWBAM                                                                                                               } from '../../modules/local/align_raw_bam/star/main'
include { ALIGN_BAM_CON                                                                                                             } from '../../modules/local/umi_align_bam/main'
include { ALIGN_BAM_RAW                                                                                                             } from '../../modules/local/umi_align_bam/main'
include { PRESEQ_CCURVE                                                                                                             } from '../../modules/local/preseq/ccurve/main'
include { PRESEQ_LCEXTRAP                                                                                                           } from '../../modules/local/preseq/lcextrap/main'
include { UMI_READ_COUNTS                                                                                                           } from '../../modules/local/umi_read_counts/main'
include { FASTQC_CONSENSUS                                                                                                          } from '../../modules/local/fastqc_consensus/main'
include { FGBIO_FASTQTOBAM                                                                                                          } from '../../modules/nf-core/fgbio/fastqtobam/main'
include { FGBIO_SORTCONBAM                                                                                                          } from '../../modules/local/fgbio/sortconbam/main.nf'
include { FGBIO_CORRECTUMIS                                                                                                         } from '../../modules/local/fgbio/correctumis/main'
include { SAMTOOLS_FASTQ_FIL                                                                                                        } from '../../modules/local/samtools/fastq_fil/main'
include { SAMTOOLS_FASTQ_RAW                                                                                                        } from '../../modules/local/samtools/fastq_raw/main'
include { COLLECT_UMI_METRICS                                                                                                       } from '../../modules/local/collect_umi_metrics/main'
include { COLLECTHSMETRICS_CON                                                                                                      } from '../../modules/local/picard/collecthsmetrics/main'
include { COLLECTHSMETRICS_RAW                                                                                                      } from '../../modules/local/picard/collecthsmetrics/main'
include { SAMTOOLS_COLLATEFASTQ                                                                                                     } from '../../modules/local/samtools/collatefastq/main'
include { FGBIO_GROUPREADSBYUMI                                                                                                     } from '../../modules/local/fgbio/groupreadsbyumi/main'
include { SAMTOOLS_SORT_INDEX_FIN                                                                                                   } from '../../modules/local/samtools/sort_index/main'
include { SAMTOOLS_SORT_INDEX_ORI                                                                                                   } from '../../modules/local/samtools/sort_index/main'
include { FGBIO_FILTERCONSENSUSREADS                                                                                                } from '../../modules/local/fgbio/filterconsensusreads/main'
include { FGBIO_COLLECTDUPLEXSEQMETRICS                                                                                             } from '../../modules/local/fgbio/collectduplexseqmetrics/main'
include { PICARD_COLLECTMULTIPLEMETRICS                                                                                             } from '../../modules/local/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTRNASEQMETRICS_CON                                                                                           } from '../../modules/local/picard/collectrnaseqmetrics/main'
include { PICARD_COLLECTRNASEQMETRICS_RAW                                                                                           } from '../../modules/local/picard/collectrnaseqmetrics/main'
include { FGBIO_CALLMOLECULARCONSENSUSREADS                                                                                         } from '../../modules/nf-core/fgbio/callmolecularconsensusreads/main'
include { FGBIO_ERRORRATEBYREADPOSITION_CON                                                                                         } from '../../modules/local/fgbio/errorratebyreadposition/main'
include { FGBIO_ERRORRATEBYREADPOSITION_RAW                                                                                         } from '../../modules/local/fgbio/errorratebyreadposition/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                          IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML                                                                                                    } from '../nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                            RUN SUBWORKFLOW 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow UMIPROCESSING {

    take:
    ch_fai
    ch_gtf
    ch_bwa2
    ch_dict
    ch_fasta
    ch_fastqs
    ch_refflat
    ch_star_index
    ch_rrna_intervals

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run fgbio FastqToBam
    //
    FGBIO_FASTQTOBAM(ch_fastqs)
    ch_versions = ch_versions.mix(FGBIO_FASTQTOBAM.out.versions.first())
    ch_ubam = FGBIO_FASTQTOBAM.out.bam

    //
    // MODULE: Align with bwa mem but avoid sort
    //
    sort = false
    ALIGN_BAM_RAW(ch_ubam, ch_fasta, ch_fai, ch_dict, ch_star_index, ch_gtf, params.seq_platform, params.seq_center, sort)
    ch_versions = ch_versions.mix(ALIGN_BAM_RAW.out.versions.first())
    ch_raw_bam = ALIGN_BAM_RAW.out.bam

    //
    // MODULE: Run fgbio correctumis
    //
    FGBIO_CORRECTUMIS(ch_raw_bam, params.correct_max_mismatch, params.correct_min_distance, params.correct_min_corrected)
    ch_multiqc_files = ch_multiqc_files.mix(FGBIO_CORRECTUMIS.out.metrics)
    ch_versions = ch_versions.mix(FGBIO_CORRECTUMIS.out.versions.first())
    ch_bam_fcu = FGBIO_CORRECTUMIS.out.bam

    //
    // MODULE: Run SamToools Sort & Index
    //
    SAMTOOLS_SORT_INDEX_ORI(ch_bam_fcu, ch_fasta, params.fai)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_INDEX_ORI.out.versions.first())
    ch_bam_fcu_sort = SAMTOOLS_SORT_INDEX_ORI.out.bam
    ch_bam_fcu_indx = SAMTOOLS_SORT_INDEX_ORI.out.bai
    ch_bam_fcu_stix = SAMTOOLS_SORT_INDEX_ORI.out.bam_bai

    //
    // MODULE: Run Picard Tool CollectMultipleMetrics
    //
    PICARD_COLLECTMULTIPLEMETRICS(ch_bam_fcu_stix, ch_fasta, ch_fai)
    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())

    //
    // MODULE: Run Picard's Collect RNAseq Metrics for raw BAM files
    //
    PICARD_COLLECTRNASEQMETRICS_RAW(ch_bam_fcu_sort, ch_bam_fcu_indx, ch_refflat, ch_rrna_intervals)
    ch_versions = ch_versions.mix(PICARD_COLLECTRNASEQMETRICS_RAW.out.versions.first())

    //
    // MODULE: Run Picard's Collect HS Metrics for raw BAM files
    //
    COLLECTHSMETRICS_RAW(ch_bam_fcu_stix, ch_fasta, ch_fai, ch_dict, params.optional_bait_intervals, params.optional_target_intervals, params.seq_library)
    ch_versions = ch_versions.mix(COLLECTHSMETRICS_RAW.out.versions.first())

    //
    // MODULE: Run ErrorRateByReadPosition 
    //
    FGBIO_ERRORRATEBYREADPOSITION_RAW(ch_bam_fcu_sort, ch_fasta, ch_fai, ch_dict, params.known_sites, params.known_sites_tbi, params.rrna_intervals)
    ch_versions = ch_versions.mix(FGBIO_ERRORRATEBYREADPOSITION_RAW.out.versions.first())

    //
    // MODULE: Run SamBlaster
    //
    SAMBLASTER(ch_bam_fcu)
    ch_versions = ch_versions.mix(SAMBLASTER.out.versions.first())

    //
    // MODULE: Run fgbio GroupReadsByUmi
    //
    FGBIO_GROUPREADSBYUMI(SAMBLASTER.out.bam, params.group_strategy, params.group_edits, params.group_include_secondary, params.group_allow_inter_contig, params.group_include_supplementary, params.group_min_map_q, params.group_include_non_pf_reads, params.group_mark_duplicates)
    ch_multiqc_files = ch_multiqc_files.mix(FGBIO_GROUPREADSBYUMI.out.histogram.map{it[1]}.collect())
    ch_versions = ch_versions.mix(FGBIO_GROUPREADSBYUMI.out.versions.first())
    ch_grouped_family_sizes = FGBIO_GROUPREADSBYUMI.out.histogram
    ch_bam_grouped = FGBIO_GROUPREADSBYUMI.out.bam

    //
    // MODULE: Run fgbio CollectDuplexSeqMetrics
    //
    FGBIO_COLLECTDUPLEXSEQMETRICS(ch_bam_grouped)
    ch_multiqc_files = ch_multiqc_files.mix(FGBIO_COLLECTDUPLEXSEQMETRICS.out.metrics.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FGBIO_COLLECTDUPLEXSEQMETRICS.out.pdf.map{it[1]}.collect())   
    ch_versions = ch_versions.mix(FGBIO_COLLECTDUPLEXSEQMETRICS.out.versions.first())

    //
    // MODULE: Run fgbio CallMolecularConsensusReads
    //
    FGBIO_CALLMOLECULARCONSENSUSREADS(ch_bam_grouped, params.call_min_reads, params.call_min_baseq)
    ch_versions = ch_versions.mix(FGBIO_CALLMOLECULARCONSENSUSREADS.out.versions.first())
    ch_called_bam = FGBIO_CALLMOLECULARCONSENSUSREADS.out.bam

    //
    // MODULE: Run fgbio SortBam
    //
    FGBIO_SORTCONBAM(ch_called_bam)
    ch_versions = ch_versions.mix(FGBIO_SORTCONBAM.out.versions.first())
    ch_called_bam_sorted = FGBIO_SORTCONBAM.out.bam

    //
    // MODULE: Run FgBIO FilterConsensusReads to produce the "Final", "Duplex Consensus" & "Simplex Consensus" BAM files
    //
    FGBIO_FILTERCONSENSUSREADS(ch_called_bam_sorted, params.fasta, params.fai, params.filter_min_reads, params.filter_min_base_quality, params.filter_max_base_error_rate, params.filter_max_read_error_rate, params.filter_max_no_call_fraction)
    ch_versions = ch_versions.mix(FGBIO_FILTERCONSENSUSREADS.out.versions.first())
    ch_filtered_bam = FGBIO_FILTERCONSENSUSREADS.out.bam_bai

    //
    // MODULE: Align with STAR, Zipp with FGbio and Sort & Index with Samtools
    //
    ALIGN_BAM_CON(ch_filtered_bam, ch_fasta, ch_fai, ch_dict, ch_star_index, ch_gtf, params.seq_platform, params.seq_center, sort)
    ch_versions = ch_versions.mix(ALIGN_BAM_CON.out.versions.first())
    ch_bam_consensus = ALIGN_BAM_CON.out.bam_bai

    //
    // MODULE: Run SamTools View to count reads accross the BAM files
    //
    UMI_READ_COUNTS(ch_ubam, ch_bam_fcu, ch_bam_grouped, ch_called_bam, ch_filtered_bam, ch_bam_consensus)
    ch_versions = ch_versions.mix(UMI_READ_COUNTS.out.versions.first())

    //
    // MODULE: Run SamTools View to count reads accross the BAM files
    //
    COLLECT_UMI_METRICS(ch_bam_consensus)
    ch_versions = ch_versions.mix(COLLECT_UMI_METRICS.out.versions.first())
    ch_cons_family_sizes = COLLECT_UMI_METRICS.out.cons_family_sizes

    //
    // MODULE: Run Preseq CCurve
    //
    PRESEQ_CCURVE(ch_cons_family_sizes)
    ch_versions = ch_versions.mix(PRESEQ_CCURVE.out.versions.first())

    //
    // MODULE: Run Preseq LCExtrap
    //
    PRESEQ_LCEXTRAP(ch_grouped_family_sizes)
    ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())

    //
    // MODULE: Run ErrorRateByReadPosition in Final BAM
    //
    FGBIO_ERRORRATEBYREADPOSITION_CON(ch_bam_consensus, ch_fasta, ch_fai, ch_dict, params.known_sites, params.known_sites_tbi, params.interval_list)
    ch_versions = ch_versions.mix(FGBIO_ERRORRATEBYREADPOSITION_CON.out.versions)

    //
    // MODULE: Run Picard's Collect RNAseq Metrics for consensus BAM files
    //
    PICARD_COLLECTRNASEQMETRICS_CON(ch_bam_consensus, ch_refflat, ch_rrna_intervals)
    ch_versions = ch_versions.mix(PICARD_COLLECTRNASEQMETRICS_CON.out.versions.first())

    //
    // MODULE: Run Picard's Collect HS Metrics for consensus BAM files
    //
    COLLECTHSMETRICS_CON(ch_bam_consensus, ch_fasta, ch_fai, ch_dict, params.optional_bait_intervals, params.optional_target_intervals, params.seq_library)
    ch_versions = ch_versions.mix(COLLECTHSMETRICS_CON.out.versions.first())

    //
    // MODULE: Extract FastQ reads from BAM
    //
    SAMTOOLS_COLLATEFASTQ(ch_bam_consensus, ch_fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_COLLATEFASTQ.out.versions)
    ch_consensus_reads = SAMTOOLS_COLLATEFASTQ.out.fastq

    //
    // MODULE: Run FastQC
    //
    FASTQC_CONSENSUS(ch_consensus_reads)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_CONSENSUS.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC_CONSENSUS.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_cversions }

    emit:

    ubam            = ch_ubam
    raw_bam         = ch_raw_bam
    versions        = ch_cversions
    group_bam       = ch_bam_grouped
    multiqc_files   = ch_multiqc_files
    finalized_bam   = ch_bam_consensus
    reads_finalized = ch_consensus_reads

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                       THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/