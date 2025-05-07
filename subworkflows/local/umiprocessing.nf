/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                       UMIPROCESSING SUBWORKFLOW                                                    
*******************************************************************************************************************************************************************************************************

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                             IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMBLASTER                                                                                                                } from '../../modules/local/samblaster/main' // <- In use
include { STAR_FILBAM                                                                                                               } from '../../modules/local/align_fil_bam/star/main'
include { STAR_RAWBAM                                                                                                               } from '../../modules/local/align_raw_bam/star/main'
include { ALIGN_BAM_FIN                                                                                                             } from '../../modules/local/umi_align_bam/main' // <- In use
include { ALIGN_BAM_ORI                                                                                                             } from '../../modules/local/umi_align_bam/main' // <- In use
include { FASTQ_CONSENSUS                                                                                                           } from '../../modules/local/fastqc_consensus/main' // <- In use
include { FGBIO_FASTQTOBAM                                                                                                          } from '../../modules/nf-core/fgbio/fastqtobam/main'
include { FGBIO_SORTCONBAM                                                                                                          } from '../../modules/local/fgbio/sortconbam/main.nf' // <- In use
include { FGBIO_CORRECTUMIS                                                                                                         } from '../../modules/local/fgbio/correctumis/main' // <- New in use
include { SAMTOOLS_FASTQ_FIL                                                                                                        } from '../../modules/local/samtools/fastq_fil/main'
include { SAMTOOLS_FASTQ_RAW                                                                                                        } from '../../modules/local/samtools/fastq_raw/main'
include { COLLECTHSMETRICS_DUP                                                                                                      } from '../../modules/local/picard/collecthsmetrics/main' // <- In use
include { COLLECTHSMETRICS_FIN                                                                                                      } from '../../modules/local/picard/collecthsmetrics/main' // <- In use
include { COLLECTHSMETRICS_ORI                                                                                                      } from '../../modules/local/picard/collecthsmetrics/main' // <- In use
include { COLLECTHSMETRICS_SIM                                                                                                      } from '../../modules/local/picard/collecthsmetrics/main' // <- In use
include { SAMTOOLS_COLLATEFASTQ                                                                                                     } from '../../modules/local/samtools/collatefastq/main'
include { FGBIO_GROUPREADSBYUMI                                                                                                     } from '../../modules/local/fgbio/groupreadsbyumi/main'
include { SAMTOOLS_SORT_INDEX_FIN                                                                                                   } from '../../modules/local/samtools/sort_index/main' // <- In use
include { SAMTOOLS_SORT_INDEX_ORI                                                                                                   } from '../../modules/local/samtools/sort_index/main' // <- In use
include { FGBIO_FILTERCONSENSUSREADS                                                                                                } from '../../modules/local/fgbio/filterconsensusreads/main'
include { PICARD_COLLECTRNASEQMETRICS                                                                                               } from '../../modules/local/picard/collectrnaseqmetrics/main' // <- In Beta
include { FGBIO_COLLECTDUPLEXSEQMETRICS                                                                                             } from '../../modules/nf-core/fgbio/collectduplexseqmetrics/main'
include { PICARD_COLLECTMULTIPLEMETRICS                                                                                             } from '../../modules/local/picard/collectmultiplemetrics/main' // <- In use
include { FGBIO_CALLMOLECULARCONSENSUSREADS                                                                                         } from '../../modules/nf-core/fgbio/callmolecularconsensusreads/main'
include { FGBIO_ERRORRATEBYREADPOSITION_FIN                                                                                         } from '../../modules/local/fgbio/errorratebyreadposition/main' // <- In use
include { FGBIO_ERRORRATEBYREADPOSITION_ORI                                                                                         } from '../../modules/local/fgbio/errorratebyreadposition/main' // <- In use

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
    ch_star_index

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
    ALIGN_BAM_ORI(ch_ubam, ch_fasta, ch_fai, ch_dict, ch_star_index, ch_gtf, params.seq_platform, params.seq_center, sort)
    ch_versions = ch_versions.mix(ALIGN_BAM_ORI.out.versions.first())
    ch_raw_bam = ALIGN_BAM_ORI.out.bam

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
    // MODULE: Collect RNAseq metrics
    //
    PICARD_COLLECTRNASEQMETRICS(ch_bam_fcu_stix, params.rrna_refflat, params.rrna_intervals)
    ch_versions = ch_versions.mix(PICARD_COLLECTRNASEQMETRICS.out.versions)
    ch_rnaseq_metrics = Channel.empty().mix(PICARD_COLLECTRNASEQMETRICS.out.metrics)

    //
    // MODULE: Run ErrorRateByReadPosition 
    //
    FGBIO_ERRORRATEBYREADPOSITION_ORI(ch_bam_fcu_sort, ch_fasta, ch_fai, ch_dict, params.known_sites, params.known_sites_tbi, params.rrna_intervals)
    ch_versions = ch_versions.mix(FGBIO_ERRORRATEBYREADPOSITION_ORI.out.versions.first())

    //
    // MODULE: Run Picard's Collect HS Metrics for raw BAM files
    //
    COLLECTHSMETRICS_ORI(ch_bam_fcu_sort, ch_bam_fcu_indx, ch_fasta, ch_fai, ch_dict, params.hsmetrics_baits, params.hsmetrics_trgts)
    ch_versions = ch_versions.mix(COLLECTHSMETRICS_ORI.out.versions.first())
    ch_coverage_raw  = COLLECTHSMETRICS_ORI.out.coverage
    ch_hsmetrics_raw = COLLECTHSMETRICS_ORI.out.hsmetrics

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
    ch_bam_grouped = FGBIO_GROUPREADSBYUMI.out.bam

    //
    // MODULE: Run fgbio CollectDuplexSeqMetrics
    //
    FGBIO_COLLECTDUPLEXSEQMETRICS(ch_bam_grouped, params.rrna_intervals)
    ch_multiqc_files = ch_multiqc_files.mix(FGBIO_COLLECTDUPLEXSEQMETRICS.out.metrics.map{it[1]}.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FGBIO_COLLECTDUPLEXSEQMETRICS.out.pdf.map{it[1]}.collect())   
    ch_versions = ch_versions.mix(FGBIO_COLLECTDUPLEXSEQMETRICS.out.versions.first())

    //
    // MODULE: Run fgbio CallMolecularConsensusReads
    //
    FGBIO_CALLMOLECULARCONSENSUSREADS(ch_bam_grouped, params.call_min_reads, params.call_min_baseq)
    ch_versions = ch_versions.mix(FGBIO_CALLMOLECULARCONSENSUSREADS.out.versions.first())
    ch_consensus_bam = FGBIO_CALLMOLECULARCONSENSUSREADS.out.bam

    //
    // MODULE: Run fgbio SortBam
    //
    FGBIO_SORTCONBAM(ch_consensus_bam)
    ch_versions = ch_versions.mix(FGBIO_SORTCONBAM.out.versions.first())
    ch_consensus_bam_sorted = FGBIO_SORTCONBAM.out.bam

    //
    // MODULE: Run FgBIO FilterConsensusReads to produce the "Final", "Duplex Consensus" & "Simplex Consensus" BAM files
    //
    FGBIO_FILTERCONSENSUSREADS(ch_consensus_bam_sorted, params.fasta, params.fai, params.filter_min_reads, params.filter_min_base_quality, params.filter_max_base_error_rate, params.filter_max_read_error_rate, params.filter_max_no_call_fraction)
    ch_versions = ch_versions.mix(FGBIO_FILTERCONSENSUSREADS.out.versions.first())
    ch_bam_bai_final_fil = FGBIO_FILTERCONSENSUSREADS.out.suplex_bam_bai
    ch_bam_bai_duplex_fil = FGBIO_FILTERCONSENSUSREADS.out.duplex_bam_bai
    ch_bam_bai_simplex_fil = FGBIO_FILTERCONSENSUSREADS.out.simplex_bam_bai

//    //
//    // MODULE: Align with BWA mem
//    //
//    ALIGN_BAM_FIN(ch_bam_bai_final_fil, ch_bam_bai_simplex_fil, ch_bam_bai_duplex_fil, ch_star_index, ch_fai, ch_fasta, ch_dict, ch_gtf, params.star_seq_platform, params.star_seq_center)
//    ch_versions = ch_versions.mix(ALIGN_BAM_FIN.out.versions.first())
//    ch_bam_fin = ALIGN_BAM_FIN.out.bam
//    ch_bam_duplex = ALIGN_BAM_FIN.out.duplex_bam
//    ch_bam_simplex = ALIGN_BAM_FIN.out.simplex_bam

    //
    // MODULE: Align with BWA mem
    //
    ALIGN_BAM_FIN(ch_bam_bai_final_fil, ch_bam_bai_duplex_fil, ch_bam_bai_simplex_fil, ch_fasta, ch_fai, ch_dict, ch_bwa2)
    ch_versions = ch_versions.mix(ALIGN_BAM_FIN.out.versions.first())
    ch_bam_fin = ALIGN_BAM_FIN.out.bam
    ch_bam_duplex = ALIGN_BAM_FIN.out.duplex_bam
    ch_bam_simplex = ALIGN_BAM_FIN.out.simplex_bam

    //
    // MODULE: Run SamToools Sort & Index
    //
    SAMTOOLS_SORT_INDEX_FIN(ch_bam_fin, ch_fai, ch_fasta, ch_bam_duplex, ch_bam_simplex)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_INDEX_FIN.out.versions)
    ch_bam_fin_sort = SAMTOOLS_SORT_INDEX_FIN.out.bam
    ch_bam_fin_indx = SAMTOOLS_SORT_INDEX_FIN.out.bai
    ch_bam_fin_stix = SAMTOOLS_SORT_INDEX_FIN.out.bam_bai
    ch_bam_dup_stix = SAMTOOLS_SORT_INDEX_FIN.out.bam_duplex
    ch_bam_sim_stix = SAMTOOLS_SORT_INDEX_FIN.out.bam_simplex

    //
    // MODULE: Run ErrorRateByReadPosition in Final BAM
    //
    FGBIO_ERRORRATEBYREADPOSITION_FIN(ch_bam_fin_sort, ch_fasta, ch_fai, ch_dict, params.known_sites, params.known_sites_tbi, params.interval_list)
    ch_versions = ch_versions.mix(FGBIO_ERRORRATEBYREADPOSITION_FIN.out.versions)

    //
    // MODULE: Run Picard's Collect HS Metrics for consensus BAM files
    //
    COLLECTHSMETRICS_FIN(ch_bam_fin_stix, ch_fasta, ch_fai, ch_dict, params.hsmetrics_baits, params.hsmetrics_trgts)
    ch_versions = ch_versions.mix(COLLECTHSMETRICS_FIN.out.versions.first())
    ch_coverage_con  = COLLECTHSMETRICS_FIN.out.coverage
    ch_hsmetrics_con = COLLECTHSMETRICS_FIN.out.hsmetrics

    //
    // MODULE: Run Picard's Collect HS Metrics for consensus BAM files
    //
    COLLECTHSMETRICS_DUP(ch_bam_dup_stix, ch_fasta, ch_fai, ch_dict, params.hsmetrics_baits, params.hsmetrics_trgts)
    ch_versions = ch_versions.mix(COLLECTHSMETRICS_DUP.out.versions.first())
    ch_coverage_con  = COLLECTHSMETRICS_DUP.out.coverage
    ch_hsmetrics_con = COLLECTHSMETRICS_DUP.out.hsmetrics

    //
    // MODULE: Run Picard's Collect HS Metrics for consensus BAM files
    //
    COLLECTHSMETRICS_SIM(ch_bam_sim_stix, ch_fasta, ch_fai, ch_dict, params.hsmetrics_baits, params.hsmetrics_trgts)
    ch_versions = ch_versions.mix(COLLECTHSMETRICS_DUP.out.versions.first())
    ch_coverage_con  = COLLECTHSMETRICS_SIM.out.coverage
    ch_hsmetrics_con = COLLECTHSMETRICS_SIM.out.hsmetrics

    //
    // MODULE: Extract FastQ reads from BAM
    //
    SAMTOOLS_COLLATEFASTQ(ch_bam_fin_stix, ch_fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_COLLATEFASTQ.out.versions)
    ch_consensus_reads = SAMTOOLS_COLLATEFASTQ.out.fastq

    //
    // MODULE: Run FastQC
    //
    FASTQ_CONSENSUS(ch_consensus_reads)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_CONSENSUS.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQ_CONSENSUS.out.versions)

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

    emit:

    ubam            = ch_ubam
    raw_bam         = ch_raw_bam
    versions        = ch_collated_versions
    group_bam       = ch_bam_grouped
    multiqc_files   = ch_multiqc_files
    finalized_bam   = ch_bam_fin_stix
    reads_finalized = ch_consensus_reads

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                       THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/