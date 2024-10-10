/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                        COLLECT STATS AND FQ FILES SUBWORKFLOW                                                    
*****************************************************************************************************************************************************
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/ 

include { STAR_ALIGN                                                                } from '../modules/nf-core/star/align/main'
include { QUALIMAP_RNASEQ } from '../modules/nf-core/qualimap/rnaseq/main' 

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                 IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML                                                    } from '../nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                   RUN SUBWORKFLOW 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNAINN_COLL_STATS_FQ {

    take:

    ch_versions,
    ch_multiqc_files

    main:

//
//
//
//    //
//    // MODULE: Collect RNAseq metrics
//    //
//    PICARD_COLLECTRNASEQMETRICS(ch_bam_star_arriba_indexed, params.rrna_refflat, params.rrna_intervals)
//    ch_versions = ch_versions.mix(PICARD_COLLECTRNASEQMETRICS.out.versions)
//    ch_rnaseq_metrics = Channel.empty().mix(PICARD_COLLECTRNASEQMETRICS.out.metrics)
//
//    //
//    // MODULE: Collect Insert Size metrics
//    //
//    PICARD_COLLECTINSERTSIZEMETRICS(ch_bam_star_fusion_indexed)
//    ch_versions = ch_versions.mix(PICARD_COLLECTINSERTSIZEMETRICS.out.versions)
//    ch_insertsize_metrics = Channel.empty().mix(PICARD_COLLECTINSERTSIZEMETRICS.out.metrics)
//
//    //
//    // MODULE: Run Samtools Stat
//    //
//    SAMTOOLS_STATS(ch_bam_star_arriba_indexed, ch_fasta)
//    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)
//    ch_samtools_stats = SAMTOOLS_STATS.out.stats
//
//    //
//    // MODULE: Run Samtools FlagStats
//    //
//    SAMTOOLS_FLAGSTAT(ch_bam_star_arriba_indexed)
//    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)
//    ch_samtools_flagstat = SAMTOOLS_FLAGSTAT.out.flagstat
//
//    //
//    // MODULE: Run Samtools IdxStats
//    //
//    SAMTOOLS_IDXSTATS(ch_bam_star_arriba_indexed)
//    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions)
//    ch_idxstats = SAMTOOLS_IDXSTATS.out.idxstats



    emit:

    multiqc_files   = ch_multiqc_files.transpose().map{it[1]}
    versions        = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                       THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/