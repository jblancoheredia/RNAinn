/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                             VARIANT CALLING SUBWORKFLOW                                                    
*****************************************************************************************************************************************************
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/ 

include { TABIX_TABIX                                                               } from '../../modules/nf-core/tabix/tabix/main' 
include { SNPEFF_SNPEFF                                                             } from '../../modules/local/snpeff/snpeff/main'    
include { SAMTOOLS_SORT                                                             } from '../../modules/nf-core/samtools/sort/main'
include { STAR_ALIGNV                                                               } from '../../modules/local/star/align/alignv/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_VC1                                      } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_VC2                                      } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_VC1                                      } from '../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_VC2                                      } from '../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_VC3                                      } from '../../modules/nf-core/samtools/stats/main'
include { GATK4_APPLYBQSR                                                           } from '../../modules/local/gatk4/applybqsr/main' 
include { GATK4_MERGEVCFS                                                           } from '../../modules/nf-core/gatk4/mergevcfs/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_VC1                                  } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_VC2                                  } from '../../modules/nf-core/tabix/bgziptabix/main'
include { SAMTOOLS_IDXSTATS as SAMTOOLS_IDXSTATS_VC1                                } from '../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_IDXSTATS as SAMTOOLS_IDXSTATS_VC2                                } from '../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_VC1                                } from '../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_VC2                                } from '../../modules/nf-core/samtools/flagstat/main'
include { GATK4_MARKDUPLICATES                                                      } from '../../modules/nf-core/gatk4/markduplicates/main'
include { GATK4_HAPLOTYPECALLER                                                     } from '../../modules/local/gatk4/haplotypecaller/main'
include { GATK4_SPLITNCIGARREADS                                                    } from '../../modules/local/gatk4/splitncigarreads/main'
include { GATK4_BASERECALIBRATOR                                                    } from '../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_VARIANTFILTRATION                                                   } from '../../modules/nf-core/gatk4/variantfiltration/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                 IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML                                                    } from '../nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                 INITIALIZE CHANNELS 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_reports                              =                                           Channel.empty()
ch_versions                             =                                           Channel.empty()
ch_snpeff_db                            = params.snpeff_db         ?:               Channel.empty()
ch_snpeff_cache                         = params.snpeff_cache      ?                Channel.fromPath(params.snpeff_cache).collect()  : []
ch_multiqc_files                        =                                           Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                   RUN SUBWORKFLOW 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANTDSCVRY {

    take:

    ch_bed
    ch_fai
    ch_gft
    ch_dict
    ch_dbsnp
    ch_fasta
    ch_dbsnp_tbi
    ch_star_index
    ch_known_indels
    ch_rrna_intervals
    ch_hsmetrics_trgts
    ch_reads_finalized
    ch_known_indels_tbi
    ch_gatk_interval_list

    main:

    //
    // MODULE: Run STAR
    //
    STAR_ALIGNV(ch_reads_finalized, ch_star_index, ch_gft, params.star_seq_platform, params.star_seq_center)
    ch_versions = ch_versions.mix(STAR_ALIGNV.out.versions)
    ch_star_bam = STAR_ALIGNV.out.bam

    //
    // MODULE: Run SamTools Sort
    //
    SAMTOOLS_SORT(ch_star_bam, ch_fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    //
    // MODULE: Run SamTools Index
    //
    SAMTOOLS_INDEX_VC1( SAMTOOLS_SORT.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_VC1.out.versions)

    SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_INDEX_VC1.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX_VC1.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }
        .set { ch_bam_bai }

    //
    // MODULE: Run SamTools Stats
    //
    SAMTOOLS_STATS_VC1(ch_bam_bai, ch_fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_STATS_VC1.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_VC1.out.stats.map{it[1]}.collect())

    //
    // MODULE: Run SamTools FlagStat
    //
    SAMTOOLS_FLAGSTAT_VC1(ch_bam_bai)
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT_VC1.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_FLAGSTAT_VC1.out.flagstat.map{it[1]}.collect())

    //
    // MODULE: Run SamTools IdxStats
    //
    SAMTOOLS_IDXSTATS_VC1(ch_bam_bai)
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS_VC1.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_IDXSTATS_VC1.out.idxstats.map{it[1]}.collect())

    //
    // SUBWORKFLOW: Mark duplicates with GATK4
    //
    GATK4_MARKDUPLICATES(ch_bam_bai, params.fasta, params.fai)
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(GATK4_MARKDUPLICATES.out.metrics.map{it[1]}.collect())

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    SAMTOOLS_INDEX_VC2(GATK4_MARKDUPLICATES.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_VC2.out.versions)
    ch_genome_bam_bai = GATK4_MARKDUPLICATES.out.bam.join(SAMTOOLS_INDEX_VC2.out.bai)

    //
    // MODULE: Run SamTools Stats
    //
    SAMTOOLS_STATS_VC2(ch_genome_bam_bai, ch_fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_STATS_VC2.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_VC2.out.stats.map{it[1]}.collect())

    //
    // MODULE: Run SamTools FlagStat
    //
    SAMTOOLS_FLAGSTAT_VC2(ch_genome_bam_bai)
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT_VC2.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_FLAGSTAT_VC2.out.flagstat.map{it[1]}.collect())

    //
    // MODULE: Run SamTools IdxStats
    //
    SAMTOOLS_IDXSTATS_VC2(ch_genome_bam_bai)
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS_VC2.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_IDXSTATS_VC2.out.idxstats.map{it[1]}.collect())

    //
    // Module: Run GATK4 SplitNCigarReads
    //
    GATK4_SPLITNCIGARREADS(ch_genome_bam_bai, params.gatk_interval_list, ch_fasta, ch_fai, ch_dict)
    ch_splitncigar_bam_bai = GATK4_SPLITNCIGARREADS.out.bam
    ch_versions = ch_versions.mix(GATK4_SPLITNCIGARREADS.out.versions)

    //
    // MODULE: BaseRecalibrator from GATK4 -> Generates a recalibration table based on various co-variates
    // 

    ch_interval_list_flat = ch_gatk_interval_list.map{ meta, bed -> [bed] }.flatten()
    ch_splitncigar_bam_bai.combine(ch_interval_list_flat)
        .map{ meta, bam, bai, interval -> [ meta, bam, bai, interval]
    }.set{ch_splitncigar_bam_bai_interval}

    GATK4_BASERECALIBRATOR(
        ch_splitncigar_bam_bai_interval,
        params.fasta,
        params.fai,
        params.dict,
        ch_dbsnp,
        ch_dbsnp_tbi
    )
    ch_bqsr_table   = GATK4_BASERECALIBRATOR.out.table

    ch_reports  = ch_reports.mix(ch_bqsr_table.map{ meta, table -> table})
    ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions.first().ifEmpty(null))

    ch_splitncigar_bam_bai_bqsr_table_interval_combined = ch_splitncigar_bam_bai_interval.join(ch_bqsr_table)
    ch_applybqsr_bam_bai_interval = ch_splitncigar_bam_bai_bqsr_table_interval_combined.map { meta, bam, bai, interval, bqsr_table ->
        [meta, bam, bai, bqsr_table, interval]
    }

//    ch_bam_recalibrated_qc = Channel.empty()

    //
    // MODULE: ApplyBaseRecalibrator from GATK4
    //
    GATK4_APPLYBQSR(ch_applybqsr_bam_bai_interval, params.fasta, params.fai, params.dict)
    ch_bam_recalibrated = GATK4_APPLYBQSR.out.bam
    ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions)

    //
    // MODULE: Run SamTools Stats
    //
    SAMTOOLS_STATS_VC3(ch_bam_recalibrated, ch_fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_STATS_VC3.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_VC3.out.stats.map{it[1]}.collect())

    //
    // MODULE: HaplotypeCaller from GATK4 (Calls germline SNPs and indels via local re-assembly of haplotypes.)
    //
    // Run haplotyper even in the absence of dbSNP files
    if (!params.dbsnp){
        ch_dbsnp = []
        ch_dbsnp_tbi = []
    }

    GATK4_HAPLOTYPECALLER(ch_fai, ch_dict, ch_dbsnp, ch_fasta, ch_dbsnp_tbi, ch_hsmetrics_trgts, ch_bam_recalibrated)
    ch_versions  = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.ifEmpty(null))
    ch_haplotypecaller_raw = GATK4_HAPLOTYPECALLER.out.vcf
    ch_haplotypecaller_tbi = GATK4_HAPLOTYPECALLER.out.tbi
    ch_haplotypecaller_raw_tbi_combined = ch_haplotypecaller_raw.join(ch_haplotypecaller_tbi)
    ch_haplotypecaller_vcf_tbi = ch_haplotypecaller_raw_tbi_combined.map { meta, vcf, tbi ->
        [meta, vcf, tbi]
    }

    //
    // MODULE: VariantFiltration from GATK4 (Filter variant calls based on certain criteria.)
    // 
    GATK4_VARIANTFILTRATION(ch_haplotypecaller_vcf_tbi, ch_fasta, ch_fai, ch_dict)
    ch_versions  = ch_versions.mix(GATK4_VARIANTFILTRATION.out.versions.ifEmpty(null))
    ch_final_vcf = GATK4_VARIANTFILTRATION.out.vcf

//    //
//    // MODULE: Run SNPEFF
//    //
//    SNPEFF_SNPEFF(ch_final_vcf, params.snpeff_db, params.snpeff_cache)
//    ch_snpeff_vcf_annotated = SNPEFF_SNPEFF.out.vcf
//    ch_versions = ch_versions.mix(SNPEFF_SNPEFF.out.versions)
//    ch_multiqc_files = ch_multiqc_files.mix(SNPEFF_SNPEFF.out.summary_html.map{it[1]}.collect())
//
//    //
//    // MODULES: Run Tabix in bgzip mode
//    //
//    TABIX_BGZIPTABIX_VC1(ch_snpeff_vcf_annotated)
//    ch_snpeff_gz_tbi = TABIX_BGZIPTABIX_VC1.out.gz_tbi
//    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX_VC1.out.versions)
//
//    //
//    // MODULE: Run anotator EnsemblVEP
//    //
//    ENSEMBLVEP_VEP(ch_final_vcf, vep_genome, vep_species, vep_cache_version, vep_cache)
//    ch_versions = ch_versions.mix(ENSEMBLVEP.out.versions)
//    ch_multiqc_files = ch_multiqc_files.mix(ENSEMBLVEP_VEP.out.report.map{it[1]}.collect())
//
//    //
//    // MODULES: Run Tabix in bgzip mode
//    //
//    TABIX_BGZIPTABIX_VC2(ENSEMBLVEP_VEP.out.vcf)
//    ch_vep_gz_tbi = TABIX_BGZIPTABIX_VC1.out.gz_tbi
//    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX_VC2.out.versions)
//
//    //
//    // MODULE: Run GATK MergeVCFs
//    //
//    ch_ssnpeff_vep_combined = ch_snpeff_gz_tbi.join(ch_vep_gz_tbi)
//    ch_vcf_2_merge = ch_ssnpeff_vep_combined.map { meta, snpeff_vcf, vep_vcf ->
//        [meta, snpeff_vcf, vep_vcf]
//    }
//    GATK4_MERGEVCFS(ch_vcf_2_merge)
//    ch_versions = ch_versions.mix(GATK4_MERGEVCFS.out.versions)

    emit:

    versions        = ch_versions
    multiqc_files   = ch_multiqc_files
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                       THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/