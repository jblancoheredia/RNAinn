/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                             VARIANT CALLING SUBWORKFLOW                                                    
*****************************************************************************************************************************************************
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/ 

include { STAR_ALIGNV                                                               } from '../../modules/local/star/alignv/main'
include { TABIX_TABIX                                                               } from '../../modules/nf-core/tabix/tabix/main' 
include { SNPEFF_SNPEFF                                                             } from '../../modules/nf-core/snpeff/snpeff/main'    
include { SAMTOOLS_SORT                                                             } from '../../modules/nf-core/samtools/sort/main'
include { ENSEMBLVEP_VEP                                                            } from '../../modules/nf-core/ensemblvep/vep/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_VC1                                      } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_VC2                                      } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_VC3                                      } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_VC1                                      } from '../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_VC2                                      } from '../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_VC3                                      } from '../../modules/nf-core/samtools/stats/main'
include { GATK4_APPLYBQSR                                                           } from '../../modules/nf-core/gatk4/applybqsr/main' 
include { GATK4_MERGEVCFS                                                           } from '../../modules/nf-core/gatk4/mergevcfs/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_VC1                                  } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_VC2                                  } from '../../modules/nf-core/tabix/bgziptabix/main'
include { SAMTOOLS_IDXSTATS as SAMTOOLS_IDXSTATS_VC1                                } from '../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_IDXSTATS as SAMTOOLS_IDXSTATS_VC2                                } from '../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_VC1                                } from '../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_VC2                                } from '../../modules/nf-core/samtools/flagstat/main'
include { GATK4_MARKDUPLICATES                                                      } from '../../modules/nf-core/gatk4/markduplicates/main'
include { GATK4_HAPLOTYPECALLER                                                     } from '../../modules/nf-core/gatk4/haplotypecaller/main'
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
ch_vep_cache                            = params.vep_cache         ?                Channel.fromPath(params.vep_cache).collect()     : []
ch_snpeff_db                            = params.snpeff_db         ?:               Channel.empty()
ch_vep_genome                           = params.vep_genome        ?:               Channel.empty()
ch_vep_species                          = params.vep_species       ?:               Channel.empty()
ch_snpeff_cache                         = params.snpeff_cache      ?                Channel.fromPath(params.snpeff_cache).collect()  : []
ch_multiqc_files                        =                                           Channel.empty()
ch_vep_cache_version                    = params.vep_cache_version ?:               Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                   RUN SUBWORKFLOW 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CALL_VARIANTS {

    take:

    ch_bed
    ch_fai
    ch_gft
    ch_dict
    ch_dbsnp
    ch_fasta
    ch_dbsnp_tbi
    ch_reads_all
    ch_star_index
    ch_known_indels
    ch_rrna_intervals
    ch_known_indels_tbi
    ch_gatk_interval_list

    main:

    //
    // MODULE: Run STAR
    //
    STAR_ALIGNV(ch_reads_all, ch_star_index, ch_gft, params.star_seq_platform, params.star_seq_center)
    ch_versions = ch_versions.mix(STAR_ALIGNV.out.versions)

    //
    // MODULE: Run SamTools Sort
    //
    SAMTOOLS_SORT(STAR_ALIGNV.out.bam, ch_fasta)
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

    ch_bam_variant_calling = Channel.empty()
    ch_bqsr_table   = Channel.empty()
    ch_interval_list_recalib = ch_gatk_interval_list.map{ meta, bed -> [bed] }.flatten()
    ch_splitncigar_bam_bai.combine(ch_interval_list_recalib)
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

    ch_bam_applybqsr       = ch_splitncigar_bam_bai.join(ch_bqsr_table, by: [0])
    ch_bam_recalibrated_qc = Channel.empty()

    ch_interval_list_applybqsr = ch_gatk_interval_list.map{ meta, bed -> [bed] }.flatten()
    ch_bam_applybqsr.combine(ch_interval_list_applybqsr)
        .map{ meta, bam, bai, table, interval -> [ meta, bam, bai, table, interval]
    }.set{ch_applybqsr_bam_bai_interval}

    //
    // MODULE: ApplyBaseRecalibrator from GATK4
    //
    GATK4_APPLYBQSR(ch_applybqsr_bam_bai_interval, params.fasta, params.fai, params.dict)
    bam_recalibrated = GATK4_APPLYBQSR.out.bam
    ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions)
//
//    //
//    // MODULE: Run SamTools Index
//    //
//    SAMTOOLS_INDEX_VC3(bam_recalibrated)
//    ch_bam_variant_calling = bam_recalibrated
//        .join(SAMTOOLS_INDEX_VC4.out.bai, by: [0], remainder: true)
//        .join(SAMTOOLS_INDEX_VC4.out.csi, by: [0], remainder: true)
//        .map{meta, bam, bai, csi ->
//            if (bai) [meta, bam, bai]
//            else [meta, bam, csi]
//        }
//    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_VC4.out.versions)
//
//    //
//    // MODULE: Run SamTools Stats
//    //
//    SAMTOOLS_STATS_VC3(ch_bam_variant_calling, ch_fasta)
//    ch_versions = ch_versions.mix(SAMTOOLS_STATS_VC3.out.versions)
//    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_VC3.out.stats.map{it[1]}.collect())
//
//    // Run haplotyper even in the absence of dbSNP files
//    if (!params.dbsnp){
//        ch_dbsnp = []
//        ch_dbsnp_tbi = []
//    }
//
//    ch_haplotypecaller_vcf = Channel.empty()
//    ch_haplotypecaller_interval_bam = ch_bam_variant_calling.combine(ch_gatk_interval_list)
//        .map{ meta, bam, bai, interval_list ->
//            new_meta = meta.clone()
//            new_meta.id = meta.id + "_" + interval_list.baseName
//            new_meta.sample = meta.id
//            [new_meta, bam, bai, interval_list]
//        }
//
//    //
//    // MODULE: HaplotypeCaller from GATK4 (Calls germline SNPs and indels via local re-assembly of haplotypes.)
//    // 
//    GATK4_HAPLOTYPECALLER(ch_haplotypecaller_interval_bam, ch_fasta, ch_fai, ch_dict, ch_dbsnp, ch_dbsnp_tbi)
//    ch_haplotypecaller_raw = GATK4_HAPLOTYPECALLER.out.vcf
//        .map{ meta, vcf ->
//            meta.id = meta.sample
//            [meta, vcf]}
//        .groupTuple()
//    ch_versions  = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.ifEmpty(null))
//
//    //
//    // MODULE: MergeVCFS from GATK4 (Merge multiple VCF files into one VCF)
//    //
//    GATK4_MERGEVCFS(ch_haplotypecaller_raw, ch_dict)
//    ch_haplotypecaller_vcf = GATK4_MERGEVCFS.out.vcf
//    ch_versions  = ch_versions.mix(GATK4_MERGEVCFS.out.versions.ifEmpty(null))
//
//    //
//    // MODULE: Index the VCF using TABIX
//    //
//    TABIX_TABIX(ch_haplotypecaller_vcf)
//    ch_haplotypecaller_vcf_tbi = ch_haplotypecaller_vcf
//        .join(TABIX_TABIX.out.tbi, by: [0], remainder: true)
//        .join(TABIX_TABIX.out.csi, by: [0], remainder: true)
//        .map{meta, vcf, tbi, csi ->
//            if (tbi) [meta, vcf, tbi]
//            else [meta, vcf, csi]
//        }
//    ch_versions     = ch_versions.mix(TABIX_TABIX.out.versions.ifEmpty(null))
//    ch_final_vcf    = ch_haplotypecaller_vcf
//
//    //
//    // MODULE: VariantFiltration from GATK4 (Filter variant calls based on certain criteria.)
//    // 
//    GATK4_VARIANTFILTRATION(ch_haplotypecaller_vcf_tbi, ch_fasta, ch_fai, ch_dict)
//    ch_filtered_vcf = GATK4_VARIANTFILTRATION.out.vcf
//    ch_final_vcf    = ch_filtered_vcf
//    ch_versions     = ch_versions.mix(GATK4_VARIANTFILTRATION.out.versions.ifEmpty(null))

//    //
//    // MODULE: Run SNPEFF
//    //
//    SNPEFF_SNPEFF(ch_final_vcf, snpeff_db, snpeff_cache)
//    ch_versions = ch_versions.mix(SNPEFF.out.versions)
//    ch_multiqc_files = ch_multiqc_files.mix(SNPEFF_SNPEFF.out.summary_html.map{it[1]}.collect())
//
//    //
//    // MODULES: Run Tabix in bgzip mode
//    //
//    TABIX_BGZIPTABIX_VC1(SNPEFF_SNPEFF.out.vcf)
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
//    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX_VC2.out.versions)

    emit:

    versions        = ch_versions
    multiqc_files   = ch_multiqc_files
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                       THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/