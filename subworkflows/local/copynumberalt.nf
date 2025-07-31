/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                       COPYNUMBERALT SUBWORKFLOW                                                    
*******************************************************************************************************************************************************************************************************

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                             IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ONCOCNV                                                                                                                   } from '../../modules/local/oncocnv/main'
include { FACETS_CNV                                                                                                                } from '../../modules/local/facets/main'
include { CNVKIT_CALL                                                                                                               } from '../../modules/nf-core/cnvkit/call/main'
include { CNVKIT_BATCH                                                                                                              } from '../../modules/nf-core/cnvkit/batch/main'
include { CNVKIT_EXPORT                                                                                                             } from '../../modules/nf-core/cnvkit/export/main'
include { SEQUENZA_FITS                                                                                                             } from '../../modules/local/sequenza/fits/main'
include { SEQUENZA_SEQZ                                                                                                             } from '../../modules/local/sequenza/seqz/main'
include { BCFTOOLS_MPILEUP                                                                                                          } from '../../modules/local/bcftools/mpileup/main'
include { SAMTOOLS_MPILEUP                                                                                                          } from '../../modules/local/samtools/mpileup/main'
include { CNVKIT_REFERENCE                                                                                                          } from '../../modules/nf-core/cnvkit/reference/main'
include { CNVKIT_ANTITARGET                                                                                                         } from '../../modules/nf-core/cnvkit/antitarget/main'
include { CNVKIT_GENEMETRICS                                                                                                        } from '../../modules/nf-core/cnvkit/genemetrics/main'
include { CONTROLFREEC_OT_FREEC                                                                                                     } from '../../modules/local/controlfreec_ot/freec/main'
include { CONTROLFREEC_OT_FREEC2BED                                                                                                 } from '../../modules/local/controlfreec_ot/freec2bed/main' 
include { CONTROLFREEC_OT_FREEC2CIRCOS                                                                                              } from '../../modules/local/controlfreec_ot/freec2circos/main'
include { CONTROLFREEC_OT_ASSESSSIGNIFICANCE                                                                                        } from '../../modules/local/controlfreec_ot/assesssignificance/main'

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

workflow COPYNUMBERALT {

    take:
    ch_fai
    ch_fasta
    ch_raw_bam
    ch_raw_bai
    ch_intervals
    ch_targets_bed
    ch_consensus_bam
    ch_cnvkit_reference
    ch_cnvkit_antitarget

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

//    //
//    // MODULE: Run OncoCNV
//    //
//    ONCOCNV(ch_consensus_bam, ch_fasta, ch_targets_bed, ch_fai, params.normal_bai, params.normal_bam)
//    ch_versions = ch_versions.mix(ONCOCNV.out.versions)
//    ch_oncocnv_png      = ONCOCNV.out.png
//    ch_oncocnv_profile  = ONCOCNV.out.profile
//    ch_oncocnv_summary  = ONCOCNV.out.summary

    //
    // MODULE: Run BCFtools Mpileup
    //
    BCFTOOLS_MPILEUP(ch_consensus_bam, ch_fai, ch_fasta, ch_intervals)
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_MPILEUP.out.stats.collect{it[1]})
    ch_bcf_mpileup = BCFTOOLS_MPILEUP.out.vcf
    ch_bcf_mpileup_tbi = BCFTOOLS_MPILEUP.out.tbi

    //
    // MODULE: Run SamTools Mpileup
    //
    SAMTOOLS_MPILEUP(ch_consensus_bam, ch_intervals, params.fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_MPILEUP.out.versions)
    ch_sam_mpileup = SAMTOOLS_MPILEUP.out.pup
    ch_sam_mpileup_tbi = SAMTOOLS_MPILEUP.out.tbi

    //
    // MODULE: Run ControlFreec Freec 
    //
    CONTROLFREEC_OT_FREEC(ch_sam_mpileup, ch_consensus_bam, params.fasta, params.chr_fai, params.known_sites, params.known_sites_tbi, params.by_chr_dir, params.cf_mappability, params.intervals, params.cf_control_mpileup, params.cf_coeff, params.cf_contamination, params.cf_contamination_adjustment, params.cf_ploidy)
    ch_versions = ch_versions.mix(CONTROLFREEC_OT_FREEC.out.versions)
    ch_cfot_baf   = CONTROLFREEC_OT_FREEC.out.BAF
    ch_cfot_cnvs  = CONTROLFREEC_OT_FREEC.out.CNV
    ch_cfot_ratio = CONTROLFREEC_OT_FREEC.out.ratio

    //
    // MODULE: Run ControlFreec Assess Significance
    //
    CONTROLFREEC_OT_ASSESSSIGNIFICANCE(ch_cfot_cnvs, ch_cfot_ratio)
    ch_versions = ch_versions.mix(CONTROLFREEC_OT_ASSESSSIGNIFICANCE.out.versions)

    //
    // MODULE: Run ControlFreec Freec2BED
    //
    CONTROLFREEC_OT_FREEC2BED(ch_cfot_ratio)
    ch_versions = ch_versions.mix(CONTROLFREEC_OT_FREEC2BED.out.versions)

    //
    // MODULE: Run ControlFreec Freec2Circos
    //
    CONTROLFREEC_OT_FREEC2CIRCOS(ch_cfot_ratio)
    ch_versions = ch_versions.mix(CONTROLFREEC_OT_FREEC2CIRCOS.out.versions)

//    //
//    // MODULE: Run CNVkit AntiTarget Module just once per panel
//    //
//    CNVKIT_ANTITARGET(ch_targets_bed)
//    ch_versions = ch_versions.mix(CNVKIT_ANTITARGET.out.versions)
//    ch_antitargets = CNVKIT_ANTITARGET.out.bed

//    //
//    // MODULE: Run CNVkit Reference Module just once per panel
//    //
//    CNVKIT_REFERENCE(ch_fasta, ch_targets_bed, ch_cnvkit_antitarget)
//    ch_versions = ch_versions.mix(CNVKIT_REFERENCE.out.versions)
//    ch_cnvkit_reference = CNVKIT_REFERENCE.out.cnn

    //
    // MODULE: Run CNVKIT Batch
    //
    generate_pon = false
    CNVKIT_BATCH(ch_consensus_bam, ch_fasta, ch_fai, ch_cnvkit_antitarget, ch_cnvkit_reference, generate_pon)
    ch_versions = ch_versions.mix(CNVKIT_BATCH.out.versions.first())
    ch_cnvkit_call_input = CNVKIT_BATCH.out.cns.map{ meta, cns -> [meta, cns[2], []]}

    //
    // MODULE: Run CNVkit Call
    //
    CNVKIT_CALL(ch_cnvkit_call_input)
    ch_versions = ch_versions.mix(CNVKIT_CALL.out.versions.first())
    ch_cnv_calls_raw = CNVKIT_CALL.out.cns

    //
    // MODULE: Run CNVkit Export
    //
    CNVKIT_EXPORT(CNVKIT_CALL.out.cns)
    ch_versions = ch_versions.mix(CNVKIT_EXPORT.out.versions.first())
    ch_cnv_calls_export = CNVKIT_EXPORT.out.output

    //
    // MODULE: Run CNVkit GeneMetrics
    //
    ch_genemetrics_input = CNVKIT_BATCH.out.cnr.join(CNVKIT_BATCH.out.cns).map{ meta, cnr, cns -> [meta, cnr, cns[2]]}
    CNVKIT_GENEMETRICS(ch_genemetrics_input)
    ch_versions = ch_versions.mix(CNVKIT_GENEMETRICS.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(CNVKIT_GENEMETRICS.out.tsv)

    //
    // MODULE: Run FACETS 
    //
    FACETS_CNV(ch_consensus_bam, params.normal_bam, params.normal_bai, params.common_vcf, params.common_vcf_tbi)
    ch_versions = ch_versions.mix(FACETS_CNV.out.versions.first())

//    //
//    // MODULE: Run Sequenzautils BAM2seqz
//    //
//    SEQUENZA_SEQZ(ch_consensus_bam, ch_fasta, ch_fai, ch_sam_mpileup, ch_sam_mpileup_tbi, params.normal_mpileup, params.normal_mpileup_tbi, params.wigfile50, 50)
//    ch_versions = ch_versions.mix(SEQUENZA_SEQZ.out.versions.first())
//    ch_seqz = SEQUENZA_SEQZ.out.seqz
//
//    //
//    // MODULE: Run Sequenzautils BAM2seqz
//    //
//    SEQUENZA_FITS(ch_seqz)
//    ch_versions = ch_versions.mix(SEQUENZA_FITS.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'software_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    emit:

    versions            = ch_collated_versions
    sam_mpileup         = ch_sam_mpileup
    bcf_mpileup         = ch_bcf_mpileup
    multiqc_files       = ch_multiqc_files

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                                THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
