/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                              FUSION_SPLICE SUBWORKFLOW                                                    
*****************************************************************************************************************************************************
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP                                                                         } from '../../modules/local/fastp/main'
include { ARRIBA                                                                        } from '../../modules/nf-core/arriba/main'
include { FUSVIZ                                                                        } from '../../modules/local/fusviz/main'
include { CAT_CAT                                                                       } from '../../modules/nf-core/cat/cat/main'
include { IGVREPORTS                                                                    } from '../../modules/nf-core/igvreports/main' 
include { STARFUSION                                                                    } from '../../modules/local/starfusion/detect/main'
include { VCF_COLLECT                                                                   } from '../../modules/local/vcf_collect/main'
include { STAR_ARRIBA                                                                   } from '../../modules/local/arriba/star/main'
include { STAR_FUSION                                                                   } from '../../modules/local/starfusion/star/main'
include { ARRIBA_INDEX                                                                  } from '../../modules/local/arriba/index/main'
include { FUSIONREPORT                                                                  } from '../../modules/local/fusionreport/detect/main'
include { FUSIONCATCHER                                                                 } from '../../modules/local/fusioncatcher/detect/main'
include { FUSIONINSPECTOR                                                               } from '../../modules/local/fusioninspector/main'
include { PORTCULLIS_FULL                                                               } from '../../modules/nf-core/portcullis/full/main'
include { STARFUSION_INDEX                                                              } from '../../modules/local/starfusion/index/main'
include { ARRIBA_VISUALISATION                                                          } from '../../modules/local/arriba/visualisation/main'
include { AGAT_CONVERTSPGFF2TSV                                                         } from '../../modules/nf-core/agat/convertspgff2tsv/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                 IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FUSIONREPORT_WORKFLOW                                                         } from '../local/fusionreport_workflow'
include { softwareVersionsToYAML                                                        } from '../nf-core/utils_nfcore_pipeline'
include { FUSIONINSPECTOR_WORKFLOW                                                      } from '../local/fusioninspector_workflow'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                   RUN SUBWORKFLOW 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FUSION_SPLICE {

    take:

    ch_fai
    ch_gtf
    ch_fasta
    ch_chrgtf
    ch_refflat
    ch_versions
    ch_hgnc_ref
    ch_hgnc_date
    ch_reads_all
    ch_star_index
    ch_samplesheet
    ch_multiqc_files
    ch_rrna_intervals
    ch_fusionreport_ref
    ch_arriba_ref_blocklist
    ch_arriba_ref_known_fusions
    ch_arriba_ref_protein_domains

    main:

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
    // Combine the BAM & Arriba's Output By ID
    //
    ch_arribaviz_input = ch_bam_star_arriba_indexed
        .join(ch_arriba_fusions)
        .map { id, meta_bam, bam, bai, meta_tsv, tsv ->
            tuple(
                meta_bam, 
                meta_bam,  bam, bai,
                meta_tsv, tsv
            )
        }

    //
    // MODULE: Run Arriba visualisation tool
    //
    ARRIBA_VISUALISATION (ch_arribaviz_input, params.gtf, params.arriba_ref_protein_domains, params.arriba_ref_cytobands)
    ch_arriba_visualisation = ARRIBA_VISUALISATION.out.pdf
    ch_versions = ch_versions.mix(ARRIBA_VISUALISATION.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(ch_arriba_visualisation.collect{it[1]}.ifEmpty([]))

    //
    // MODULE: Run Portcullis
    //
    PORTCULLIS_FULL (ch_bam_star_arriba, params.bed, params.fai, params.fasta)
    ch_versions = ch_versions.mix(PORTCULLIS_FULL.out.versions)
    ch_portcullis_log = PORTCULLIS_FULL.out.log
    ch_portcullis_bed = PORTCULLIS_FULL.out.pass_junctions_bed
    ch_portcullis_tab = PORTCULLIS_FULL.out.pass_junctions_tab

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
    // MODULE: Run FastP to filter by min lenght since FusionCatcher can handle much
    //
    FASTP(ch_reads_all)
    ch_versions = ch_versions.mix(FASTP.out.versions)
    ch_reads_4fusioncatcher = FASTP.out.reads

    //
    // WORKFLOW: Run FusionCatcher
    //
    FUSIONCATCHER(ch_reads_4fusioncatcher, params.fusioncatcher_ref)
    ch_versions = ch_versions.mix(FUSIONCATCHER.out.versions)
    ch_fusioncatcher_fusions = FUSIONCATCHER.out.fusions

    //
    // WORKFLOW: Run FusionReport
    //
    FUSIONREPORT_WORKFLOW (ch_reads_all, ch_fusionreport_ref, ch_arriba_fusions, ch_starfusion_fusions, ch_fusioncatcher_fusions)
    ch_versions = ch_versions.mix(FUSIONREPORT_WORKFLOW.out.versions)

    //
    // WORKFLOW: Run FusionInspector
    //
    FUSIONINSPECTOR_WORKFLOW (
        ch_reads_all,
        FUSIONREPORT_WORKFLOW.out.fusion_list,
        FUSIONREPORT_WORKFLOW.out.fusion_list_filtered,
        FUSIONREPORT_WORKFLOW.out.report,
        FUSIONREPORT_WORKFLOW.out.csv,
        ch_bam_star_fusion_indexed,
        ch_hgnc_ref,
        ch_hgnc_date
    )
    ch_versions = ch_versions.mix(FUSIONINSPECTOR_WORKFLOW.out.versions)
    ch_fusioninspectortsv = FUSIONINSPECTOR_WORKFLOW.out.fusioninspectortsv

//    //
//    // MODULE: Run IGV Reports
//    //
//    IGVREPORTS()
//    ch_versions = ch_versions.mix(IGVREPORTS.out.versions)
//

    //
    // Join annotated SVs with BAM pairs based on patient
    //
    ch_fusviz_input = ch_bam_star_fusion_indexed
        .join(ch_fusioninspectortsv)
        .map { id, meta_bam, bam, bai, meta_tsv, tsv ->
            tuple(
                meta_bam, 
                meta_bam,  bam, bai,
                meta_tsv, tsv
            )
        }

    //
    // MODULE: Run FusViz
    //
    FUSVIZ(ch_fusviz_input, params.annotations, params.cytobands, params.fusviz_chr, params.protein_domains)
    ch_fusviz_pdf = FUSVIZ.out.pdf
    ch_versions = ch_versions.mix(FUSVIZ.out.versions)

    emit:

    versions        = ch_versions
    multiqc_files   = ch_multiqc_files 

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                       THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/