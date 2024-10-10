/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                              FUSION_SPLICE SUBWORKFLOW                                                    
*****************************************************************************************************************************************************
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/ 

include { ARRIBA                                                                    } from '../../modules/nf-core/arriba/main'
include { STARFUSION                                                                } from '../../modules/local/starfusion/detect/main'
include { PORTCULLIS                                                                } from '../../modules/local/portcullis/main'
include { STAR_ARRIBA                                                               } from '../../modules/local/arriba/star/main'
include { STAR_FUSION                                                               } from '../../modules/local/starfusion/star/main'
include { ARRIBA_INDEX                                                              } from '../../modules/local/arriba/index/main'
include { FUSIONREPORT                                                              } from '../../modules/local/fusionreport/detect/main'
include { FUSIONCATCHER                                                             } from '../../modules/local/fusioncatcher/detect/main'
include { STARFUSION_INDEX                                                          } from '../../modules/local/starfusion/index/main'
include { ARRIBA_VISUALISATION                                                      } from '../../modules/local/arriba/visualisation/main'


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

workflow RNAINN_FUSION_SPLICE {

    take:

    ch_gtf
    ch_fasta
    ch_chrgtf
    ch_versions
    ch_reads_all
    ch_star_index
    ch_samplesheet
    ch_multiqc_files
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
    ch_portcullis_log = PORTCULLIS.out.log
    ch_portcullis_bed = PORTCULLIS.out.bed
    ch_portcullis_tab = PORTCULLIS.out.tab
    ch_portcullis_txt = PORTCULLIS.out.txt

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
    // WORKFLOW: Run Fusion-Report
    //
    reads_fusions = ch_reads_all
        .join(ch_arriba_fusions, remainder: true)
        .join(ch_starfusion_fusions, remainder: true)
        .join(ch_fusioncatcher_fusions, remainder: true)

    FUSIONREPORT(reads_fusions, ch_fusionreport_ref, params.tools_cutoff)
    ch_fusionreport_list_filtered = FUSIONREPORT.out.fusion_list_filtered
    ch_versions = ch_versions.mix(FUSIONREPORT.out.versions)
    ch_fusionreport_list = FUSIONREPORT.out.fusion_list
    ch_fusionreport_report = FUSIONREPORT.out.report
    ch_fusionreport_csv = FUSIONREPORT.out.csv

    emit:

    versions        = ch_versions
    multiqc_files   = ch_multiqc_files 

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                       THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/