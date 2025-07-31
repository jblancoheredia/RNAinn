/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                             FUSION_INSPECTOR SUBWORKFLOW                                                    
*****************************************************************************************************************************************************
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CAT_CAT                                                                     } from '../../modules/nf-core/cat/cat/main'
include { VCF_COLLECT                                                                 } from '../../modules/local/vcf_collect/main'
include { FUSIONINSPECTOR                                                             } from '../../modules/local/fusioninspector/main'
include { AGAT_CONVERTSPGFF2TSV                                                       } from '../../modules/nf-core/agat/convertspgff2tsv/main'
include { FUSIONINSPECTOR_VISUALISATION                                               } from '../../modules/local/fusioninspector/visualisation/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                   RUN SUBWORKFLOW 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FUSIONINSPECTOR_WORKFLOW {
    take:
        reads
        fusion_list
        fusion_list_filtered
        fusionreport_out
        fusionreport_csv
        bam_sorted_indexed
        ch_hgnc_ref
        ch_hgnc_date

    main:
        ch_versions = Channel.empty()
        ch_fusioninspector_visualisation = Channel.empty()
        index ="${params.starfusion_ref}"

        ch_fusion_list = ( params.tools_cutoff > 1 ? fusion_list_filtered : fusion_list )
        .branch{
            no_fusions: it[1].size() == 0
            fusions: it[1].size() > 0
        }

        if (params.allowlist)  {
            ch_allowlist = ch_fusion_list.fusions.combine(Channel.value(file(params.allowlist, checkIfExists:true)))
                            .map { meta, fusions, allowlist -> [ meta, [fusions, allowlist] ] }

            CAT_CAT(ch_allowlist)
            ch_versions = ch_versions.mix(CAT_CAT.out.versions)
            ch_reads_fusion = reads.join(CAT_CAT.out.file_out )
        }
        else {
            ch_reads_fusion = reads.join(ch_fusion_list.fusions )
        }

        FUSIONINSPECTOR( ch_reads_fusion, index)
        ch_versions = ch_versions.mix(FUSIONINSPECTOR.out.versions)

        AGAT_CONVERTSPGFF2TSV(FUSIONINSPECTOR.out.out_gtf)
        ch_versions = ch_versions.mix(AGAT_CONVERTSPGFF2TSV.out.versions)

        fusion_data = FUSIONINSPECTOR.out.tsv_coding_effect.join(AGAT_CONVERTSPGFF2TSV.out.tsv).join(fusionreport_out).join(fusionreport_csv)
        VCF_COLLECT(fusion_data, ch_hgnc_ref, ch_hgnc_date)
        ch_versions = ch_versions.mix(VCF_COLLECT.out.versions)

        FUSIONINSPECTOR_VISUALISATION(bam_sorted_indexed, FUSIONINSPECTOR.out.tsv, params.gtf, params.arriba_ref_protein_domains, params.arriba_ref_cytobands)
        ch_versions = ch_versions.mix(FUSIONINSPECTOR_VISUALISATION.out.versions)
        ch_fusioninspector_visualisation = FUSIONINSPECTOR_VISUALISATION.out.pdf

    emit:
        ch_fusioninspector_visualisation
        versions                = ch_versions
        fusioninspectortsv      = FUSIONINSPECTOR.out.tsv
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                       THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/