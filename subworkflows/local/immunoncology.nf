/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                       IMMUNONCOLOGY SUBWORKFLOW                                                    
*******************************************************************************************************************************************************************************************************

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                             IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HLAHD                                                                                                                     } from '../../modules/local/hlahd/main'
include { HLAIMGT                                                                                                                   } from '../../modules/local/hlaimgt/main'
include { OPTITYPE                                                                                                                  } from '../../modules/local/optitype/main'
include { HLALA_TYPING                                                                                                              } from '../../modules/local/hlala/typing/main'  
include { ALIGN_HLA_CHR6                                                                                                            } from '../../modules/local/bwamem2/chr6/main'
include { ALIGN_HLA_IMGT                                                                                                            } from '../../modules/local/bwamem2/imgt/main'

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

workflow IMMUNONCOLOGY {

    take:
    ch_bam
    ch_bai
    ch_fastqs
    ch_fin_fastqs

    main:
    ch_versions = Channel.empty()
//    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run BWAMEM2 to align to IMGT-HLA reference
    //
    ALIGN_HLA_IMGT(ch_fastqs, params.bwa2_imgt, params.imgt_fa, params.imgt_fai)
    ch_versions = ch_versions.mix(ALIGN_HLA_IMGT.out.versions)
    ch_imgt_bam = ALIGN_HLA_IMGT.out.bam
    ch_imgt_bai = ALIGN_HLA_IMGT.out.bai

    //
    // MODULE: Run BWAMEM2 to align to IMGT-HLA reference
    //
    ALIGN_HLA_CHR6(ch_fastqs, params.bwa2_chr6, params.chr6_fa, params.chr6_fai)
    ch_versions = ch_versions.mix(ALIGN_HLA_CHR6.out.versions)
    ch_chr6_bam = ALIGN_HLA_CHR6.out.bam
    ch_chr6_bai = ALIGN_HLA_CHR6.out.bai

    //
    // MODULE: HLA-HD
    //
    HLAHD(ch_fastqs, params.read_length)
    ch_hlahd = HLAHD.out.hla_type
    ch_versions = ch_versions.mix(HLAHD.out.versions)

    //
    // MODULE: Optitype
    //
    OPTITYPE(ch_fastqs, params.seq_type, params.optitype_image)
    ch_versions = ch_versions.mix(OPTITYPE.out.versions)

    //
    // MODULE: Build IMGT table
    //
    HLAIMGT(ch_imgt_bam, ch_imgt_bai, params.bwa2_imgt)
    ch_imgt_tsv = HLAIMGT.out.tsv
    ch_versions = ch_versions.mix(HLAIMGT.out.versions)

//    //
//    // MODULE: Run HLA-LA 
//    //
//    HLALA_TYPING(ch_bam, ch_bai, params.hlala_graph)
//    ch_hlahd = HLALA_TYPING.out.results
//    ch_versions = ch_versions.mix(HLALA_TYPING.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'software_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    emit:

    versions        = ch_collated_versions
//    multiqc_files   = ch_multiqc_files 

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                                THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

