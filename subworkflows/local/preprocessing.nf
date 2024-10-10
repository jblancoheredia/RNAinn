/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                        PREPROCESSING-UMI-HANDLING SUBWORKFLOW                                                    
*****************************************************************************************************************************************************
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP                                                                 } from '../../modules/nf-core/fastp/main'
include { FASTQC                                                                } from '../../modules/nf-core/fastqc/main'
include { MULTIQC                                                               } from '../../modules/nf-core/multiqc/main'
include { CAT_FASTQ                                                             } from '../../modules/nf-core/cat/fastq/main'
include { FGBIO_SORT                                                            } from '../../modules/local/fgbio/sorts/main'
include { STAR_FILBAM                                                           } from '../../modules/local/align_fil_bam/star/main'
include { STAR_RAWBAM                                                           } from '../../modules/local/align_raw_bam/star/main'
include { FGBIO_ZIPFILBAM                                                       } from '../../modules/local/fgbio/zipfilbam/main'
include { FGBIO_ZIPRAWBAM                                                       } from '../../modules/local/fgbio/ziprawbam/main'
include { FGBIO_FASTQTOBAM                                                      } from '../../modules/nf-core/fgbio/fastqtobam/main'
include { SAMTOOLS_FASTQ_FIL                                                    } from '../../modules/local/samtools/fastq_fil/main'
include { SAMTOOLS_FASTQ_RAW                                                    } from '../../modules/local/samtools/fastq_raw/main'
include { SAMTOOLS_COLLATEFASTQ                                                 } from '../../modules/nf-core/samtools/collatefastq/main'
include { FGBIO_GROUPREADSBYUMI                                                 } from '../../modules/nf-core/fgbio/groupreadsbyumi/main'
include { FGBIO_FILTERCONSENSUSREADS                                            } from '../../modules/nf-core/fgbio/filterconsensusreads/main'
include { FGBIO_COLLECTDUPLEXSEQMETRICS                                         } from '../../modules/nf-core/fgbio/collectduplexseqmetrics/main'
include { FGBIO_CALLMOLECULARCONSENSUSREADS                                     } from '../../modules/nf-core/fgbio/callmolecularconsensusreads/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                 IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML                                                } from '../nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                   RUN SUBWORKFLOW 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNAINN_PREPROCESSING {

    take:
    ch_gtf
    ch_fasta
    ch_star_index
    ch_samplesheet
    ch_fastq_single
    ch_fastq_multiple

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    
    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (ch_fastq_multiple)
    .reads
    .mix(ch_fastq_single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC (ch_cat_fastq)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    //
    // MODULE: Run fgbio FastqToBam
    //
    FGBIO_FASTQTOBAM(ch_cat_fastq)
    ch_versions = ch_versions.mix(FGBIO_FASTQTOBAM.out.versions.first())

    //
    // MODULE: Run fgbio SortBam
    //
    FGBIO_SORT(FGBIO_FASTQTOBAM.out.bam, 'unmapped_sorted.bam')
    ch_versions = ch_versions.mix(FGBIO_SORT.out.versions.first())
    ch_ubam_sorted = FGBIO_SORT.out.ubam

    //
    // MODULE: Run SamTools fastq
    //
    SAMTOOLS_FASTQ_RAW(ch_ubam_sorted)
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ_RAW.out.versions.first())

    //
    // MODULE: Run STAR on Raw Data
    //
    STAR_RAWBAM(SAMTOOLS_FASTQ_RAW.out.fastq, ch_star_index, ch_gtf, params.star_seq_platform, params.star_seq_center)
    ch_versions = ch_versions.mix(STAR_RAWBAM.out.versions.first())
    ch_raw_bam = STAR_RAWBAM.out.bam

    //
    // MODULE: Run fgbio ZipperBAMs
    //
    FGBIO_ZIPRAWBAM(ch_ubam_sorted, ch_raw_bam, ch_fasta, params.dict)
    ch_versions = ch_versions.mix(FGBIO_ZIPRAWBAM.out.versions.first())

    //
    // MODULE: Run fgbio GroupReadsByUMI
    //
    FGBIO_GROUPREADSBYUMI(FGBIO_ZIPRAWBAM.out.bam, params.groupreadsbyumi_strategy)
    ch_multiqc_files = ch_multiqc_files.mix(FGBIO_GROUPREADSBYUMI.out.histogram.collect{it[1]})
    ch_versions = ch_versions.mix(FGBIO_GROUPREADSBYUMI.out.versions.first())
    ch_bam_grouped = FGBIO_GROUPREADSBYUMI.out.bam

    //
    // MODULE: Run fgbio CallMolecularConsensusReads
    //
    FGBIO_COLLECTDUPLEXSEQMETRICS(ch_bam_grouped, [])
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
    // MODULE: Run fgbio FilterConsensusReads
    //
    FGBIO_FILTERCONSENSUSREADS(ch_consensus_bam, ch_fasta, params.filter_min_reads, params.filter_min_baseq, params.filter_max_base_error_rate)
    ch_versions = ch_versions.mix(FGBIO_FILTERCONSENSUSREADS.out.versions.first())

    //
    // MODULE: Run SamTools fastq
    //
    SAMTOOLS_FASTQ_FIL(FGBIO_FILTERCONSENSUSREADS.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ_FIL.out.versions.first())

    //
    // MODULE: Run STAR on Consensus Data
    //
    STAR_FILBAM(SAMTOOLS_FASTQ_FIL.out.fastq, ch_star_index, ch_gtf, params.star_seq_platform, params.star_seq_center)
    ch_versions = ch_versions.mix(STAR_FILBAM.out.versions.first())
    ch_consensus_bam_filtered = STAR_FILBAM.out.bam

    //
    // MODULE: Run fgbio ZipperBAMs
    //
    FGBIO_ZIPFILBAM(ch_consensus_bam, ch_consensus_bam_filtered, ch_fasta, params.dict)
    ch_versions = ch_versions.mix(FGBIO_ZIPFILBAM.out.versions.first())
    ch_filzip_bam = FGBIO_ZIPFILBAM.out.bam
    
    //
    // MODULE: Extract FastQ reads from BAM
    //
    SAMTOOLS_COLLATEFASTQ(ch_filzip_bam, ch_fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_COLLATEFASTQ.out.versions)
    ch_consensus_reads_filtered = SAMTOOLS_COLLATEFASTQ.out.fastq

    //
    // MODULE: Run fastp
    //
    FASTP(ch_consensus_reads_filtered, params.adapter_fasta, false, false, false)
    ch_reads_all = FASTP.out.reads
    ch_fastp_html = FASTP.out.html
    ch_fastp_json = FASTP.out.json
    ch_versions = ch_versions.mix(FASTP.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(ch_fastp_html.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_fastp_json.collect{it[1]}.ifEmpty([]))

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

    raw_bam         = ch_raw_bam
    versions        = ch_versions
    all_reads       = ch_reads_all
    raw_reads       = ch_cat_fastq
    group_bam       = ch_bam_grouped
    ubam_sorted     = ch_ubam_sorted
    multiqc_files   = ch_multiqc_files 
    consensus_bam   = ch_consensus_bam
    cons_filt_bam   = ch_consensus_bam_filtered
    reads_cons_fil  = ch_consensus_reads_filtered

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                       THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/