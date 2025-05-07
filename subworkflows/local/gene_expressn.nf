/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                             GENE_EXPRESSION SUBWORKFLOW                                                    
*****************************************************************************************************************************************************
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/ 

include { STRINGTIE                                                                 } from '../../modules/local/stringtie/main'
include { STAR_INDEX                                                                } from '../../modules/local/star/index/main'
include { STAR_ALIGNX                                                               } from '../../modules/local/star/alignx/main'
include { DESEQ2_QC_RSEM                                                            } from '../../modules/local/rsem/deseq2_qc/main'
include { KALLISTO_QUANT                                                            } from '../../modules/local/kallisto/quant/main'
include { STRINGTIE_MERGE                                                           } from '../../modules/local/stringtie/merge/main'
include { KALLISTO_TX2GENE                                                          } from '../../modules/local/kallisto/tx2gene/main'
include { SALMON_QUANT_BAM                                                          } from '../../modules/local/salmon/quant_bam/main'
include { SALMON_QUANT_FQS                                                          } from '../../modules/local/salmon/quant_fqs/main'
include { KALLISTO_SE_GENE                                                          } from '../../modules/local/kallisto/se_gene/main'
include { RSEM_MERGE_COUNTS                                                         } from '../../modules/local/rsem/merge_counts/main'
include { KALLISTO_TXIMPORT                                                         } from '../../modules/local/kallisto/tximport/main'
include { SALMON_BAM_SE_GENE                                                        } from '../../modules/local/salmon/quant_bam/se_gene/main'
include { SALMON_FQS_SE_GENE                                                        } from '../../modules/local/salmon/quant_fqs/se_gene/main'
include { SALMON_BAM_TX2GENE                                                        } from '../../modules/local/salmon/quant_bam/tx2gene/main'
include { SALMON_FQS_TX2GENE                                                        } from '../../modules/local/salmon/quant_fqs/tx2gene/main'
include { SAMTOOLS_STATS_RSEM                                                       } from '../../modules/local/rsem/samtools/stats/main'
include { SALMON_BAM_TXIMPORT                                                       } from '../../modules/local/salmon/quant_bam/tximport/main'
include { SALMON_FQS_TXIMPORT                                                       } from '../../modules/local/salmon/quant_fqs/tximport/main'
include { DESEQ2_QC_SALMON_BAM                                                      } from '../../modules/local/salmon/quant_bam/deseq2_qc/main'
include { DESEQ2_QC_SALMON_FQS                                                      } from '../../modules/local/salmon/quant_fqs/deseq2_qc/main'
include { SUBREAD_FEATURECOUNTS                                                     } from '../../modules/nf-core/subread/featurecounts/main' 
include { MULTIQC_CUSTOM_BIOTYPE                                                    } from '../../modules/local/multiqc_custom_biotype/main'
include { SAMTOOLS_FLAGSTAT_RSEM                                                    } from '../../modules/local/rsem/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS_RSEM                                                    } from '../../modules/local/rsem/samtools/idxstats/main'
include { KALLISTO_SE_TRANSCRIPT                                                    } from '../../modules/local/kallisto/se_transcript/main'
include { KALLISTO_SE_GENE_SCALED                                                   } from '../../modules/local/kallisto/se_gene_scaled/main'
include { RSEM_CALCULATEEXPRESSION                                                  } from '../../modules/local/rsem/calculateexpression/main'
include { SALMON_BAM_SE_TRANSCRIPT                                                  } from '../../modules/local/salmon/quant_bam/se_transcript/main'
include { SALMON_FQS_SE_TRANSCRIPT                                                  } from '../../modules/local/salmon/quant_fqs/se_transcript/main'
include { SALMON_BAM_SE_GENE_SCALED                                                 } from '../../modules/local/salmon/quant_bam/se_gene_scaled/main'
include { SALMON_FQS_SE_GENE_SCALED                                                 } from '../../modules/local/salmon/quant_fqs/se_gene_scaled/main'
include { KALLISTO_SE_GENE_LENGTH_SCALED                                            } from '../../modules/local/kallisto/se_gene_length_scaled/main'
include { SALMON_BAM_SE_GENE_LENGTH_SCALED                                          } from '../../modules/local/salmon/quant_bam/se_gene_length_scaled/main'
include { SALMON_FQS_SE_GENE_LENGTH_SCALED                                          } from '../../modules/local/salmon/quant_fqs/se_gene_length_scaled/main'

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

workflow GENE_EXPRESSN {

    take:

    ch__gtf
    ch_fasta
    ch_rsem_ref
    ch_versions
    ch_reads_all
    ch_star_index
    ch_samplesheet
    ch_salmon_index
    ch_multiqc_files
    ch_kallisto_index
    ch_pca_header_multiqc
    ch_biotypes_header_multiqc
    ch_clustering_header_multiqc

    main:

    ch_gtf = Channel.fromPath(params.gtf)
    ch_samplesheet = Channel.value(file(params.input, checkIfExists: true))

    //
    // MODULE: Run Kallisto pseudoalignment
    //
    KALLISTO_QUANT (ch_reads_all, ch_kallisto_index, params.gtf, params.chromosomes, [], [])
    ch_pseudo_kallisto_results = KALLISTO_QUANT.out.results
    ch_pseudo_kallisto_multiqc = KALLISTO_QUANT.out.log
    ch_versions = ch_versions.mix(KALLISTO_QUANT.out.versions)

//    //
//    // MODULE: Run TX2GENE for Kallisto
//    //
//    KALLISTO_TX2GENE(ch_gtf, ch_pseudo_kallisto_results, 'kallisto', 'gene_id', 'gene_name')
//    ch_versions = ch_versions.mix(KALLISTO_TX2GENE.out.versions)
//    ch_kallisto_tx2gene = KALLISTO_TX2GENE.out.tx2gene
//
//    //
//    // MODULE: Run TxImport for Kallisto
//    //
//    KALLISTO_TXIMPORT(ch_pseudo_kallisto_results, KALLISTO_TX2GENE.out.tx2gene, 'kallisto')
//    ch_versions = ch_versions.mix(KALLISTO_TXIMPORT.out.versions)
//
//    //
//    // MODULES: Run SummarizedExperiment for Kallisto
//    //
//    KALLISTO_SE_GENE(KALLISTO_TXIMPORT.out.counts_gene.concat(KALLISTO_TXIMPORT.out.tpm_gene).groupTuple(), KALLISTO_TX2GENE.out.tx2gene, ch_samplesheet.map { [ [:], it ] })
//    ch_versions = ch_versions.mix(KALLISTO_SE_GENE.out.versions)
//
//    KALLISTO_SE_GENE_LENGTH_SCALED(KALLISTO_TXIMPORT.out.counts_gene_length_scaled.concat(KALLISTO_TXIMPORT.out.tpm_gene).groupTuple(), KALLISTO_TX2GENE.out.tx2gene, ch_samplesheet.map { [ [:], it ] })
//
//    KALLISTO_SE_GENE_SCALED(KALLISTO_TXIMPORT.out.counts_gene_scaled.concat(KALLISTO_TXIMPORT.out.tpm_gene).groupTuple(), KALLISTO_TX2GENE.out.tx2gene, ch_samplesheet.map { [ [:], it ] })
//
//    KALLISTO_SE_TRANSCRIPT(KALLISTO_TXIMPORT.out.counts_transcript.concat(KALLISTO_TXIMPORT.out.tpm_transcript).groupTuple(), KALLISTO_TX2GENE.out.tx2gene, ch_samplesheet.map { [ [:], it ] })

    //
    // MODULE: Run Salmon pseudoalignment
    //
    SALMON_QUANT_FQS(ch_reads_all, params.salmon_index, params.lib_type)
    ch_versions = ch_versions.mix(SALMON_QUANT_FQS.out.versions)
    ch_salmon_fqs_results = SALMON_QUANT_FQS.out.results
    ch_salmon_fqs_multiqc = SALMON_QUANT_FQS.out.log

    //
    // MODULE: Run TX2GENE for Salmon pseudoalignment
    //
    SALMON_FQS_TX2GENE(ch_gtf, ch_salmon_fqs_results, 'salmon', 'gene_id', 'gene_name')
    ch_versions = ch_versions.mix(SALMON_FQS_TX2GENE.out.versions)

    //
    // MODULE: Run TxImport for Salmon pseudoalignment
    //
    SALMON_FQS_TXIMPORT(ch_salmon_fqs_results, SALMON_FQS_TX2GENE.out.tx2gene, 'salmon')
    ch_versions = ch_versions.mix(SALMON_FQS_TXIMPORT.out.versions)

    //
    // MODULES: Run SummarizedExperiment for Salmon pseudoalignment
    //
    SALMON_FQS_SE_GENE(SALMON_FQS_TXIMPORT.out.counts_gene.concat(SALMON_FQS_TXIMPORT.out.tpm_gene).groupTuple(), SALMON_FQS_TX2GENE.out.tx2gene, ch_samplesheet.map { [ [:], it ] })
    ch_versions = ch_versions.mix(SALMON_FQS_SE_GENE.out.versions)

    SALMON_FQS_SE_GENE_LENGTH_SCALED(SALMON_FQS_TXIMPORT.out.counts_gene_length_scaled.concat(SALMON_FQS_TXIMPORT.out.tpm_gene).groupTuple(), SALMON_FQS_TX2GENE.out.tx2gene, ch_samplesheet.map { [ [:], it ] })

    SALMON_FQS_SE_GENE_SCALED(SALMON_FQS_TXIMPORT.out.counts_gene_scaled.concat(SALMON_FQS_TXIMPORT.out.tpm_gene).groupTuple(), SALMON_FQS_TX2GENE.out.tx2gene, ch_samplesheet.map { [ [:], it ] })

    SALMON_FQS_SE_TRANSCRIPT(SALMON_FQS_TXIMPORT.out.counts_transcript.concat(SALMON_FQS_TXIMPORT.out.tpm_transcript).groupTuple(), SALMON_FQS_TX2GENE.out.tx2gene, ch_samplesheet.map { [ [:], it ] })

    //
    // MODULE: Run Dseq2 QC for Salmon STAR-aligned
    //
    DESEQ2_QC_SALMON_FQS(SALMON_FQS_TXIMPORT.out.counts_gene_length_scaled.map { it[1] }, ch_pca_header_multiqc, ch_clustering_header_multiqc)
    ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_SALMON_FQS.out.pca_multiqc.collect())
    ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_SALMON_FQS.out.dists_multiqc.collect())
    ch_versions = ch_versions.mix(DESEQ2_QC_SALMON_FQS.out.versions)

    //
    // MODULE: Run STAR
    //
    STAR_ALIGNX(ch_reads_all, ch_star_index, ch__gtf, params.star_seq_platform, params.star_seq_center)
    ch_star_sam = STAR_ALIGNX.out.sam
    ch_star_wig = STAR_ALIGNX.out.wig
    ch_star_bam = STAR_ALIGNX.out.bam
    ch_star_tab = STAR_ALIGNX.out.tab
    ch_star_fastq = STAR_ALIGNX.out.fastq
    ch_star_log_out = STAR_ALIGNX.out.log_out
    ch_star_junction = STAR_ALIGNX.out.junction
    ch_star_bedgraph = STAR_ALIGNX.out.bedgraph
    ch_star_log_final = STAR_ALIGNX.out.log_final
    ch_star_bam_sorted = STAR_ALIGNX.out.bam_sorted
    ch_star_log_progress = STAR_ALIGNX.out.log_progress
    ch_star_bam_unsorted = STAR_ALIGNX.out.bam_unsorted
    ch_star_spl_junc_tab = STAR_ALIGNX.out.spl_junc_tab
    ch_versions = ch_versions.mix(STAR_ALIGNX.out.versions)
    ch_star_bam_transcript = STAR_ALIGNX.out.bam_transcript
    ch_star_read_per_gene_tab = STAR_ALIGNX.out.read_per_gene_tab

    //
    // MODULE: Index BAM file
    //
    STAR_INDEX(ch_star_bam_sorted)
    ch_bai = STAR_INDEX.out.bai
    ch_versions = ch_versions.mix(STAR_INDEX.out.versions)
    ch_bam_star_indexed = STAR_ALIGNX.out.bam_sorted.join(STAR_INDEX.out.bai)

    //
    // WORKFLOW: Run StringTie
    //
    STRINGTIE(ch_bam_star_indexed, params.chrgtf)
    ch_versions = ch_versions.mix(STRINGTIE.out.versions)
    ch_stringtie_gtf = STRINGTIE.out.transcript_gtf

    //
    // MODULE: Merge StringTie GTF files
    //
    STRINGTIE_MERGE (ch_stringtie_gtf, params.chrgtf)
    ch_versions = ch_versions.mix(STRINGTIE_MERGE.out.versions)
    ch_stringtie_merged_gtf = STRINGTIE_MERGE.out.gtf

    //
    // MODULE: Run Salmon with STAR alignment
    //
    SALMON_QUANT_BAM(ch_star_bam_transcript, params.transcripts, params.lib_type)
    ch_versions = ch_versions.mix(SALMON_QUANT_BAM.out.versions)
    ch_salmon_bam_results = SALMON_QUANT_BAM.out.results
    ch_salmon_bam_multiqc = SALMON_QUANT_BAM.out.log

    //
    // MODULE: Run TX2GENE for Salmon with STAR alignment
    //
    SALMON_BAM_TX2GENE(ch_gtf, ch_salmon_bam_results, 'salmon', 'gene_id', 'gene_name')
    ch_versions = ch_versions.mix(SALMON_BAM_TX2GENE.out.versions)

    //
    // MODULE: Run TxImport for Salmon with STAR alignment
    //
    SALMON_BAM_TXIMPORT(ch_salmon_bam_results, SALMON_BAM_TX2GENE.out.tx2gene, 'salmon')
    ch_versions = ch_versions.mix(SALMON_BAM_TXIMPORT.out.versions)

    //
    // MODULES: Run SummarizedExperiment for Salmon with STAR alignment
    //
    SALMON_BAM_SE_GENE(SALMON_BAM_TXIMPORT.out.counts_gene.concat(SALMON_BAM_TXIMPORT.out.tpm_gene).groupTuple(), SALMON_BAM_TX2GENE.out.tx2gene, ch_samplesheet.map { [ [:], it ] }) 
    ch_versions = ch_versions.mix(SALMON_BAM_SE_GENE.out.versions)

    SALMON_BAM_SE_GENE_LENGTH_SCALED(SALMON_BAM_TXIMPORT.out.counts_gene_length_scaled.concat(SALMON_BAM_TXIMPORT.out.tpm_gene).groupTuple(), SALMON_BAM_TX2GENE.out.tx2gene, ch_samplesheet.map { [ [:], it ] })

    SALMON_BAM_SE_GENE_SCALED(SALMON_BAM_TXIMPORT.out.counts_gene_scaled.concat(SALMON_BAM_TXIMPORT.out.tpm_gene).groupTuple(), SALMON_BAM_TX2GENE.out.tx2gene, ch_samplesheet.map { [ [:], it ] })

    SALMON_BAM_SE_TRANSCRIPT(SALMON_BAM_TXIMPORT.out.counts_transcript.concat(SALMON_BAM_TXIMPORT.out.tpm_transcript).groupTuple(), SALMON_BAM_TX2GENE.out.tx2gene, ch_samplesheet.map { [ [:], it ] })

    //
    // MODULE: Run Dseq2 QC for Salmon STAR-aligned
    //
    DESEQ2_QC_SALMON_BAM(SALMON_BAM_TXIMPORT.out.counts_gene_length_scaled.map { it[1] }, ch_pca_header_multiqc, ch_clustering_header_multiqc)
    ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_SALMON_BAM.out.pca_multiqc.collect())
    ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_SALMON_BAM.out.dists_multiqc.collect())
    ch_versions = ch_versions.mix(DESEQ2_QC_SALMON_BAM.out.versions)

    //
    // MODULE: Run RSEM CalculateExpression
    //
    RSEM_CALCULATEEXPRESSION(ch_star_bam_transcript, params.rsem_ref)
    ch_versions = ch_versions.mix(RSEM_CALCULATEEXPRESSION.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(RSEM_CALCULATEEXPRESSION.out.stat.collect{it[1]})
    ch_bam_bai_transcript = RSEM_CALCULATEEXPRESSION.out.bam_transcript_sorted.join(RSEM_CALCULATEEXPRESSION.out.bai_transcript_sorted)

    //
    // MODULE: Run SamTools Stats on RSEM genome sorted indexed BAM file
    //
    SAMTOOLS_STATS_RSEM(ch_bam_bai_transcript, ch_fasta)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_RSEM.out.stats.collect{it[1]})
    ch_versions = ch_versions.mix(SAMTOOLS_STATS_RSEM.out.versions)

    //
    // MODULE: Run SamTools FlagStat on RSEM genome sorted indexed BAM file
    //
    SAMTOOLS_FLAGSTAT_RSEM ( ch_bam_bai_transcript )
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_FLAGSTAT_RSEM.out.flagstat.collect{it[1]})
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT_RSEM.out.versions)

    //
    // MODULE: Run SamTools IdxStats on RSEM genome sorted indexed BAM file
    //
    SAMTOOLS_IDXSTATS_RSEM ( ch_bam_bai_transcript )
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_IDXSTATS_RSEM.out.idxstats.collect{it[1]})
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS_RSEM.out.versions)

    //
    // MODULE: Run 
    //
    RSEM_MERGE_COUNTS(RSEM_CALCULATEEXPRESSION.out.counts_gene, RSEM_CALCULATEEXPRESSION.out.counts_transcript)
    ch_versions = ch_versions.mix(RSEM_MERGE_COUNTS.out.versions)

    //
    // MODULE: Run Dseq2 QC for RSEM
    //
    DESEQ2_QC_RSEM(RSEM_MERGE_COUNTS.out.counts_gene, ch_pca_header_multiqc, ch_clustering_header_multiqc)
    ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_RSEM.out.pca_multiqc.collect())
    ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_RSEM.out.dists_multiqc.collect())
    ch_versions = ch_versions.mix(DESEQ2_QC_RSEM.out.versions)

    //
    // SUBMODULE: Run SubReadFeatureCounts
    //
    def biotype = params.featurecounts_group_type
    ch_star_bam
        .combine(ch_gtf)
        .set { ch_featurecounts }

    SUBREAD_FEATURECOUNTS(ch_featurecounts)
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())
    ch_subread_feature_counts = SUBREAD_FEATURECOUNTS.out.counts

    MULTIQC_CUSTOM_BIOTYPE(ch_subread_feature_counts, ch_biotypes_header_multiqc)
    ch_multiqc_files = ch_multiqc_files.mix(MULTIQC_CUSTOM_BIOTYPE.out.tsv.collect{it[1]})
    ch_versions = ch_versions.mix(MULTIQC_CUSTOM_BIOTYPE.out.versions.first())

    emit:

    versions                                = ch_versions
    multiqc_files                           = ch_multiqc_files 
//    kallisto_tpm_gene                       = KALLISTO_TXIMPORT.out.tpm_gene
//    kallisto_counts_gene                    = KALLISTO_TXIMPORT.out.counts_gene
//    kallisto_lengths_gene                   = KALLISTO_TXIMPORT.out.lengths_gene
//    kallisto_counts_gene_length_scaled      = KALLISTO_TXIMPORT.out.counts_gene_length_scaled
//    kallisto_counts_gene_scaled             = KALLISTO_TXIMPORT.out.counts_gene_scaled
//    kallisto_tpm_transcript                 = KALLISTO_TXIMPORT.out.tpm_transcript
//    kallisto_counts_transcript              = KALLISTO_TXIMPORT.out.counts_transcript
//    kallisto_lengths_transcript             = KALLISTO_TXIMPORT.out.lengths_transcript
//    kallisto_merged_gene_rds                = KALLISTO_SE_GENE.out.rds
//    kallisto_merged_gene_rds_length_scaled  = KALLISTO_SE_GENE_LENGTH_SCALED.out.rds
//    kallisto_merged_gene_rds_scaled         = KALLISTO_SE_GENE_SCALED.out.rds
//    kallisto_merged_transcript_rds          = KALLISTO_SE_TRANSCRIPT.out.rds
//    salmon_fqs_tpm_gene                         = SALMON_FQS_TXIMPORT.out.tpm_gene
//    salmon_fqs_counts_gene                      = SALMON_FQS_TXIMPORT.out.counts_gene
//    salmon_fqs_lengths_gene                     = SALMON_FQS_TXIMPORT.out.lengths_gene
//    salmon_fqs_counts_gene_length_scaled        = SALMON_FQS_TXIMPORT.out.counts_gene_length_scaled
//    salmon_fqs_counts_gene_scaled               = SALMON_FQS_TXIMPORT.out.counts_gene_scaled
//    salmon_fqs_tpm_transcript                   = SALMON_FQS_TXIMPORT.out.tpm_transcript
//    salmon_fqs_counts_transcript                = SALMON_FQS_TXIMPORT.out.counts_transcript
//    salmon_fqs_lengths_transcript               = SALMON_FQS_TXIMPORT.out.lengths_transcript
//    salmon_fqs_merged_gene_rds                  = SALMON_FQS_SE_GENE.out.rds
//    salmon_fqs_merged_gene_rds_length_scaled    = SALMON_FQS_SE_GENE_LENGTH_SCALED.out.rds
//    salmon_fqs_merged_gene_rds_scaled           = SALMON_FQS_SE_GENE_SCALED.out.rds
//    salmon_fqs_merged_transcript_rds            = SALMON_FQS_SE_TRANSCRIPT.out.rds
//    salmon_bam_tpm_gene                         = SALMON_BAM_TXIMPORT.out.tpm_gene
//    salmon_bam_counts_gene                      = SALMON_BAM_TXIMPORT.out.counts_gene
//    salmon_bam_lengths_gene                     = SALMON_BAM_TXIMPORT.out.lengths_gene
//    salmon_bam_counts_gene_length_scaled        = SALMON_BAM_TXIMPORT.out.counts_gene_length_scaled
//    salmon_bam_counts_gene_scaled               = SALMON_BAM_TXIMPORT.out.counts_gene_scaled
//    salmon_bam_tpm_transcript                   = SALMON_BAM_TXIMPORT.out.tpm_transcript
//    salmon_bam_counts_transcript                = SALMON_BAM_TXIMPORT.out.counts_transcript
//    salmon_bam_lengths_transcript               = SALMON_BAM_TXIMPORT.out.lengths_transcript
//    salmon_bam_merged_gene_rds                  = SALMON_BAM_SE_GENE.out.rds
//    salmon_bam_merged_gene_rds_length_scaled    = SALMON_BAM_SE_GENE_LENGTH_SCALED.out.rds
//    salmon_bam_merged_gene_rds_scaled           = SALMON_BAM_SE_GENE_SCALED.out.rds
//    salmon_bam_merged_transcript_rds            = SALMON_BAM_SE_TRANSCRIPT.out.rds
//    subread_feature_counts                  = ch_subread_feature_counts
//    merged_counts_gene       = RSEM_MERGE_COUNTS.out.counts_gene              //    path: *.gene_counts.tsv
//    merged_tpm_gene          = RSEM_MERGE_COUNTS.out.tpm_gene                 //    path: *.gene_tpm.tsv
//    merged_counts_transcript = RSEM_MERGE_COUNTS.out.counts_transcript        //    path: *.transcript_counts.tsv
//    merged_tpm_transcript    = RSEM_MERGE_COUNTS.out.tpm_transcript           //    path: *.transcript_tpm.tsv
    star_bam_transcript                     = ch_star_bam_transcript
    star_bam_indexed                        = ch_bam_star_indexed
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                       THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/