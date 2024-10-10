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
include { STAR_ALIGN                                                                } from '../../modules/nf-core/star/align/main'
include { DESEQ2_QC_RSEM                                                            } from '../../modules/local/deseq2_qc_rsem/main'
include { SALMON_TX2GENE                                                            } from '../../modules/local/salmon/tx2gene/main'
include { SALMON_SE_GENE                                                            } from '../../modules/local/salmon/se_gene/main'
include { KALLISTO_QUANT                                                            } from '../../modules/nf-core/kallisto/quant/main'
include { SALMON_TXIMPORT                                                           } from '../../modules/local/salmon/tximport/main'
include { STRINGTIE_MERGE                                                           } from '../../modules/nf-core/stringtie/merge/main'
include { KALLISTO_TX2GENE                                                          } from '../../modules/local/kallisto/tx2gene/main'
include { SALMON_QUANT_FQS                                                          } from '../../modules/local/salmon/quant_fqs/main'
include { RSEM_MERGE_COUNTS                                                         } from '../../modules/local/rsem_merge_counts/main'
include { KALLISTO_TXIMPORT                                                         } from '../../modules/local/kallisto/tximport/main'
include { KALLISTO_SE_GENE                                                          } from '../../modules/local/kallisto/se_gene/main'
include { SALMON_SE_TRANSCRIPT                                                      } from '../../modules/local/salmon/se_transcript/main'
include { SALMON_SE_GENE_SCALED                                                     } from '../../modules/local/salmon/se_gene_scaled/main'
include { SUBREAD_FEATURECOUNTS                                                     } from '../../modules/nf-core/subread/featurecounts/main' 
include { MULTIQC_CUSTOM_BIOTYPE                                                    } from '../../modules/local/multiqc_custom_biotype'
include { KALLISTO_SE_TRANSCRIPT                                                    } from '../../modules/local/kallisto/se_transcript/main'
include { KALLISTO_SE_GENE_SCALED                                                   } from '../../modules/local/kallisto/se_gene_scaled/main'
include { RSEM_CALCULATEEXPRESSION                                                  } from '../../modules/nf-core/rsem/calculateexpression'
include { SALMON_SE_GENE_LENGTH_SCALED                                              } from '../../modules/local/salmon/se_gene_length_scaled/main'
include { KALLISTO_SE_GENE_LENGTH_SCALED                                            } from '../../modules/local/kallisto/se_gene_length_scaled/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                 IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { biotypeInGtf                                                              } from './utils_nfcore_rnainn_pipeline'
include { QUANTIFY_RSEM                                                             } from './quantify_rsem/main'
include { softwareVersionsToYAML                                                    } from '../nf-core/utils_nfcore_pipeline'
include { BAM_SORT_STATS_SAMTOOLS                                                   } from '../nf-core/bam_sort_stats_samtools'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                   RUN SUBWORKFLOW 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNAINN_GENE_EXPRESSN {

    take:

    ch_gtf
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

    //
    // MODULE: Run Kallisto pseudoalignment
    //
    KALLISTO_QUANT (ch_reads_all, ch_kallisto_index, params.gtf, params.chromosomes, [], [])
    ch_pseudo_kallisto_results = KALLISTO_QUANT.out.results
    ch_pseudo_kallisto_multiqc = KALLISTO_QUANT.out.log
    ch_versions = ch_versions.mix(KALLISTO_QUANT.out.versions)

    //
    // MODULE: Run TX2GENE for Kallisto
    //
    KALLISTO_TX2GENE(ch_gtf, ch_pseudo_kallisto_results, 'kallisto', 'gene_id', 'transcript_id')
    ch_versions = ch_versions.mix(KALLISTO_TX2GENE.out.versions)

    //
    // MODULE: Run TxImport for Kallisto
    //
    KALLISTO_TXIMPORT(ch_pseudo_kallisto_results.collect{ it[1] }.map { [ ['id': 'all_samples'], it ] }, KALLISTO_TX2GENE.out.tx2gene, 'kallisto')
    ch_versions = ch_versions.mix(KALLISTO_TXIMPORT.out.versions)

    //
    // MODULES: Run SummarizedExperiment for Kallisto
    //
    KALLISTO_SE_GENE(KALLISTO_TXIMPORT.out.counts_gene.concat(KALLISTO_TXIMPORT.out.tpm_gene).groupTuple(), KALLISTO_TX2GENE.out.tx2gene, ch_samplesheet)
    ch_versions = ch_versions.mix(KALLISTO_SE_GENE.out.versions)

    KALLISTO_SE_GENE_LENGTH_SCALED(KALLISTO_TXIMPORT.out.counts_gene_length_scaled.concat(KALLISTO_TXIMPORT.out.tpm_gene).groupTuple(), KALLISTO_TX2GENE.out.tx2gene, ch_samplesheet)

    KALLISTO_SE_GENE_SCALED(KALLISTO_TXIMPORT.out.counts_gene_scaled.concat(KALLISTO_TXIMPORT.out.tpm_gene).groupTuple(), KALLISTO_TX2GENE.out.tx2gene, ch_samplesheet)

    KALLISTO_SE_TRANSCRIPT(KALLISTO_TXIMPORT.out.counts_transcript.concat(KALLISTO_TXIMPORT.out.tpm_transcript).groupTuple(), KALLISTO_TX2GENE.out.tx2gene, ch_samplesheet)

    //
    // MODULE: Run Salmon pseudoalignment
    //
    SALMON_QUANT_FQS(ch_reads_all, params.salmon_index, params.lib_type)
    ch_versions = ch_versions.mix(SALMON_QUANT_FQS.out.versions)
    ch_pseudo_salmon_results = SALMON_QUANT_FQS.out.results.collect{ it[1] }.map { [ [:], it ] }
    ch_pseudo_salmon_multiqc = ch_pseudo_salmon_results

    //
    // MODULE: Run TX2GENE for Salmon
    //
    SALMON_TX2GENE(ch_gtf, ch_pseudo_salmon_results, 'salmon', 'gene_id', 'transcript_id')
    ch_versions = ch_versions.mix(SALMON_TX2GENE.out.versions)

    //
    // MODULE: Run TxImport for Salmon
    //
    SALMON_TXIMPORT(ch_pseudo_salmon_results.collect{ it[1] }.map { [ ['id': 'all_samples'], it ] }, SALMON_TX2GENE.out.tx2gene, 'salmon')
    ch_versions = ch_versions.mix(SALMON_TXIMPORT.out.versions)

    //
    // MODULES: Run SummarizedExperiment for Salmon
    //
    SALMON_SE_GENE(SALMON_TXIMPORT.out.counts_gene.concat(SALMON_TXIMPORT.out.tpm_gene).groupTuple(), SALMON_TX2GENE.out.tx2gene, ch_samplesheet)
    ch_versions = ch_versions.mix(SALMON_SE_GENE.out.versions)

    SALMON_SE_GENE_LENGTH_SCALED(SALMON_TXIMPORT.out.counts_gene_length_scaled.concat(SALMON_TXIMPORT.out.tpm_gene).groupTuple(), SALMON_TX2GENE.out.tx2gene, ch_samplesheet)

    SALMON_SE_GENE_SCALED(SALMON_TXIMPORT.out.counts_gene_scaled.concat(SALMON_TXIMPORT.out.tpm_gene).groupTuple(), SALMON_TX2GENE.out.tx2gene, ch_samplesheet)

    SALMON_SE_TRANSCRIPT(SALMON_TXIMPORT.out.counts_transcript.concat(SALMON_TXIMPORT.out.tpm_transcript).groupTuple(), SALMON_TX2GENE.out.tx2gene, ch_samplesheet)

//    //
//    // SUBWORKFLOW: Run 
//    //
//    QUANTIFY_RSEM(ch_reads_all, params.rsem_ref, ch_fasta)
//    ch_genome_bam       = QUANTIFY_RSEM.out.bam
//    ch_genome_bam_index = QUANTIFY_RSEM.out.bai
//    ch_star_log         = QUANTIFY_RSEM.out.logs
//    ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_RSEM.out.stats.collect{it[1]})
//    ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_RSEM.out.flagstat.collect{it[1]})
//    ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_RSEM.out.idxstats.collect{it[1]})
//    ch_multiqc_files = ch_multiqc_files.mix(ch_star_log.collect{it[1]})
//    ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_RSEM.out.stat.collect{it[1]})
//    ch_versions = ch_versions.mix(QUANTIFY_RSEM.out.versions)
//
//    //
//    // MODULE: Run Dseq2 QC
//    //
//    DESEQ2_QC_RSEM(QUANTIFY_RSEM.out.merged_counts_gene, ch_pca_header_multiqc, ch_clustering_header_multiqc)
//    ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_RSEM.out.pca_multiqc.collect())
//    ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_RSEM.out.dists_multiqc.collect())
//    ch_versions = ch_versions.mix(DESEQ2_QC_RSEM.out.versions)

//    //
//    // MODULE: Run STAR
//    //
//    STAR_ALIGN(ch_reads_all, ch_star_index, ch_gtf, params.star_seq_platform, params.star_seq_center)
//    ch_star_sam = STAR_ALIGN.out.sam
//    ch_star_wig = STAR_ALIGN.out.wig
//    ch_star_bam = STAR_ALIGN.out.bam
//    ch_star_tab = STAR_ALIGN.out.tab
//    ch_star_fastq = STAR_ALIGN.out.fastq
//    ch_star_log_out = STAR_ALIGN.out.log_out
//    ch_star_junction = STAR_ALIGN.out.junction
//    ch_star_bedgraph = STAR_ALIGN.out.bedgraph
//    ch_star_log_final = STAR_ALIGN.out.log_final
//    ch_star_bam_sorted = STAR_ALIGN.out.bam_sorted
//    ch_star_log_progress = STAR_ALIGN.out.log_progress
//    ch_star_bam_unsorted = STAR_ALIGN.out.bam_unsorted
//    ch_star_spl_junc_tab = STAR_ALIGN.out.spl_junc_tab
//    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)
//    ch_star_bam_transcript = STAR_ALIGN.out.bam_transcript
//    ch_star_read_per_gene_tab = STAR_ALIGN.out.read_per_gene_tab
//
//    //
//    // MODULE: Index BAM file
//    //
//    STAR_INDEX(ch_star_bam_sorted)
//    ch_bai = STAR_INDEX.out.bai
//    ch_versions = ch_versions.mix(STAR_INDEX.out.versions)
//    ch_bam_star_indexed = STAR_ALIGN.out.bam_sorted.join(STAR_INDEX.out.bai)
//
//    //
//    // WORKFLOW: Run StringTie
//    //
//    STRINGTIE(ch_bam_star_indexed, params.chrgtf)
//    ch_versions = ch_versions.mix(STRINGTIE.out.versions)
//    STRINGTIE
//        .out
//        .transcript_gtf
//        .map { it -> it[1] }
//        .set { stringtie_gtf }
//    ch_versions = ch_versions.mix(STRINGTIE.out.versions)
//
//    //
//    // MODULE: Merge StringTie GTF files
//    //
//    STRINGTIE_MERGE (stringtie_gtf, params.chrgtf)
//    ch_versions = ch_versions.mix(STRINGTIE_MERGE.out.versions)
//    ch_stringtie_gtf = STRINGTIE_MERGE.out.gtf

//    //
//    // SUBMODULE: Run SubReadFeatureCounts
//    //
//    def biotype = params.gencode ? "gene_type" : params.featurecounts_group_type
//    ch_gtf.map{biotypeInGtf(it, biotype)}.set {biotype_in_gtf}
//    // Prevent any samples from running if GTF file doesn't have a valid biotype
//    ch_star_bam
//        .combine(ch_gtf)
//        .combine(biotype_in_gtf)
//        .filter { it[-1] }
//        .map { it[0..<it.size()-1] }
//        .set { ch_featurecounts }
//
//    SUBREAD_FEATURECOUNTS(ch_featurecounts)
//    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())
//    ch_subread_feature_counts = SUBREAD_FEATURECOUNTS.out.counts
//
//    MULTIQC_CUSTOM_BIOTYPE(ch_subread_feature_counts, ch_biotypes_header_multiqc)
//    ch_multiqc_files = ch_multiqc_files.mix(MULTIQC_CUSTOM_BIOTYPE.out.tsv.collect{it[1]})
//    ch_versions = ch_versions.mix(MULTIQC_CUSTOM_BIOTYPE.out.versions.first())

    emit:

    versions                                = ch_versions
    multiqc_files                           = ch_multiqc_files 
    kallisto_tpm_gene                       = KALLISTO_TXIMPORT.out.tpm_gene
    kallisto_counts_gene                    = KALLISTO_TXIMPORT.out.counts_gene
    kallisto_lengths_gene                   = KALLISTO_TXIMPORT.out.lengths_gene
    kallisto_counts_gene_length_scaled      = KALLISTO_TXIMPORT.out.counts_gene_length_scaled
    kallisto_counts_gene_scaled             = KALLISTO_TXIMPORT.out.counts_gene_scaled
    kallisto_tpm_transcript                 = KALLISTO_TXIMPORT.out.tpm_transcript
    kallisto_counts_transcript              = KALLISTO_TXIMPORT.out.counts_transcript
    kallisto_lengths_transcript             = KALLISTO_TXIMPORT.out.lengths_transcript
    kallisto_merged_gene_rds                = KALLISTO_SE_GENE.out.rds
    kallisto_merged_gene_rds_length_scaled  = KALLISTO_SE_GENE_LENGTH_SCALED.out.rds
    kallisto_merged_gene_rds_scaled         = KALLISTO_SE_GENE_SCALED.out.rds
    kallisto_merged_transcript_rds          = KALLISTO_SE_TRANSCRIPT.out.rds
    salmon_tpm_gene                         = SALMON_TXIMPORT.out.tpm_gene
    salmon_counts_gene                      = SALMON_TXIMPORT.out.counts_gene
    salmon_lengths_gene                     = SALMON_TXIMPORT.out.lengths_gene
    salmon_counts_gene_length_scaled        = SALMON_TXIMPORT.out.counts_gene_length_scaled
    salmon_counts_gene_scaled               = SALMON_TXIMPORT.out.counts_gene_scaled
    salmon_tpm_transcript                   = SALMON_TXIMPORT.out.tpm_transcript
    salmon_counts_transcript                = SALMON_TXIMPORT.out.counts_transcript
    salmon_lengths_transcript               = SALMON_TXIMPORT.out.lengths_transcript
    salmon_merged_gene_rds                  = SALMON_SE_GENE.out.rds
    salmon_merged_gene_rds_length_scaled    = SALMON_SE_GENE_LENGTH_SCALED.out.rds
    salmon_merged_gene_rds_scaled           = SALMON_SE_GENE_SCALED.out.rds
    salmon_merged_transcript_rds            = SALMON_SE_TRANSCRIPT.out.rds
//    subread_feature_counts                  = ch_subread_feature_counts
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                       THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/