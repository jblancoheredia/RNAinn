/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                       FILTER_AND_SET_STRANDEDNESS SUBWORKFLOW                                                    
*****************************************************************************************************************************************************
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/ 

include { SORTMERNA                                                                 } from '../modules/nf-core/sortmerna/main' 
include { STAR_ALIGN                                                                } from '../modules/nf-core/star/align/main'
include { BBMAP_BBSPLIT                                                             } from '../modules/nf-core/bbmap/bbsplit/main' 


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                 IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML                                                    } from '../nf-core/utils_nfcore_pipeline'
include { FASTQ_SUBSAMPLE_FQ_SALMON                                                 } from '../nf-core/fastq_subsample_fq_salmon/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                               IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Import libraries
//
import groovy.json.JsonSlurper

//
// Function to determine library type by comparing type counts.
//
def calculateStrandedness(forwardFragments, reverseFragments, unstrandedFragments, stranded_threshold=0.8, unstranded_threshold=0.1) {
    def totalFragments = forwardFragments + reverseFragments + unstrandedFragments
    def totalStrandedFragments = forwardFragments + reverseFragments
    def library_strandedness = 'undetermined'
    if (totalStrandedFragments > 0) {
        def forwardProportion = forwardFragments / (totalStrandedFragments as double)
        def reverseProportion = reverseFragments / (totalStrandedFragments as double)
        def proportionDifference = Math.abs(forwardProportion - reverseProportion)
        if (forwardProportion >= stranded_threshold) {
            strandedness = 'forward'
        } else if (reverseProportion >= stranded_threshold) {
            strandedness = 'reverse'
        } else if (proportionDifference <= unstranded_threshold) {
            strandedness = 'unstranded'
        }
    }
    return [
        inferred_strandedness: strandedness,
        forwardFragments: (forwardFragments / (totalFragments as double)) * 100,
        reverseFragments: (reverseFragments / (totalFragments as double)) * 100,
        unstrandedFragments: (unstrandedFragments / (totalFragments as double)) * 100
    ]
}

//
// Function that parses Salmon quant 'lib_format_counts.json' output file to get inferred strandedness
//
def getSalmonInferredStrandedness(json_file, stranded_threshold = 0.8, unstranded_threshold = 0.1) {
    def libCounts = new JsonSlurper().parseText(json_file.text)
    def forwardKeys = ['SF', 'ISF', 'MSF', 'OSF']
    def reverseKeys = ['SR', 'ISR', 'MSR', 'OSR']
    def unstrandedKeys = ['IU', 'U', 'MU']
    def forwardFragments = forwardKeys.collect { libCounts[it] ?: 0 }.sum()
    def reverseFragments = reverseKeys.collect { libCounts[it] ?: 0 }.sum()
    def unstrandedFragments = unstrandedKeys.collect { libCounts[it] ?: 0 }.sum()
    return calculateStrandedness(forwardFragments, reverseFragments, unstrandedFragments, stranded_threshold, unstranded_threshold)
}

def pass_trimmed_reads = [:]

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                   RUN SUBWORKFLOW 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNAINN_FILTER_STRAND {

    take:

    ch_versions,
    ch_multiqc_files

    //
    // MODULE: Run BBTools to remove contaminant reads
    //
    BBMAP_BBSPLIT(ch_reads_all, ch_bbsplit_index, [], [ [], [] ], false)
    BBMAP_BBSPLIT.out.primary_fastq.set{ch_reads_all}
    ch_versions = ch_versions.mix(BBMAP_BBSPLIT.out.versions.first())

    //
    // MODULE: Run SortMeRNA
    //
    SORTMERNA(ch_reads_all, ch_sortmerna_fastas, ch_sortmerna_index)
    SORTMERNA.out.reads.set{ch_reads_all}
    ch_multiqc_files = ch_multiqc_files.mix(SORTMERNA.out.log)
    ch_versions = ch_versions.mix(SORTMERNA.out.versions.first())

    //
    // BRANCH: Run if strandedness was specified as 'auto'  
    //
    ch_reads_all
        .branch {
            meta, fastq ->
                auto_strand : meta.strandedness == 'auto'
                    return [ meta, fastq ]
                known_strand: meta.strandedness != 'auto'
                    return [ meta, fastq ]
        }
        .set{ch_strand_fastq}
    
    //
    // SUBWORKFLOW: Sub-sample FastQ files and pseudoalign with Salmon to auto-infer strandedness
    //
    // Return empty channel if ch_strand_fastq.auto_strand is empty so salmon index isn't created
    //
    ch_fasta
        .combine(ch_strand_fastq.auto_strand)
        .map { it.first() }
        .first()
        .set { ch_genome_fasta }

    FASTQ_SUBSAMPLE_FQ_SALMON (
        ch_strand_fastq.auto_strand,
        ch_genome_fasta,
        ch_transcript_fasta,
        ch_gtf,
        ch_salmon_index,
        make_salmon_index
    )
    ch_versions = ch_versions.mix(FASTQ_SUBSAMPLE_FQ_SALMON.out.versions)

    FASTQ_SUBSAMPLE_FQ_SALMON
        .out
        .lib_format_counts
        .join(ch_strand_fastq.auto_strand)
        .map {
            meta, json, reads ->
                def salmon_strand_analysis = getSalmonInferredStrandedness(json, stranded_threshold=stranded_threshold, unstranded_threshold=unstranded_threshold)
                strandedness = salmon_strand_analysis.inferred_strandedness
                if (strandedness == 'undetermined') {
                    strandedness = 'unstranded'
                }
                return [ meta + [ strandedness: strandedness, salmon_strand_analysis: salmon_strand_analysis ], reads ]
        }
        .mix(ch_strand_fastq.known_strand)
        .set { ch_strand_inferred_fastq }

    emit:

    reads           = ch_strand_inferred_fastq
    trim_read_count = ch_trim_read_count

    multiqc_files   = ch_multiqc_files.transpose().map{it[1]}
    versions        = ch_versions                     // channel: [ versions.yml ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                       THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/