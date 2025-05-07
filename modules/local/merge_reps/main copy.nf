process MERGE_REPS {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta),  path(vcfs)
    tuple val(meta1), path(calls)
    tuple val(meta2), path(repeatseqs)

    output:
    tuple val(meta), path("${meta.id}_repeatseq.vcf"),        emit: vcf
    tuple val(meta), path("${meta.id}_repeatseq.calls"),      emit: calls
    tuple val(meta), path("${meta.id}_repeatseq.repeatseq"),  emit: repeatseq

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat \$(head -8 ${vcfs[0]}) > ${prefix}_repeatseq.vcf
    for vcf in ${vcfs[*]}; do
        grep -v "^#" \$vcf >> ${prefix}_repeatseq.vcf
    done

    # Merge calls files
    cat ${calls[*]} > ${prefix}_repeatseq.calls

    # Merge repeatseq files
    cat ${repeatseqs[*]} > ${prefix}_repeatseq.repeatseq
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_repeatseq.vcf
    touch ${prefix}_repeatseq.calls
    touch ${prefix}_repeatseq.repeatseq
    """
}
