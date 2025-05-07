process REPEATSEQ {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/repeatseq:0.8.2' :
        'blancojmskcc/repeatseq:0.8.2' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fai)
    path(rep_regions)

    output:
    tuple val(meta),  path("*.vcf"),        emit: vcf
    tuple val(meta),  path("*.calls"),      emit: calls
    tuple val(meta2), path("*.repeatseq"),  emit: repeatseq
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    repeatseq \\
        -repeatseq \\
        -calls \\
        ${bam} \\
        ${fasta} \\
        ${rep_regions}

    mv *.vcf ${prefix}_repeatseq.vcf
    mv *.calls ${prefix}_repeatseq.calls
    mv *.repeatseq ${prefix}_repeatseq.repeatseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatseq: \$(repeatseq -v 2>&1 | grep "RepeatSeq v" | sed 's/.*RepeatSeq v//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_repeatseq.vcf
    touch ${prefix}_repeatseq.calls
    touch ${prefix}_repeatseq.repeatseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatseq: \$(repeatseq -v 2>&1 | grep "RepeatSeq v" | sed 's/.*RepeatSeq v//')
    END_VERSIONS
    """
}    
