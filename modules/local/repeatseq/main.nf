process REPEATSEQ {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/repeatseq:0.8.2' :
        'blancojmskcc/repeatseq:0.8.2' }"

    input:
    tuple val(meta),  val(chunk), path(bam)
    tuple val(meta0), val(chunk), path(bai)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fai)
    path(rep_regions)

    output:
    tuple val(meta), path("*_repeatseq.vcf"), optional: true,       emit: vcf
    tuple val(meta), path("*_repeatseq.calls"), optional: true,     emit: calls
    tuple val(meta), path("*_repeatseq.repeatseq"), optional: true, emit: repeatseq
    path "versions.yml",                                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk = bam.name.replace(meta.id + "_", "").replace(".bam", "")
    """
    repeatseq \\
        -repeatseq \\
        -calls \\
        ${bam} \\
        ${fasta} \\
        ${rep_regions}
    
    if [ \$? -eq 0 ]; then
        mv *.vcf ${prefix}_${chunk}_repeatseq.vcf
        mv *.calls ${prefix}_${chunk}_repeatseq.calls
        mv *.repeatseq ${prefix}_${chunk}_repeatseq.repeatseq
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatseq: \$(repeatseq -v 2>&1 | grep "RepeatSeq v" | sed 's/.*RepeatSeq v//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk = bam.name.replace(meta.id + "_", "").replace(".bam", "")
    """
    touch ${prefix}_${chunk}_repeatseq.vcf
    touch ${prefix}_${chunk}_repeatseq.calls
    touch ${prefix}_${chunk}_repeatseq.repeatseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatseq: \$(repeatseq -v 2>&1 | grep "RepeatSeq v" | sed 's/.*RepeatSeq v//')
    END_VERSIONS
    """
}