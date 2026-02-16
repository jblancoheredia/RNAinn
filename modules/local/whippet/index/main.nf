process WHIPPET_INDEX {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/whippet:1.6.2':
        'blancojmskcc/whippet:1.6.2' }"

    input:
    tuple val(meta) , path(bam), path(bai)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(gtf)
    
    output:
    tuple val(meta), path("*.jls")         , emit: jls
    tuple val(meta), path("*.exons.tab.gz"), emit: graph
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    whippet-index \\
        --fasta ${fasta} \\
        --bam ${bam} \\
        --gtf ${gtf} \\
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whippet-index: 1.6.2
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args
    
    touch graph.exons.tab.gz
    touch graph.jls

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whippet-index: 1.6.2
    END_VERSIONS
    """
}
