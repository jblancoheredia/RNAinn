process WHIPPET_INDEX {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/whippet:1.6.2':
        'blancojmskcc/whippet:1.6.2' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta1), path(gtf)
    
    output:
    tuple val(meta), path("*.jls")                                                          , emit: jls
    tuple val(meta), path("*.exons.tab.gz")                                                 , emit: graph
    tuple val("${task.process}"), val('whippet'), eval("whippet --version"), topic: versions, emit: versions_whippet

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    whippet-index \\
        --fasta ${fasta} \\
        --gtf ${gtf} \\
        $args \\
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args
    
    touch graph.exons.tab.gz
    touch graph.jls
    """
}
