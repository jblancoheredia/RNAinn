process GTFCOLLAPSER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/rnainn_gtfcollapser:1.0.0':
        'blancojmskcc/rnainn_gtfcollapser:1.0.0' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.core.gtf"), emit: core_gtf
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def output = "${gtf.simpleName}.core.gtf"

    """
    GTFcollapser \\
        ${gtf} \\
        ${output} \\
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtfcollapser: 1.0.0
    END_VERSIONS
    """

    stub:
    def output = "${gtf.simpleName}.core.gtf"

    """
    touch ${output}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtfcollapser: 1.0.0
    END_VERSIONS
    """
}
