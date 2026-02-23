process GTFCOLLAPSER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.70--np112py36_0':
        'biocontainers/biopython:1.70--np112py36_0' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.core.gtf.gz"), emit: core_gtf
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def output = gtf.getName().replaceFirst(/\.gtf\.gz$/, '.core.gtf.gz')

    """
    collapse_annotation.py \\
        ${gtf} \\
        ${output} \\
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtfcollapser: 1.0.0
    END_VERSIONS
    """

    stub:
    def output = gtf.getName().replaceFirst(/\.gtf\.gz$/, '.core.gtf.gz')

    """
    touch ${output}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtfcollapser: 1.0.0
    END_VERSIONS
    """
}
