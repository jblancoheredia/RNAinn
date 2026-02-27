process GTFCOLLAPSER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'docker://community.wave.seqera.io/library/pip_bx-python_pandas:677de3d2fb4d8a9e':
        'community.wave.seqera.io/library/pip_bx-python_pandas:677de3d2fb4d8a9e' }"

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
