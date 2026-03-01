process WHIPPET_QUANT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/whippet:1.6.2':
        'blancojmskcc/whippet:1.6.2' }"

    input:
    tuple val(meta), path(inputs)

    output:
    tuple val(meta), path("*.psi.gz"), emit: psi
    tuple val(meta), path("*.tpm.gz"), emit: tmp
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def jls    = inputs[2]
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads1 = inputs[0]
    def reads2 = inputs[1]
    """
    export JULIA_DEPOT_PATH="${PWD}/.julia:/opt/julia_depot"
    
    whippet-quant \\
        ${reads1} \\
        ${reads2} \\
        -x ${jls} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whippet-index: 1.6.2
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """  
    touch ${prefix}.psi.gz
    touch ${prefix}.tmp.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whippet-index: 1.6.2
    END_VERSIONS
    """
}
