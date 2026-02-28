process WHIPPET_INDEX {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/whippet:1.6.2':
        'blancojmskcc/whippet:1.6.2' }"

    input:
    tuple val(meta) , path(bam, stageAs:'bam_input'), path(bai, stageAs:'bai_input')
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
    def prefix = meta ? "${meta.id}_" : ""
    def bam_input = bam ? "--bam ${bam}" : ""
    """
    export JULIA_DEPOT_PATH="${PWD}/.julia:/opt/julia_depot"
    
    whippet-index \\
        --fasta ${fasta} \\
        --gtf ${gtf} \\
        ${bam_input} \\
        -x ${prefix}graph.jls \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whippet-index: 1.6.2
    END_VERSIONS
    """
    stub:
    def prefix = meta ? "${meta.id}" : ""

    """  
    touch ${prefix}_graph.jls.exons.tab.gz
    touch ${prefix}_graph.jls

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whippet-index: 1.6.2
    END_VERSIONS
    """
}