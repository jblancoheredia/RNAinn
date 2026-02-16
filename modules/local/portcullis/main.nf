process PORTCULLIS {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::portcullis=1.2.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/portcullis:1.2.4--py38haf070c8_0':
        'biocontainers/portcullis:1.2.4--py38haf070c8_0' }"

    input:
    tuple val(meta), path(bam)
    path(bed)
    path(fasta)

    output:
    tuple val(meta), path("*.bed")                                   , emit: bed
    tuple val(meta), path("*.log")                                   , emit: log
    tuple val(meta), path("*.tab")                                   , emit: tab
    tuple val(meta), path("*.spliced.bam"), path("*.spliced.bam.bai"), emit: bam_bai
    path "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def log_file = "${meta.id}.portcullis.log"
    """
    portcullis \\
        full \\
        ${args} \\
        -t ${task.cpus} \\
        -o ${prefix} \\
        -r ${bed} \\
        ${fasta} \\
        ${bam} > ${log_file}

    mv ${prefix}/3-filt/portcullis_filtered.pass.junctions.bed ${prefix}_portcullis_junctions.bed
    mv ${prefix}/3-filt/portcullis_filtered.pass.junctions.tab ${prefix}_portcullis_junctions.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        portcullis: \$(portcullis --version |& sed '1!d ; s/portcullis //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.spliced.bam
    touch ${prefix}.portcullis.log
    touch ${prefix}.spliced.bam.bai
    touch ${prefix}_portcullis_junctions.bed
    touch ${prefix}_portcullis_junctions.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        portcullis: \$(portcullis --version |& sed '1!d ; s/portcullis //')
    END_VERSIONS
    """
}
