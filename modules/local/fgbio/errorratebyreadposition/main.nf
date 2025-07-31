process FGBIO_ERRORRATEBYREADPOSITION {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::fgbio=2.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.2.1--hdfd78af_0' :
        'quay.io/biocontainers/fgbio:2.2.1--hdfd78af_0' }"

    input:
    tuple val(meta),  path(bam), path(bai)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fasta_fai)
    tuple val(meta3), path(fasta_dict)
    path(dbSNP)
    path(dbSNP_idx)
    path(intervals)

    output:
    tuple val(meta), path("*.error_rate_by_read_position.txt"), emit: txt
    tuple val(meta), path("*.error_rate_by_read_position.pdf"), emit: pdf
    path "versions.yml",                                        emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_gb = 8
    if (!task.memory) {
        log.info '[fgbio ErrorRateByReadPosition] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else {
        mem_gb = task.memory.giga
    }

    """
    fgbio \\
        -Xmx${mem_gb}g \\
        --tmp-dir=. \\
        --async-io=true \\
        --compression=1 \\
        ErrorRateByReadPosition \\
        -i ${bam} \\
        -o ${prefix} \\
        -r ${fasta} \\
        -v ${dbSNP} \\
        -l ${intervals} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.error_rate_by_read_position.txt
    touch  ${prefix}.error_rate_by_read_position.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}

process FGBIO_ERRORRATEBYREADPOSITION_RAW {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::fgbio=2.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.2.1--hdfd78af_0' :
        'quay.io/biocontainers/fgbio:2.2.1--hdfd78af_0' }"

    input:
    tuple val(meta),  path(bam)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fasta_fai)
    tuple val(meta3), path(fasta_dict)
    path(dbSNP)
    path(dbSNP_idx)
    path(intervals)

    output:
    tuple val(meta), path("*.error_rate_by_read_position.txt"), emit: txt
    tuple val(meta), path("*.error_rate_by_read_position.pdf"), emit: pdf
    path "versions.yml",                                        emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_gb = 8
    if (!task.memory) {
        log.info '[fgbio ErrorRateByReadPosition] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else {
        mem_gb = task.memory.giga
    }

    """
    fgbio \\
        -Xmx${mem_gb}g \\
        --tmp-dir=. \\
        --async-io=true \\
        --compression=1 \\
        ErrorRateByReadPosition \\
        -i ${bam} \\
        -o ${prefix}_raw \\
        -r ${fasta} \\
        -v ${dbSNP} \\
        -l ${intervals} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}_raw.error_rate_by_read_position.txt
    touch  ${prefix}_raw.error_rate_by_read_position.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}

process FGBIO_ERRORRATEBYREADPOSITION_CON {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::fgbio=2.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.2.1--hdfd78af_0' :
        'quay.io/biocontainers/fgbio:2.2.1--hdfd78af_0' }"

    input:
    tuple val(meta),  path(bam), path(bai)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fasta_fai)
    tuple val(meta3), path(fasta_dict)
    path(dbSNP)
    path(dbSNP_idx)
    path(intervals)

    output:
    tuple val(meta), path("*.error_rate_by_read_position.txt"), emit: txt
    tuple val(meta), path("*.error_rate_by_read_position.pdf"), emit: pdf
    path "versions.yml",                                        emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_gb = 8
    if (!task.memory) {
        log.info '[fgbio ErrorRateByReadPosition] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else {
        mem_gb = task.memory.giga
    }

    """
    fgbio \\
        -Xmx${mem_gb}g \\
        --tmp-dir=. \\
        --async-io=true \\
        --compression=1 \\
        ErrorRateByReadPosition \\
        -i ${bam} \\
        -o ${prefix}_fin \\
        -r ${fasta} \\
        -v ${dbSNP} \\
        -l ${intervals} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}_fin.error_rate_by_read_position.txt
    touch  ${prefix}_fin.error_rate_by_read_position.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}
