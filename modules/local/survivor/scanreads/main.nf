process SURVIVOR_SCAN_READS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/survivor:1.0.7':
        'blancojmskcc/survivor:1.0.7' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val(minrl) // Min read lengt (-1 to disable)

    output:
    tuple val(meta), path("*.error_profile.tsv"), emit: error_profile
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools \\
        view \\
        --threads ${task.cpus} \\
        ${bam} \\
        | \\
        SURVIVOR \\
        scanreads \\
        ${minrl} \\
        ${prefix}.error_profile.tsv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: \$(echo \$(SURVIVOR 2>&1 | grep "Version" | sed 's/^Version: //'))
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.error_profile.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: \$(echo \$(SURVIVOR 2>&1 | grep "Version" | sed 's/^Version: //'))
    END_VERSIONS
    """
}

process SURVIVOR_SCAN_READS_RAW {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/survivor:1.0.7':
        'blancojmskcc/survivor:1.0.7' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val(minrl) // Min read lengt (-1 to disable)

    output:
    tuple val(meta), path("*.raw_error_profile.tsv"), emit: error_profile
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools \\
        view \\
        --threads ${task.cpus} \\
        ${bam} \\
        | \\
        SURVIVOR \\
        scanreads \\
        ${minrl} \\
        ${prefix}.raw_error_profile.tsv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: \$(echo \$(SURVIVOR 2>&1 | grep "Version" | sed 's/^Version: //'))
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.raw_error_profile.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: \$(echo \$(SURVIVOR 2>&1 | grep "Version" | sed 's/^Version: //'))
    END_VERSIONS
    """
}

process SURVIVOR_SCAN_READS_CON {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/survivor:1.0.7':
        'blancojmskcc/survivor:1.0.7' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val(minrl) // Min read lengt (-1 to disable)

    output:
    tuple val(meta), path("*.con_error_profile.tsv"), emit: error_profile
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools \\
        view \\
        --threads ${task.cpus} \\
        ${bam} \\
        | \\
        SURVIVOR \\
        scanreads \\
        ${minrl} \\
        ${prefix}.con_error_profile.tsv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: \$(echo \$(SURVIVOR 2>&1 | grep "Version" | sed 's/^Version: //'))
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.con_error_profile.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: \$(echo \$(SURVIVOR 2>&1 | grep "Version" | sed 's/^Version: //'))
    END_VERSIONS
    """
}
