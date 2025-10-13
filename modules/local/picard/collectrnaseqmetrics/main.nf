process PICARD_COLLECTRNASEQMETRICS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::picard=3.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.1.0--hdfd78af_0' :
        'quay.io/biocontainers/picard:3.1.0--hdfd78af_0' }"

    input:
    tuple val(meta) , path(bam), path(bai)
    tuple val(meta2), path(rrna_intervals)
    tuple val(meta3), path(refflat)

    output:
    tuple val(meta), path("*rna_metrics.txt")   , emit: metrics
    path  "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def strandedness = ''
    // def strandedness = '--STRAND_SPECIFICITY FIRST_READ_TRANSCRIPTION_STRAND'
    if ("${meta.strandedness}" == 'forward') {
        strandedness = '--STRAND_SPECIFICITY FIRST_READ_TRANSCRIPTION_STRAND'
    } else if ("${meta.strandedness}" == 'reverse') {
        strandedness = '--STRAND_SPECIFICITY SECOND_READ_TRANSCRIPTION_STRAND'
    } else {
        strandedness = '--STRAND_SPECIFICITY NONE'
    }

    def rrna = rrna_intervals == [] ? '' : "--RIBOSOMAL_INTERVALS ${rrna_intervals}"
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CollectRnaMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    def bam_path = "${bam}"
    def bai_path = "${bai}"
    """
    picard \\
        -Xmx${avail_mem}M \\
        CollectRnaSeqMetrics \\
        --TMP_DIR ./tmp \\
        ${strandedness} \\
        ${rrna} \\
        --REF_FLAT ${refflat} \\
        --INPUT ${bam_path} \\
        --OUTPUT ${prefix}_rna_metrics.txt \\
        --VALIDATION_STRINGENCY SILENT \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectRnaMetrics --version 2>&1 | grep -o 'Version.*' | cut -f2- -d:)
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_rna_metrics.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectRnaMetrics --version 2>&1 | grep -o 'Version.*' | cut -f2- -d:)
    END_VERSIONS
    """
}

process PICARD_COLLECTRNASEQMETRICS_RAW {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::picard=3.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.1.0--hdfd78af_0' :
        'quay.io/biocontainers/picard:3.1.0--hdfd78af_0' }"

    input:
    tuple val(meta) , path(bam)
    tuple val(meta2), path(bai)
    tuple val(meta3), path(refflat)
    tuple val(meta4), path(rrna_intervals)

    output:
    tuple val(meta), path("*.raw.rna_metrics.txt"), emit: metrics
    path  "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def strandedness = ''
    // def strandedness = '--STRAND_SPECIFICITY FIRST_READ_TRANSCRIPTION_STRAND'
    if ("${meta.strandedness}" == 'forward') {
        strandedness = '--STRAND_SPECIFICITY FIRST_READ_TRANSCRIPTION_STRAND'
    } else if ("${meta.strandedness}" == 'reverse') {
        strandedness = '--STRAND_SPECIFICITY SECOND_READ_TRANSCRIPTION_STRAND'
    } else {
        strandedness = '--STRAND_SPECIFICITY NONE'
    }

    def rrna = rrna_intervals == [] ? '' : "--RIBOSOMAL_INTERVALS ${rrna_intervals}"
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CollectRnaMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    def bam_path = "${bam}"
    def bai_path = "${bai}"
    """
    picard \\
        -Xmx${avail_mem}M \\
        CollectRnaSeqMetrics \\
        --TMP_DIR ./tmp \\
        ${strandedness} \\
        ${rrna} \\
        --REF_FLAT ${refflat} \\
        --INPUT ${bam_path} \\
        --OUTPUT ${prefix}.raw.rna_metrics.txt \\
        --VALIDATION_STRINGENCY SILENT \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectRnaMetrics --version 2>&1 | grep -o 'Version.*' | cut -f2- -d:)
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.raw.rna_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectRnaMetrics --version 2>&1 | grep -o 'Version.*' | cut -f2- -d:)
    END_VERSIONS
    """
}

process PICARD_COLLECTRNASEQMETRICS_CON {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::picard=3.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.1.0--hdfd78af_0' :
        'quay.io/biocontainers/picard:3.1.0--hdfd78af_0' }"

    input:
    tuple val(meta) , path(bam), path(bai)
    tuple val(meta2), path(refflat)
    tuple val(meta3), path(rrna_intervals)

    output:
    tuple val(meta), path("*.con.rna_metrics.txt"), emit: metrics
    path  "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def strandedness = ''
    // def strandedness = '--STRAND_SPECIFICITY FIRST_READ_TRANSCRIPTION_STRAND'
    if ("${meta.strandedness}" == 'forward') {
        strandedness = '--STRAND_SPECIFICITY FIRST_READ_TRANSCRIPTION_STRAND'
    } else if ("${meta.strandedness}" == 'reverse') {
        strandedness = '--STRAND_SPECIFICITY SECOND_READ_TRANSCRIPTION_STRAND'
    } else {
        strandedness = '--STRAND_SPECIFICITY NONE'
    }

    def rrna = rrna_intervals == [] ? '' : "--RIBOSOMAL_INTERVALS ${rrna_intervals}"
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CollectRnaMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    def bam_path = "${bam}"
    def bai_path = "${bai}"
    """
    picard \\
        -Xmx${avail_mem}M \\
        CollectRnaSeqMetrics \\
        --TMP_DIR ./tmp \\
        ${strandedness} \\
        ${rrna} \\
        --REF_FLAT ${refflat} \\
        --INPUT ${bam_path} \\
        --OUTPUT ${prefix}.con.rna_metrics.txt \\
        --VALIDATION_STRINGENCY SILENT \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectRnaMetrics --version 2>&1 | grep -o 'Version.*' | cut -f2- -d:)
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.con.rna_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectRnaMetrics --version 2>&1 | grep -o 'Version.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
