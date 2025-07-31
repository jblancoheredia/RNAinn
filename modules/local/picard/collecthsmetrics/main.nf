process COLLECTHSMETRICS_CON {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.3.0--hdfd78af_0' :
        'quay.io/biocontainers/picard:3.3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    path(bait_intervals, stageAs: "baits/*")
    path(target_intervals, stageAs: 'targets/*')
    val(library)

    output:
    tuple val(meta), path("*.HSmetrics.txt")         , emit: hsmetrics
    tuple val(meta), path("*.target_${library}.covg"), emit: coverage
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CollectHsMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    def bait_interval_list = bait_intervals
    def bait_intervallist_cmd = ""
    if (bait_intervals =~ /.(bed|bed.gz)$/){
        bait_interval_list = bait_intervals.toString().replaceAll(/.(bed|bed.gz)$/, ".interval_list")
        bait_intervallist_cmd = "picard -Xmx${avail_mem}M  BedToIntervalList --INPUT ${bait_intervals} --OUTPUT ${bait_interval_list} --SEQUENCE_DICTIONARY ${dict} --TMP_DIR ."
    }

    def target_interval_list = target_intervals
    def target_intervallist_cmd = ""
    if (target_intervals =~ /.(bed|bed.gz)$/){
        target_interval_list = target_intervals.toString().replaceAll(/.(bed|bed.gz)$/, ".interval_list")
        target_intervallist_cmd = "picard -Xmx${avail_mem}M  BedToIntervalList --INPUT ${target_intervals} --OUTPUT ${target_interval_list} --SEQUENCE_DICTIONARY ${dict} --TMP_DIR ."
    }


    """

    $bait_intervallist_cmd
    $target_intervallist_cmd

    picard \\
        -Xmx${avail_mem}M \\
        CollectHsMetrics \\
        $args \\
        $reference \\
        --BAIT_INTERVALS $bait_interval_list \\
        --TARGET_INTERVALS $target_interval_list \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.con.HSmetrics.txt \\
        --PER_TARGET_COVERAGE ${prefix}.con.target_${library}.covg


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectHsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.con.HSmetrics.txt
    touch ${prefix}.con.target_${library}.covg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectHsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}

process COLLECTHSMETRICS_RAW {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.3.0--hdfd78af_0' :
        'quay.io/biocontainers/picard:3.3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    path(bait_intervals, stageAs: "baits/*")
    path(target_intervals, stageAs: 'targets/*')
    val(library)

    output:
    tuple val(meta), path("*.HSmetrics.txt")         , emit: hsmetrics
    tuple val(meta), path("*.target_${library}.covg"), emit: coverage
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CollectHsMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    def bait_interval_list = bait_intervals
    def bait_intervallist_cmd = ""
    if (bait_intervals =~ /.(bed|bed.gz)$/){
        bait_interval_list = bait_intervals.toString().replaceAll(/.(bed|bed.gz)$/, ".interval_list")
        bait_intervallist_cmd = "picard -Xmx${avail_mem}M  BedToIntervalList --INPUT ${bait_intervals} --OUTPUT ${bait_interval_list} --SEQUENCE_DICTIONARY ${dict} --TMP_DIR ."
    }

    def target_interval_list = target_intervals
    def target_intervallist_cmd = ""
    if (target_intervals =~ /.(bed|bed.gz)$/){
        target_interval_list = target_intervals.toString().replaceAll(/.(bed|bed.gz)$/, ".interval_list")
        target_intervallist_cmd = "picard -Xmx${avail_mem}M  BedToIntervalList --INPUT ${target_intervals} --OUTPUT ${target_interval_list} --SEQUENCE_DICTIONARY ${dict} --TMP_DIR ."
    }


    """

    $bait_intervallist_cmd
    $target_intervallist_cmd

    picard \\
        -Xmx${avail_mem}M \\
        CollectHsMetrics \\
        $args \\
        $reference \\
        --BAIT_INTERVALS $bait_interval_list \\
        --TARGET_INTERVALS $target_interval_list \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.raw.HSmetrics.txt \\
        --PER_TARGET_COVERAGE ${prefix}.raw.target_${library}.covg


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectHsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.raw.HSmetrics.txt
    touch ${prefix}.raw.target_${library}.covg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectHsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
