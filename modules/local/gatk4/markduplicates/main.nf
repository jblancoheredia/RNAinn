process GATK4_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d9e7bad0f7fbc8f4458d5c3ab7ffaaf0235b59fb:7cc3d06cbf42e28c5e2ebfc7c858654c7340a9d5-0':
        'biocontainers/mulled-v2-d9e7bad0f7fbc8f4458d5c3ab7ffaaf0235b59fb:7cc3d06cbf42e28c5e2ebfc7c858654c7340a9d5-0' }"

    input:
    tuple val(meta) , path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    path "versions.yml",                                    emit: versions
    tuple val(meta), path("*.metrics"),                     emit: metrics
    tuple val(meta), path("*_md.bam"),                      emit: bam_md
    tuple val(meta), path("*_complex_metrics.txt"),         emit: complex_metrics
    tuple val(meta), path("*_deduped.bam"), path("*.bai"),  emit: bam_dp

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def input_list = bam.collect{"--INPUT $it"}.join(' ')
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    // Using samtools and not Markduplicates to compress to CRAM speeds up computation:
    // https://medium.com/@acarroll.dna/looking-at-trade-offs-in-compression-levels-for-genomics-tools-eec2834e8b94
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        MarkDuplicates \\
        $input_list \\
        --OUTPUT ${prefix}_md.bam \\
        --METRICS_FILE ${prefix}.markdup.metrics \\
        --TMP_DIR . \\
        ${reference} \\
        $args

    samtools view \\
        -b \\
        -f 3 \\
        -F 1024 \\
        -F 2048 \\
        -o ${prefix}_deduped.bam \\
        ${prefix}_md.bam

    samtools index \\
        ${prefix}_deduped.bam 

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        EstimateLibraryComplexity \\
        --INPUT ${prefix}_deduped.bam \\
        --OUTPUT ${prefix}_est_lib_complex_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.bam"
    prefix_no_suffix = task.ext.prefix ? prefix.tokenize('.')[0] : "${meta.id}"
    """
    touch ${prefix}_md.bam
    touch ${prefix}.metrics
    touch ${prefix}_deduped.bam
    touch ${prefix}_deduped.bam.bai
    touch ${prefix}_est_lib_complex_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
