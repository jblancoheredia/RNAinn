process SAMTOOLS_FASTQ_FIL {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(unmapped_bam)

    output:
    tuple val(meta), path("*_{1,2}.fastq.gz")      , emit: fastq
    path  "versions.yml"                           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = "-1 ${prefix}_consensus_1.fastq.gz -2 ${prefix}_consensus_2.fastq.gz"
    """
    samtools \\
        fastq \\
        --threads ${task.cpus-1} \\
        -0 ${prefix}_other.fastq.gz \\
        $unmapped_bam \\
        $output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = "-1 ${prefix}_consensus_1.fastq.gz -2 ${prefix}_consensus_2.fastq.gz"

    """
    touch ${prefix}_1_consensus.fastq.gz
    touch ${prefix}_2_consensus.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}