process FASTP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--h5f740d0_0' :
        'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*_filtered.fastq.gz') , emit: reads
    tuple val(meta), path('*.json')           , optional:true, emit: json
    tuple val(meta), path('*.html')           , optional:true, emit: html
    tuple val(meta), path('*.log')            , optional:true, emit: log
    tuple val(meta), path('*.fail.fastq.gz')  , optional:true, emit: reads_fail
    tuple val(meta), path('*.merged.fastq.gz'), optional:true, emit: reads_merged
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """ 
    fastp \\
        -i ${prefix}_R1.fastq.gz \\
        -I ${prefix}_R2.fastq.gz \\
        -o ${prefix}_R1_filtered.fastq.gz \\
        -O ${prefix}_R2_filtered.fastq.gz \\
        -l 30

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
    stub:
    def prefix              = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_R1_filtered.fastq.gz
    touch ${prefix}_R2_filtered.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}
