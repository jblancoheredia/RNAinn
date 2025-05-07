process DOWNSAMPLING_COUNT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.4--h577a1d6_3':
        'quay.io/biocontainers/seqtk:1.4--h577a1d6_3' }"

    input:
    tuple val(meta), path(reads)

    output:
    path "${meta.id}_total_reads.txt"   , emit: total_reads
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    total_reads=0
    for fastq in ${reads}; do
        total_reads=\$((\$total_reads + \$(zcat \$fastq | wc -l) / 4))
    done

    rounded_reads=\$(awk -v reads=\$total_reads 'BEGIN {printf "%d", int(reads/1000000)*1000000}')

    echo "\$rounded_reads" > ${prefix}_total_reads.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        downsampling_count: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    total_reads=0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        downsampling_count: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
