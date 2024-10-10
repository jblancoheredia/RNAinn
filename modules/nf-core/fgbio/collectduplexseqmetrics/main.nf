process FGBIO_COLLECTDUPLEXSEQMETRICS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-51891ad0b60843e4aade9cde2eb5d40c5ae92b80:72c944cdea5caff7f03b96034968ce2a4f1737bc-0':
        'biocontainers/mulled-v2-51891ad0b60843e4aade9cde2eb5d40c5ae92b80:72c944cdea5caff7f03b96034968ce2a4f1737bc-0' }"

    input:
    tuple val(meta), path(grouped_bam)
    path interval_list

    output:
    output:
    tuple val(meta), path("*duplex_seq_metrics*.txt"), emit: metrics
    tuple val(meta), path("*duplex_seq_metrics*.pdf"), emit: pdf
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_gb = 8
    if (!task.memory) {
        log.info '[fgbio CollectDuplexSeqMetrics] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else {
        mem_gb = task.memory.giga
    }
    """
    fgbio \\
        -Xmx${mem_gb}g \\
        --tmp-dir=. \\
        --async-io=true \\
        --compression=1 \\
        CollectDuplexSeqMetrics \\
        --input $grouped_bam \\
        --output ${prefix}.duplex_seq_metrics \\
        --duplex-umi-counts=true \\
        $interval_list \\
        $args;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.duplex_seq_metrics.duplex_family_sizes.txt
    touch ${prefix}.duplex_seq_metrics.duplex_umi_counts.txt
    touch ${prefix}.duplex_seq_metrics.duplex_yield_metrics.txt
    touch ${prefix}.duplex_seq_metrics.family_sizes.txt
    touch ${prefix}.duplex_seq_metrics.umi_counts.txt
    touch ${prefix}.duplex_seq_metrics.duplex_qc.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}
