process COLLECT_UMI_METRICS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/collect_umi_consensus_metrics:1.0.3':
        'blancojmskcc/collect_umi_consensus_metrics:1.0.3' }"

    input:
    tuple val(meta), path(cnsns_bam), path{cnsns_bai}

    output:
    tuple val(meta), path("*.cons_umi_counts.txt")        , emit: cons_umi_counts
    tuple val(meta), path("*.cons_family_sizes.txt")      , emit: cons_family_sizes
    tuple val(meta), path("*.consensus_umi_counts.txt")   , emit: consensus_umi_counts
    tuple val(meta), path("*.consensus_family_sizes.txt") , emit: consensus_family_sizes
    tuple val(meta), path("*.consensus_yield_metrics.txt"), emit: consensus_yield_metrics
    path  "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=.

    mkdir consensus/

    collect_umi_metrics -i ${cnsns_bam} -o consensus/${prefix}

    mv consensus/${prefix}.umi_counts.txt ${prefix}.cons_umi_counts.txt
    mv consensus/${prefix}.family_sizes.txt ${prefix}.cons_family_sizes.txt
    mv consensus/${prefix}.duplex_umi_counts.txt ${prefix}.consensus_umi_counts.txt
    mv consensus/${prefix}.duplex_family_sizes.txt ${prefix}.consensus_family_sizes.txt
    mv consensus/${prefix}.duplex_yield_metrics.txt ${prefix}.consensus_yield_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        collect_umi_metrics: 1.0.3
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
	touch ${prefix}.cons_umi_counts.txt
	touch ${prefix}.cons_family_sizes.txt
	touch ${prefix}.consensus_umi_counts.txt
	touch ${prefix}.consensus_family_sizes.txt
	touch ${prefix}.consensus_yield_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        collect_umi_metrics: 1.0.3
    END_VERSIONS
    """
}
