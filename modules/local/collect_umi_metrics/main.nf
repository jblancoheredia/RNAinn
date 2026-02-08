process COLLECT_UMI_METRICS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/collect_umi_consensus_metrics:1.0.3':
        'blancojmskcc/collect_umi_consensus_metrics:1.0.3' }"

    input:
    tuple val(meta), path(con_bam), path(con_bai)

    output:
    tuple val(meta), path("*.con_umi_counts.txt")   , emit: con_umi_counts
    tuple val(meta), path("*.con_family_sizes.txt") , emit: con_family_sizes
    tuple val(meta), path("*.con_ab_fractions.txt") , emit: con_ab_fractions
    tuple val(meta), path("*.con_yield_metrics.txt"), emit: con_yield_metrics
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=.

    mkdir consensus/

    collect_umi_metrics -i ${con_bam} -o consensus/${prefix}

    mv consensus/${prefix}.umi_counts.txt ${prefix}.con_umi_counts.txt
    mv consensus/${prefix}.family_sizes.txt ${prefix}.con_family_sizes.txt
    mv consensus/${prefix}.duplex_family_sizes.txt ${prefix}.con_ab_fractions.txt
    mv consensus/${prefix}.duplex_yield_metrics.txt ${prefix}.con_yield_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        collect_umi_metrics: 1.0.3
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
	touch ${prefix}.con_umi_counts.txt
	touch ${prefix}.con_family_sizes.txt
	touch ${prefix}.con_ab_fractions.txt
	touch ${prefix}.con_yield_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        collect_umi_metrics: 1.0.3
    END_VERSIONS
    """
}
