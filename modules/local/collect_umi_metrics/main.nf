process COLLECT_UMI_METRICS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/collect_umi_consensus_metrics:1.0.3':
        'blancojmskcc/collect_umi_consensus_metrics:1.0.3' }"

    input:
    tuple val(meta), path(con_bam), path{con_bai},
                     path(dup_bam), path{dup_bai},
                     path(sim_bam), path{sim_bai}

    output:
    tuple val(meta), path("*.con_umi_counts.txt")   , emit: con_umi_counts
    tuple val(meta), path("*.dup_umi_counts.txt")   , emit: dup_umi_counts
    tuple val(meta), path("*.sim_umi_counts.txt")   , emit: sim_umi_counts
    tuple val(meta), path("*.con_family_sizes.txt") , emit: con_family_sizes
    tuple val(meta), path("*.dup_family_sizes.txt") , emit: dup_family_sizes
    tuple val(meta), path("*.sim_family_sizes.txt") , emit: sim_family_sizes
    tuple val(meta), path("*.con_ab_fractions.txt") , emit: con_ab_fractions
    tuple val(meta), path("*.dup_ab_fractions.txt") , emit: dup_ab_fractions
    tuple val(meta), path("*.sim_ab_fractions.txt") , emit: sim_ab_fractions
    tuple val(meta), path("*.con_yield_metrics.txt"), emit: con_yield_metrics
    tuple val(meta), path("*.dup_yield_metrics.txt"), emit: dup_yield_metrics
    tuple val(meta), path("*.sim_yield_metrics.txt"), emit: sim_yield_metrics
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=.

    mkdir duplex/
    mkdir simplex/
    mkdir consensus/

    collect_umi_metrics -i ${dup_bam} -o duplex/${prefix}
    collect_umi_metrics -i ${sim_bam} -o simplex/${prefix}
    collect_umi_metrics -i ${con_bam} -o consensus/${prefix}

    mv duplex/${prefix}.umi_counts.txt ${prefix}.dup_umi_counts.txt
    mv duplex/${prefix}.family_sizes.txt ${prefix}.dup_family_sizes.txt
    mv duplex/${prefix}.duplex_family_sizes.txt ${prefix}.dup_ab_fractions.txt
    mv duplex/${prefix}.duplex_yield_metrics.txt ${prefix}.dup_yield_metrics.txt

    mv simplex/${prefix}.umi_counts.txt ${prefix}.sim_umi_counts.txt
    mv simplex/${prefix}.family_sizes.txt ${prefix}.sim_family_sizes.txt
    mv simplex/${prefix}.duplex_family_sizes.txt ${prefix}.sim_ab_fractions.txt
    mv simplex/${prefix}.duplex_yield_metrics.txt ${prefix}.sim_yield_metrics.txt

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
	touch ${prefix}.dup_umi_counts.txt
	touch ${prefix}.dup_family_sizes.txt
	touch ${prefix}.dup_ab_fractions.txt
	touch ${prefix}.dup_yield_metrics.txt
	touch ${prefix}.sim_umi_counts.txt
	touch ${prefix}.sim_family_sizes.txt
	touch ${prefix}.sim_ab_fractions.txt
	touch ${prefix}.sim_yield_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        collect_umi_metrics: 1.0.3
    END_VERSIONS
    """
}
