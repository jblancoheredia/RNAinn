process FUSVIZ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/rnainn_fusviz:5.0.0':
        'blancojmskcc/rnainn_fusviz:5.0.0' }"

    input:
    tuple val(meta) ,path(bam), path(bai), path(fusions)
    path(protein_domains)
    path(chromosomes)
    path(cytobands)
    val(genome)
    path(gtf)

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=\$(mktemp -d)

    if [ \$(wc -l < ${fusions}) -gt 1 ]; then
        FusViz \\
            --alignments ${bam} \\
            --annotation ${gtf}   \\
            --fusions ${fusions}    \\
            --cytobands ${cytobands}  \\
            --chromosomes ${chromosomes} \\
            --output ${prefix}_FusViz.pdf  \\
            --proteinDomains ${protein_domains} \\
            ${args}
    else
        echo \"No fusions detected. Creating empty report.\" > ${prefix}_FusViz_no_fusions.txt
        touch ${prefix}_FusViz.pdf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusviz: "5.0.0"
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_FusViz.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusviz: "5.0.0"
    END_VERSIONS
    """
}
