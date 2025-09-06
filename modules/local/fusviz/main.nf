process FUSVIZ {
    tag "$meta.patient"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/fusviz:3.0.0':
        'blancojmskcc/fusviz:3.0.0' }"

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
    def prefix = task.ext.prefix ?: "${meta.patient}"
    """
    FusViz \\
        --alignments ${bam} \\
        --annotation ${gtf}   \\
        --fusions ${fusions}    \\
        --cytobands ${cytobands}  \\
        --chromosomes ${chromosomes} \\
        --output ${prefix}_FusViz.pdf  \\
        --proteinDomains ${protein_domains} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusviz: "1.0.2"
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.patient}"
    """
    touch ${prefix}_FusViz.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusviz: "1.0.2"
    END_VERSIONS
    """
}
