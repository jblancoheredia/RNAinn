process FUSVIZ {
    tag "$meta.patient_id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/fusviz:1.0.2':
        'blancojmskcc/fusviz:1.0.2' }"

    input:
    tuple val(meta) , path(bam), path(bai)
    tuple val(meta2), path(tsv)
    path(gtf)
    val(genome)
    path(cytobands)
    path(protein_domains)

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.patient_id}"
    """
    preFusViz \\
        --sample ${prefix} \\
        --input ${tsv} \\
        --genome ${genome} \\
        --annotations ${gtf} \\
        ${args}

    FusViz \\
        --fusions=${prefix}_FusViz.tsv \\
        --alignments=${bam} \\
        --annotation=${gtf}   \\
        --cytobands=${cytobands} \\
        --output=${prefix}_FusViz.pdf \\
        --transcriptSelection=canonical \\
        --minConfidenceForCircosPlot=High \\
        --proteinDomains=${protein_domains} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusviz: "1.0.2"
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.patient_id}"
    """
    touch ${prefix}_FusViz.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusviz: "1.0.2"
    END_VERSIONS
    """
}
