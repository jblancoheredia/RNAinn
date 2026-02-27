process FUSVIZ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'docker://community.wave.seqera.io/library/pyranges_pysam_samtools_matplotlib_pruned:80797ab907957d3c':
        'community.wave.seqera.io/library/pyranges_pysam_samtools_matplotlib_pruned:80797ab907957d3c' }"

    input:
    tuple val(meta) ,path(bam), path(bai), path(fusions)
    path(protein_domains)
    path(chromosomes)
    path(cytobands)
    path(gtf)

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.png"), emit: png
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=\$(mktemp -d)

    export MPLBACKEND=Agg

    if [ \$(wc -l < ${fusions}) -gt 1 ]; then
        FusViz.py \\
            --alignments=${bam} \\
            --annotation=${gtf}   \\
            --fusions=${fusions}    \\
            --cytobands=${cytobands}  \\
            --chromosomes=${chromosomes} \\
            --proteinDomains ${protein_domains} \\
            ${args}
    else
        echo \"No fusions detected. Creating empty report.\" > ${prefix}_FusViz_no_fusions.txt
        touch ${prefix}_FusViz.pdf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusviz: \$(echo \$(FusViz.py --version 2>&1) | sed 's/^.*FusViz //;  s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_FusViz.pdf
    touch ${prefix}_FusViz.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusviz: \$(echo \$(FusViz.py --version 2>&1) | sed 's/^.*FusViz //;  s/ .*\$//')
    END_VERSIONS
    """
}
