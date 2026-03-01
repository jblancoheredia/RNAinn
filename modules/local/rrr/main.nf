process RAW_READS_RECOVERY {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'docker://community.wave.seqera.io/library/pyranges_pysam_samtools_matplotlib_pruned:80797ab907957d3c':
        'community.wave.seqera.io/library/pyranges_pysam_samtools_matplotlib_pruned:80797ab907957d3c' }"

    input:
    tuple val(meta), path(inputs)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bam    = inputs[2]
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads1 = inputs[0]
    def reads2 = inputs[1]

    """
    rrr \\
        $args \\
        -@ $task.cpus \\
        -o ${prefix}.bam \\
        $bam
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args
    
    touch ${prefix}.bam
    """
}
