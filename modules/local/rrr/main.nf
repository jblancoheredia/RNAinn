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
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    // TODO nf-core: If the module doesn't use arguments ($args), you SHOULD remove:
    //               - The definition of args `def args = task.ext.args ?: ''` above.
    //               - The use of the variable in the script `echo $args ` below.
    """
    echo $args
    
    touch ${prefix}.bam
    """
}
