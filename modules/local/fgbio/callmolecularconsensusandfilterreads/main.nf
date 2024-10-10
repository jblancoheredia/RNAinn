process FGBIO_CALLMOLECULARCONSENSUSANDFILTERREADS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.2.1--hdfd78af_0' :
        'biocontainers/fgbio:2.2.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(grouped_bam)
    path fasta
    val call_min_reads
    val call_min_baseq
    val filter_min_reads
    val filter_min_baseq

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_consensus_filtered_unmapped"
    def mem_gb = 8
    if (!task.memory) {
        log.info '[fgbio CallMolecularConsensusReads] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else {
        mem_gb = task.memory.giga
    }
    if ("$grouped_bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    fgbio \\
        -Xmx${mem_gb}g \\
        --tmp-dir=. \\
        --async-io=true \\
        --compression=1 \\
        CallMolecularConsensusReads \\
        --input $grouped_bam \\
        --output /dev/stdout \\
        --min-reads ${call_min_reads} \\
        --min-input-base-quality ${call_min_baseq} \\
        --threads ${task.cpus} \\
        |  fgbio -Xmx8g --compression 1 FilterConsensusReads \\
            --input /dev/stdin \\
            --output ${prefix}.bam \\
            --ref ${fasta} \\
            --min-reads ${filter_min_reads} \\
            --min-base-quality ${filter_min_baseq} \\
            --max-base-error-rate 0.2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_consensus_filtered_unmapped"
    if ("$grouped_bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """

}
