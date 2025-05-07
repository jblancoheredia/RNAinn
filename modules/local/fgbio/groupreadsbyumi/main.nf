process FGBIO_GROUPREADSBYUMI {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::fgbio=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.2.1--hdfd78af_0' :
        'quay.io/biocontainers/fgbio:2.2.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(mapped_bam)
    val(strategy)
    val(edits)
    val(include_secondary)
    val(allow_inter_contig)
    val(include_supplementary)
    val(min_map_q)
    val(include_nonpf_reads)
    val(mark_duplicates)

    output:
    tuple val(meta), path("*.grouped.bam")              , emit: bam
    tuple val(meta), path("*.grouped-family-sizes.txt") , emit: histogram
    path "versions.yml"                                 , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_gb = 8
    if (!task.memory) {
        log.info '[fgbio GroupReadsByUmi] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else {
        mem_gb = task.memory.giga
    }

    """
    fgbio \\
        -Xmx${mem_gb}g \\
        --tmp-dir=. \\
        --async-io=true \\
        --compression=1 \\
        GroupReadsByUmi \\
        --threads ${task.cpus} \\
        --strategy ${strategy} \\
        --edits ${edits} \\
        --mark-duplicates ${mark_duplicates} \\
        --include-secondary ${include_secondary} \\
        --include-non-pf-reads ${include_nonpf_reads} \\
        --include-supplementary ${include_supplementary} \\
        --min-map-q ${min_map_q} \\
        --allow-inter-contig ${allow_inter_contig} \\
        --input ${mapped_bam} \\
        --output ${prefix}.grouped.bam \\
        --family-size-histogram ${prefix}.grouped-family-sizes.txt \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.grouped.bam
    touch ${prefix}.grouped-family-sizes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}
