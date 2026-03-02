process FGBIO_GROUPREADSBYUMI {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::fgbio=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/fgbio:2.4.0' :
        'blancojmskcc/fgbio:2.4.0' }"

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
    tuple val(meta), path("*.deduped.bam"), path("*.deduped.bam.bai"), emit: bam_bai_deduped
    tuple val(meta), path("*.grouped-family-sizes.txt")              , emit: histogram
    tuple val(meta), path("*.grouped.bam")                           , emit: bam
    path "versions.yml"                                              , emit: versions

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

    samtools view \\
        -b \\
        -f 3 \\
        -F 1024 \\
        -F 2048 \\
        -o ${prefix}.tmp.bam \\
        ${prefix}.grouped.bam

    samtools sort \\
        -@ 8 \\
        -o ${prefix}.deduped.bam \\
        ${prefix}.tmp.bam

    samtools index \\
        -@ 8 \\
        ${prefix}.deduped.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.grouped.bam
    touch ${prefix}.deduped.bam
    touch ${prefix}.deduped.bam.bai
    touch ${prefix}.grouped-family-sizes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
