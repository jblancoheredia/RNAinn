process FGBIO_SORT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.2.1--hdfd78af_0' :
        'quay.io/biocontainers/fgbio:2.2.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    val suffix

    output:
    tuple val(meta), path("*mapped_sorted.bam")    , optional:true, emit: bam
    tuple val(meta), path("*filtered_sorted.bam")  , optional:true, emit: fil
    tuple val(meta), path("*unmapped_sorted.bam")  , optional:true, emit: ubam
    tuple val(meta), path("*raw_mapped_sorted.bam"), optional:true, emit: raw
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def suffix = task.ext.suffix ?: "${suffix}"
    def prefix = task.ext.prefix ?: "${meta.id}"

    def mem_gb = 8
    if (!task.memory) {
        log.info '[fgbio SortBam] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else if (mem_gb > task.memory.giga) {
        if (task.memory.giga < 2) {
            mem_gb = 1
        } else {
            mem_gb = task.memory.giga - 1
        }
    }

    """
    fgbio -Xmx${mem_gb}g \\
        --async-io=true \\
        --tmp-dir=. \\
        --compression=1 \\
        SortBam \\
        -i $bam \\
        $args \\
        -o ${prefix}_${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "${suffix}"

    """
    touch ${prefix}_${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}
