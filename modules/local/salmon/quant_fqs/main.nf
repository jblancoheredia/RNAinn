process SALMON_QUANT_FQS {
    tag "$meta.id"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.10.1--h7e5ed60_0' :
        'biocontainers/salmon:1.10.1--h7e5ed60_0' }"

    input:
    tuple val(meta), path(reads)
    path  index
    val   lib_type

    output:
    tuple val(meta), path("${prefix}")              , emit: results
    tuple val(meta), path("*info.json")             , emit: json_info, optional: true
    tuple val(meta), path("*lib_format_counts.json"), emit: lib_format_counts, optional: true
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"

    def reads1 = [], reads2 = []
    meta.single_end ? [reads].flatten().each{reads1 << it} : reads.eachWithIndex{ v, ix -> ( ix & 1 ? reads2 : reads1) << v }
    """
    salmon quant \\
        -i $index \\
        -o $prefix \\
        -l $lib_type \\
        -p ${task.cpus} \\
        --validateMappings \\
        -1 ${reads1.join(",")} -2 ${reads2.join(",")} \\

    if [ -f $prefix/aux_info/meta_info.json ]; then
        cp $prefix/aux_info/meta_info.json "${prefix}_meta_info.json"
    fi
    if [ -f $prefix/lib_format_counts.json ]; then
        cp $prefix/lib_format_counts.json "${prefix}_lib_format_counts.json"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}_meta_info.json
    touch ${prefix}_lib_format_counts.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
