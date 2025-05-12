process KALLISTO_QUANT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/kallisto:0.51.1' :
        'blancojmskcc/kallisto:0.51.1' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    path gtf
    path chromosomes
    val fragment_length
    val fragment_length_sd

    output:
    tuple val(meta), path("*_abundance.tsv")  , emit: quant
    tuple val(meta), path("*.run_info.json")  , emit: json_info
    tuple val(meta), path("*.log")            , emit: log
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def gtf_input = gtf ? "--gtf ${gtf}" : ''
    def chromosomes_input = chromosomes ? "--chromosomes ${chromosomes}" : ''
    def strandedness = ''
    if (!args.contains('--fr-stranded') && !args.contains('--rf-stranded')) {
        strandedness =  (meta.strandedness == 'forward') ? '--fr-stranded' :
                        (meta.strandedness == 'reverse') ? '--rf-stranded' : ''
    }

    """
    mkdir -p $prefix && kallisto quant \\
            --threads ${task.cpus} \\
            --index ${index} \\
            ${gtf_input} \\
            ${chromosomes_input} \\
            ${strandedness} \\
            -o $prefix \\
            ${reads} 2> >(tee -a ${prefix}/kallisto_quant.log >&2)

    cp ${prefix}/abundance.tsv ${prefix}_abundance.tsv
    cp ${prefix}/kallisto_quant.log ${prefix}.log
    cp ${prefix}/run_info.json ${prefix}.run_info.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(echo \$(kallisto version) | sed "s/kallisto, version //g" )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p $prefix
    touch ${prefix}.log
    touch ${prefix}.run_info.json
    touch ${prefix}_abundance.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(echo \$(kallisto version) | sed "s/kallisto, version //g" )
    END_VERSIONS
    """
}
