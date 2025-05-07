process FILTER_CONTIGS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::biopython=1.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'quay.io/biocontainers/biopython:1.81' }"

    input:
    tuple val(meta), path(contigs)
    val min_length
    val min_coverage

    output:
    tuple val(meta), path("*_filtered.fa"), emit: fasta
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.0' // Version of your script
    """
    gunzip $contigs

    filter_contigs.py \\
        ${prefix}.contigs.fa \\
        ${prefix}.contigs_filtered.fa \\
        -l $min_length \\
        -c $min_coverage

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        filter_contigs: $VERSION
    END_VERSIONS
    """
}