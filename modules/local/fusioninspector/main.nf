process FUSIONINSPECTOR {
    tag "$meta.patient"
    label 'process_high'

    conda "bioconda::dfam=3.3 bioconda::hmmer=3.3.2 bioconda::star-fusion=1.12.0 bioconda::samtools=1.9 bioconda::star=2.7.8a"
    container 'docker.io/trinityctat/starfusion:1.12.0'

    input:
    tuple val(meta), path(reads), path(fusion_list)
    path reference

    output:
    tuple val(meta), path("*.fusions.abridged.tsv.annotated")     , emit: tsv_annotated
    tuple val(meta), path("*.fusions.abridged.tsv")               , emit: tsv
    tuple val(meta), path("*.coding_effect")       , optional:true, emit: tsv_coding_effect
    tuple val(meta), path("*.igv.*")               , optional:true, emit: igv
    tuple val(meta), path("*.gtf")                 , optional:true, emit: out_gtf
    path "*"                                                      , emit: output
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.patient}"
    def fasta = meta.single_end ? "--left_fq ${reads[0]}" : "--left_fq ${reads[0]} --right_fq ${reads[1]}"
    def args = task.ext.args ?: ''
    """
    cp ${reads[0]} rnaseq_1.fq.gz
    cp ${reads[1]} rnaseq_2.fq.gz

    gunzip rnaseq_1.fq.gz
    gunzip rnaseq_2.fq.gz

    cut -f1 ${fusion_list} | grep -v '^#' | awk -F '--' '\$1 != \$2 && \$1 != \"NA\" && \$2 != \"NA\"' > valid_fusions.txt

    if [[ ! -s valid_fusions.txt ]]; then
        echo "No valid fusions. Skipping." > ${prefix}_FusionInspector.fusions.abridged.tsv
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    fusioninspector: \$(echo \$(FusionInspector --version 2>&1) | sed 's/^.*FusionInspector //; s/version:.*\$//')
END_VERSIONS
        exit 0
    fi

    FusionInspector -O . \\
                    --vis \\
                    --CPU ${task.cpus} \\
                    --annotate \\
                    --out_prefix ${prefix} \\
                    --examine_coding_effect \\
                    --fusions valid_fusions.txt \\
                    --genome_lib ctat_genome_lib \\
                    --left_fq rnaseq_1.fq --right_fq rnaseq_2.fq \\
                    ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioninspector: \$(echo \$(FusionInspector --version 2>&1) | sed 's/^.*FusionInspector //; s/version:.*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.patient}"
    """
    touch ${prefix}.FusionInspector.log
    touch ${prefix}.FusionInspector.fusions.tsv
    touch ${prefix}.FusionInspector.fusions.tsv.annotated.coding_effect
    touch ${prefix}.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioninspector: \$(echo \$(FusionInspector --version 2>&1) | sed 's/^.*FusionInspector //; s/version:.*\$//')
    END_VERSIONS
    """
}