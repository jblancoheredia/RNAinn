process RAW_READS_RECOVERY {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/rnainn_rawreadsrecovery:1.0.0':
        'blancojmskcc/rnainn_rawreadsrecovery:1.0.0' }"

    input:
    tuple val(meta), path(inputs), path(bam), path(bai)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads1 = inputs[0]
    def reads2 = inputs[1]

    """
    samtools view -F 260 ${bam} | cut -f1 | cut -f1 -d'_' | sort -u > ${prefix}.read_names.txt

    seqtk subseq ${reads1} ${prefix}.read_names.txt | gzip > ${prefix}_R1.fastq.gz
    seqtk subseq ${reads2} ${prefix}.read_names.txt | gzip > ${prefix}_R2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rrr: 1.0.0
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
        seqtk: \$(seqtk 2>&1 | grep Version | sed 's/Version: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """    
    touch ${prefix}_R1.fastq.gz
    touch ${prefix}_R2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rrr: 1.0.0
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
        seqtk: \$(seqtk 2>&1 | grep Version | sed 's/Version: //')
    END_VERSIONS
    """
}
