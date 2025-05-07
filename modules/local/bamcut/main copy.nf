process BAMCUT {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/samtools:1.21' :
        'blancojmskcc/samtools:1.21' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta),  path("*.bam"),    emit: bams
    tuple val(meta),  path("*.bai"),    emit: bais
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools idxstats ${bam} | cut -f1 | grep -v '*' > chrom_list.txt
    
    while read chrom; do
        if [[ \${chrom} == "X" || \${chrom} == "Y" ]]; then
            continue
        fi
        samtools view -b ${bam} \${chrom} > ${prefix}_\${chrom}.bam
    done < chrom_list.txt
    
    samtools view -b ${bam} X Y > ${prefix}_XY.bam
    
    rm ${bam}

    rm ${bai}
    
    for bamfile in *.bam; do
        samtools index \${bamfile}
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for bamfile in *.bam do
        touch \${bamfile}
    done

    for baifile in *.bai do
        touch \${baifile}
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}    
