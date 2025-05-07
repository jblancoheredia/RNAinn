process BAMCUT {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/samtools:1.21' :
        'blancojmskcc/samtools:1.21' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_chunk_*.sort.bam"), emit: bams
    tuple val(meta), path("${meta.id}_chunk_*.sort.bam.bai"), emit: bais
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p temp_splitBAM
    
    samtools sort -n -O SAM ${bam} | \
    awk -v n=1000000 -v FS="\t" 'BEGIN {part=0; line=n} \
    /^@/ {header = header\$0"\\n"; next;} \
    { \
        if(line>=n && \$1!=last_read) { \
            part++; \
            line=1; \
        } \
        print line==1 ? header""\$0 : \$0 | "samtools view -b -o ${prefix}_chunk_"part".bam"; \
        last_read = \$1; \
        line++; \
    }'
    
    for bamfile in ${prefix}_chunk_*.bam; do
        samtools sort -@ ${task.cpus} \${bamfile} -o \${bamfile%.bam}.sort.bam
        samtools index \${bamfile%.bam}.sort.bam
        rm \${bamfile}
    done
    
    rm ${bam} ${bai}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_chunk_1.sort.bam
    touch ${prefix}_chunk_1.sort.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}