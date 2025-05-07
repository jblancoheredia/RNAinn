process DOWNSAMPLING_SEQTK {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.4--h577a1d6_3':
        'quay.io/biocontainers/seqtk:1.4--h577a1d6_3' }"

    input:
    tuple val(meta) , path(reads)
    val(min_reads)

    output:
    tuple val(meta), path("*_downsampled.fastq.gz"), emit: downsampled_reads
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def final_reads = args.contains('final_reads=') ? args.split('final_reads=')[1].split()[0] : min_reads
    """
    echo "Final totalreads: ${final_reads}"
        
    downsample_pairs=\$((${final_reads} / 2))

    echo "Final read pairs: \${downsample_pairs}"
   
    for fastq in ${reads}; do
        base_name=\$(basename \$fastq .fastq.gz)
        seqtk sample -s212 \$fastq \$downsample_pairs | gzip > \${base_name}_downsampled.fastq.gz
    done
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        downsampling_seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for fastq in *_downsampled.fastq.gz; do \
        touch \${fastq} \
    done
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        downsampling_seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}