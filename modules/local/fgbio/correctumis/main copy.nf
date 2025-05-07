process FGBIO_CORRECTUMIS {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::fgbio=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.2.1--hdfd78af_0' :
        'quay.io/biocontainers/fgbio:2.2.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    val(max_mismatch)
    val(min_distance)
    val(min_corrected)

    output:
    tuple val(meta), path("*.mapped.fixedUMI.bam"), emit: bam
    tuple val(meta), path("*.metrics.txt"),         emit: metrics
    path "versions.yml",                            emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_gb = 8
    if (!task.memory) {
        log.info '[fgbio CorrectUmis] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else {
        mem_gb = task.memory.giga
    }

    """
    fgbio \\
        -Xmx${mem_gb}g \\
        --tmp-dir=. \\
        --async-io=true \\
        --compression=1 \\
        CorrectUmis \\
        -i ${bam} \\
        -o ${prefix}.mapped.fixedUMI.bam \\
        --max-mismatches=${max_mismatch} \\
        --min-distance=${min_distance} \\
        --min-corrected=${min_corrected} \\
        -M ${prefix}.metrics.txt \\
        -r ${prefix}.rejected.bam \\
        -t RX \\
        -u \\
        AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT \\
        CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT \\
        GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT \\
        TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.mapped.fixedUMI.bam
    touch  ${prefix}.metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}
