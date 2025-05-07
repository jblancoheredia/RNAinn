process ALIGN_BAM_ORI {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::fgbio=2.4.0 bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.9 bioconda::star=2.7.10b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/umi_aligner:2.0.0' :
        'blancojmskcc/umi_aligner:2.0.0' }"

    input:
    tuple val(meta) , path(unmapped_bam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    tuple val(meta4), path(star_index)
    tuple val(meta5), path(gtf)
    tuple val(meta6), path(dict)
    val seq_platform
    val seq_center

    output:
    tuple val(meta), path("*_ori.bam")      , emit: bam
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix          = task.ext.prefix ?: "${meta.id}"
    def suffix          = "ori"
    def seq_platform_c  = seq_platform ? "'PL:${seq_platform}'" : ""
    def seq_center_c    = seq_center ? "'CN:${seq_center}'" : ""
    def attrRG          = args.contains("--outSAMattrRGline") ? "" : "--outSAMattrRGline 'ID:${prefix}' ${seq_center_c} 'SM:${prefix}' ${seq_platform_c}"
    """
    samtools fastq -1 ${prefix}_R1_unmapped.fq -2 ${prefix}_R2_unmapped.fq -0 /dev/null -s /dev/null -n -F 0x900 -N -@${task.cpus} ${unmapped_bam}

    STAR \\
        --runThreadN 1 \\
        --genomeDir ${star_index} \\
        --genomeLoad NoSharedMemory \\
        --readFilesIn ${prefix}_R1_unmapped.fq ${prefix}_R2_unmapped.fq \\
        --readFilesCommand cat \\
        --outStd BAM_Unsorted \\
        --outSAMtype BAM Unsorted \\
        --outSAMunmapped Within KeepPairs \\
        --outBAMcompression 0 \\
        --outFilterMultimapNmax 50 \\
        --peOverlapNbasesMin 10 \\
        --peOverlapMMp 0.01 \\
        --alignEndsType Local \\
        --outSAMattributes All \\
        --alignSplicedMateMapLminOverLmate 0.5 \\
        --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --chimSegmentMin 10 \\
        --chimOutType WithinBAM SoftClip \\
        --chimJunctionOverhangMin 10 \\
        --chimScoreDropMax 30 \\
        --chimScoreJunctionNonGTAG 0 \\
        --chimScoreSeparation 1 \\
        --chimSegmentReadGapMax 3 \\
        --chimMultimapNmax 50 \\
        --sjdbGTFfile ${gtf} \\
        --outFileNamePrefix ${prefix}_${suffix}_tmp_ \\
        ${attrRG} >> ${prefix}_${suffix}_tmp.bam

    fgbio -Xmx16g \\
        --compression 0 \\
        --async-io=true \\
        SortBam \\
        --input ${prefix}_${suffix}_tmp.bam \\
        --sort-order TemplateCoordinate \\
        --output ${prefix}_${suffix}_tmp_sort.bam

    fgbio -Xmx16g \\
        --compression 0 \\
        --async-io=true \\
        ZipperBams \\
        --unmapped ${unmapped_bam} \\
        --input ${prefix}_${suffix}_tmp_sort.bam \\
        --ref ${fasta} \\
        --tags-to-reverse Consensus \\
        --tags-to-revcomp Consensus \\
        --output ${prefix}_${suffix}.bam

    rm ${prefix}_${suffix}_tmp_sort.bam ${prefix}_${suffix}_tmp.bam ${prefix}_R1_unmapped.fq ${prefix}_R2_unmapped.fq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = "ori"
    """
    touch ${prefix}_${suffix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process ALIGN_BAM_FIN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::fgbio=2.4.0 bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.9 bioconda::star=2.7.10b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/umi_aligner:2.0.0' :
        'blancojmskcc/umi_aligner:2.0.0' }"

    input:
    tuple val(meta) , path(consensus_bam)       , path(consensus_bai)
    tuple val(meta0), path(simplex_bam)         , path(simplex_bai)
    tuple val(meta1), path(duplex_bam)          , path(duplex_bai)
    tuple val(meta2), path(star_index)
    tuple val(meta3), path(fasta_fai)
    tuple val(meta4), path(fasta)
    tuple val(meta5), path(dict)
    tuple val(meta6), path(gtf)
    val seq_platform
    val seq_center

    output:
    tuple val(meta) , path("*_fin.bam")          , emit: bam
    tuple val(meta) , path("*_fin.duplex.bam")   , emit: duplex_bam
    tuple val(meta) , path("*_fin.simplex.bam")  , emit: simplex_bam
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix          = task.ext.prefix ?: "${meta.id}"
    def suffix          = "fin"
    def seq_platform_c  = seq_platform ? "'PL:${seq_platform}'" : ""
    def seq_center_c    = seq_center ? "'CN:${seq_center}'" : ""
    def attrRG          = args.contains("--outSAMattrRGline") ? "" : "--outSAMattrRGline 'ID:${prefix}' ${seq_center_c} 'SM:${prefix}' ${seq_platform_c}"
    """
    samtools fastq -1 ${prefix}_R1_unmapped.fq -2 ${prefix}_R2_unmapped.fq -0 /dev/null -s /dev/null -n -F 0x900 -N -@${task.cpus} ${consensus_bam}

    STAR \\
        --runThreadN 1 \\
        --genomeDir ${star_index} \\
        --genomeLoad NoSharedMemory \\
        --readFilesIn ${prefix}_R1_unmapped.fq ${prefix}_R2_unmapped.fq \\
        --readFilesCommand cat \\
        --outStd BAM_Unsorted \\
        --outSAMtype BAM Unsorted \\
        --outSAMunmapped Within KeepPairs \\
        --outBAMcompression 0 \\
        --outFilterMultimapNmax 50 \\
        --peOverlapNbasesMin 10 \\
        --peOverlapMMp 0.01 \\
        --alignEndsType Local \\
        --outSAMattributes All \\
        --alignSplicedMateMapLminOverLmate 0.5 \\
        --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --chimSegmentMin 10 \\
        --chimOutType WithinBAM SoftClip \\
        --chimJunctionOverhangMin 10 \\
        --chimScoreDropMax 30 \\
        --chimScoreJunctionNonGTAG 0 \\
        --chimScoreSeparation 1 \\
        --chimSegmentReadGapMax 3 \\
        --chimMultimapNmax 50 \\
        --sjdbGTFfile ${gtf} \\
        --outFileNamePrefix ${prefix}_${suffix}_tmp_ \\
        ${attrRG} >> ${prefix}_${suffix}_tmp.bam

    fgbio -Xmx16g \\
        --compression 0 \\
        --async-io=true \\
        SortBam \\
        --input ${prefix}_${suffix}_tmp.bam \\
        --sort-order TemplateCoordinate \\
        --output ${prefix}_${suffix}_tmp_sort.bam

    fgbio -Xmx16g \\
        --compression 0 \\
        --async-io=true \\
        ZipperBams \\
        --unmapped ${consensus_bam} \\
        --input ${prefix}_${suffix}_tmp_sort.bam \\
        --ref ${fasta} \\
        --tags-to-reverse Consensus \\
        --tags-to-revcomp Consensus \\
        --sort TemplateCoordinate \\
        --output ${prefix}_${suffix}.bam

    rm ${prefix}_${suffix}_tmp_sort.bam ${prefix}_${suffix}_tmp.bam ${prefix}_R1_unmapped.fq ${prefix}_R2_unmapped.fq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = "fin"
    """
    touch ${prefix}_${suffix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}