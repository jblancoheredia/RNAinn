process ALIGN_BAM_RAW {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::fgbio=2.4.0 bioconda::samtools=1.9 bioconda::star=2.7.10b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/umi_aligner:2.0.0' :
        'blancojmskcc/umi_aligner:2.0.0' }"

    input:
    tuple val(meta),  path(unmapped_bam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    tuple val(meta4), path(dict)
    tuple val(meta5), path(index)
    tuple val(meta1), path(gtf)
    val seq_platform
    val seq_center
    val sort_type

    output:
    tuple val(meta), path("*.mapped.bam")       , emit: bam
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def samtools_fastq_args = task.ext.samtools_fastq_args ?: ''
    def samtools_sort_args = task.ext.samtools_sort_args ?: ''
    def star_args = task.ext.star_args ?: "--readFilesCommand zcat --peOverlapNbasesMin 10 --peOverlapMMp 0.01 —-outSAMunmapped Within KeepPairs --outBAMcompression 0 --outFilterMultimapNmax 50 --alignEndsType Local --chimOutType WithinBAM SoftClip --outSAMattributes All --sjdbGTFfile $gtf"
    def fgbio_args = task.ext.fgbio_args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = "mapped" 
    def fgbio_mem_gb = 4
    def extra_command = ""

    if (!task.memory) {
        log.info '[fgbio ZipperBams] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else if (fgbio_mem_gb > task.memory.giga) {
        if (task.memory.giga < 2) {
            fgbio_mem_gb = 1
        } else {
            fgbio_mem_gb = task.memory.giga - 1
        }
    }
    def seq_platform    = seq_platform ? "'PL:$seq_platform'" : ""
    def seq_center      = seq_center ? "'CN:$seq_center'" : ""
    def attrRG          = args.contains("--outSAMattrRGline") ? "" : "--outSAMattrRGline 'ID:$prefix' $seq_center 'SM:$prefix' $seq_platform"
    fgbio_zipper_bams_output = "/dev/stdout"
    fgbio_zipper_bams_compression = 0
    extra_command = " | samtools sort -n"
    extra_command += " --threads "+ task.cpus
    extra_command += " -o " + prefix + "." + suffix + ".bam "
    extra_command += " -"

    """
    samtools \\
        fastq \\
        --threads ${task.cpus-1}    \\
        -0 ${prefix}_other.fastq.gz \\
        -1 ${prefix}_ori_1.fastq.gz \\
        -2 ${prefix}_ori_2.fastq.gz \\
        $unmapped_bam
        
    STAR \\
        --genomeDir $index \\
        --readFilesIn ${prefix}_ori_1.fastq.gz ${prefix}_ori_2.fastq.gz \\
        --runThreadN 1 \\
        --outFileNamePrefix ${prefix}_${suffix}_ \\
        --outSAMtype BAM Unsorted \\
        --outSAMmultNmax 1 \\
        $attrRG \\
        $star_args

    fgbio \\
        -Xmx${fgbio_mem_gb}g \\
        --compression ${fgbio_zipper_bams_compression} \\
        --async-io=true \\
        ZipperBams \\
        --unmapped ${unmapped_bam} \\
        --input ${prefix}_${suffix}_Aligned.out.bam \\
        --ref ${fasta} \\
        --output ${fgbio_zipper_bams_output} \\
        ${extra_command};

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def index_command = (sort_type != "template-coordinate") ? "touch ${prefix}.mapped.bam.bai" : ""
    """
    touch ${prefix}.mapped.bam
    ${index_command}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}

process ALIGN_BAM_CON {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::fgbio=2.4.0 bioconda::samtools=1.9 bioconda::star=2.7.10b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/umi_aligner:2.0.0' :
        'blancojmskcc/umi_aligner:2.0.0' }"

    input:
    tuple val(meta) , path(bam), path(bai), path(duplex_bam), path(duplex_bai), path(simplex_bam), path(simplex_bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    tuple val(meta4), path(dict)
    tuple val(meta5), path(index)
    tuple val(meta6), path(gtf)
    val seq_platform
    val seq_center
    val sort_type

    output:

    tuple val(meta), path("*.mapped.bam")        , path("*.mapped.bai")        , emit: con_bam_bai
    tuple val(meta), path("*.mapped_duplex.bam") , path("*.mapped_duplex.bai") , emit: dup_bam_bai
    tuple val(meta), path("*.mapped_simplex.bam"), path("*.mapped_simplex.bam"), emit: sim_bam_bai
    path "versions.yml"                                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def samtools_fastq_args = task.ext.samtools_fastq_args ?: ''
    def samtools_sort_args = task.ext.samtools_sort_args ?: ''
    def star_args = task.ext.star_args ?: "--readFilesCommand zcat --peOverlapNbasesMin 10 --peOverlapMMp 0.01 --outSAMunmapped Within KeepPairs --outBAMcompression 0 --outFilterMultimapNmax 50 --alignEndsType Local --chimOutType WithinBAM SoftClip --outSAMattributes All --sjdbGTFfile $gtf"
    def fgbio_args = task.ext.fgbio_args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fgbio_mem_gb = 4
    if (!task.memory) {
        log.info '[fgbio ZipperBams] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else if (fgbio_mem_gb > task.memory.giga) {
        if (task.memory.giga < 2) {
            fgbio_mem_gb = 1
        } else {
            fgbio_mem_gb = task.memory.giga - 1
        }
    }
    def seq_platform    = seq_platform ? "'PL:$seq_platform'" : ""
    def seq_center      = seq_center ? "'CN:$seq_center'" : ""
    def attrRG          = args.contains("--outSAMattrRGline") ? "" : "--outSAMattrRGline 'ID:$prefix' $seq_center 'SM:$prefix' $seq_platform"
    fgbio_zipper_bams_output = "/dev/stdout"
    fgbio_zipper_bams_compression = 0
    """
    samtools \\
        fastq \\
        --threads ${task.cpus-1}    \\
        -0 ${prefix}_other.fastq.gz \\
        -1 ${prefix}_con_1.fastq.gz \\
        -2 ${prefix}_con_2.fastq.gz \\
        ${bam}

    STAR \\
        --genomeDir ${index} \\
        --readFilesIn ${prefix}_con_1.fastq.gz ${prefix}_con_2.fastq.gz \\
        --runThreadN 1 \\
        --outFileNamePrefix ${prefix}_mapped_ \\
        --outSAMtype BAM Unsorted \\
        --outSAMmultNmax 1 \\
        ${attrRG} \\
        ${star_args}

    samtools \\
        sort \\
        -n \\
        --threads ${task.cpus-1} \\
        -o ${prefix}_mapped_querysort0.bam \\
        ${bam}

    samtools \\
        sort \\
        -n \\
        --threads ${task.cpus-1} \\
        -o ${prefix}_mapped_querysort1.bam \\
        ${prefix}_mapped_Aligned.out.bam

    fgbio \\
        -Xmx${fgbio_mem_gb}g \\
        --compression ${fgbio_zipper_bams_compression} \\
        --async-io=true \\
        ZipperBams \\
        --unmapped ${prefix}_mapped_querysort0.bam \\
        --input ${prefix}_mapped_querysort1.bam \\
        --ref ${fasta} \\
        --output ${fgbio_zipper_bams_output} \\
        --tags-to-reverse Consensus \\
        --tags-to-revcomp Consensus \\
        | samtools sort \\
        -@ ${task.cpus} \\
        -o ${prefix}.mapped.bam -

    samtools \\
        index \\
        -@ ${task.cpus} \\
        ${prefix}.mapped.bam \\
        ${prefix}.mapped.bai

    samtools \\
        fastq \\
        --threads ${task.cpus-1}    \\
        -0 ${prefix}_other.fastq.gz \\
        -1 ${prefix}_dup_1.fastq.gz \\
        -2 ${prefix}_dup_2.fastq.gz \\
        ${duplex_bam}

    STAR \\
        --genomeDir ${index} \\
        --readFilesIn ${prefix}_dup_1.fastq.gz ${prefix}_dup_2.fastq.gz \\
        --runThreadN 1 \\
        --outFileNamePrefix ${prefix}_mapped_duplex_ \\
        --outSAMtype BAM Unsorted \\
        --outSAMmultNmax 1 \\
        ${attrRG} \\
        ${star_args}

    samtools \\
        sort \\
        -n \\
        --threads ${task.cpus-1} \\
        -o ${prefix}_mapped_duplex_querysort0.bam \\
        ${duplex_bam}

    samtools \\
        sort \\
        -n \\
        --threads ${task.cpus-1} \\
        -o ${prefix}_mapped_duplex_querysort1.bam \\
        ${prefix}_mapped_duplex_Aligned.out.bam

    fgbio \\
        -Xmx${fgbio_mem_gb}g \\
        --compression ${fgbio_zipper_bams_compression} \\
        --async-io=true \\
        ZipperBams \\
        --unmapped ${prefix}_mapped_duplex_querysort0.bam \\
        --input ${prefix}_mapped_duplex_querysort1.bam \\
        --ref ${fasta} \\
        --output ${fgbio_zipper_bams_output} \\
        --tags-to-reverse Consensus \\
        --tags-to-revcomp Consensus \\
        | samtools sort \\
        -@ ${task.cpus} \\
        -o ${prefix}.mapped_duplex.bam -

    samtools \\
        index \\
        -@ ${task.cpus} \\
        ${prefix}.mapped_duplex.bam \\
        ${prefix}.mapped_duplex.bai

    samtools \\
        fastq \\
        --threads ${task.cpus-1}    \\
        -0 ${prefix}_other.fastq.gz \\
        -1 ${prefix}_sim_1.fastq.gz \\
        -2 ${prefix}_sim_2.fastq.gz \\
        ${simplex_bam}

    STAR \\
        --genomeDir ${index} \\
        --readFilesIn ${prefix}_sim_1.fastq.gz ${prefix}_sim_2.fastq.gz \\
        --runThreadN 1 \\
        --outFileNamePrefix ${prefix}_mapped_simplex_ \\
        --outSAMtype BAM Unsorted \\
        --outSAMmultNmax 1 \\
        ${attrRG} \\
        ${star_args}

    samtools \\
        sort \\
        -n \\
        --threads ${task.cpus-1} \\
        -o ${prefix}_mapped_simplex_querysort0.bam \\
        ${simplex_bam}

    samtools \\
        sort \\
        -n \\
        --threads ${task.cpus-1} \\
        -o ${prefix}_mapped_simplex_querysort1.bam \\
        ${prefix}_mapped_simplex_Aligned.out.bam

    fgbio \\
        -Xmx${fgbio_mem_gb}g \\
        --compression ${fgbio_zipper_bams_compression} \\
        --async-io=true \\
        ZipperBams \\
        --unmapped ${prefix}_mapped_simplex_querysort0.bam \\
        --input ${prefix}_mapped_simplex_querysort1.bam \\
        --ref ${fasta} \\
        --output ${fgbio_zipper_bams_output} \\
        --tags-to-reverse Consensus \\
        --tags-to-revcomp Consensus \\
        | samtools sort \\
        -@ ${task.cpus} \\
        -o ${prefix}.mapped_simplex.bam -

    samtools \\
        index \\
        -@ ${task.cpus} \\
        ${prefix}.mapped_simplex.bam \\
        ${prefix}.mapped_simplex.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | tail -n 1)
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mapped.bam
    touch ${prefix}.mapped.bai
    touch ${prefix}.mapped_duplex.bam
    touch ${prefix}.mapped_duplex.bai
    touch ${prefix}.mapped_simplex.bam
    touch ${prefix}.mapped_simplex.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | tail -n 1)
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}
