process MOSDEPTH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mosdepth:0.3.8--hd299d5a_0' :
        'quay.io/biocontainers/mosdepth:0.3.8--hd299d5a_0'}"

    input:
    tuple val(meta),  path(bam)
    tuple val(meta1), path(bai)
    tuple val(meta2), path(fasta)
    path(fai)
    path(bed)
    path(bed_index)

    output:
    tuple val(meta), path('*.global.dist.txt')      , emit: global_txt
    tuple val(meta), path('*.summary.txt')          , emit: summary_txt
    tuple val(meta), path('*.region.dist.txt')      , optional:true, emit: regions_txt
    tuple val(meta), path('*.per-base.d4')          , optional:true, emit: per_base_d4
    tuple val(meta), path('*.per-base.bed.gz')      , optional:true, emit: per_base_bed
    tuple val(meta), path('*.per-base.bed.gz.csi')  , optional:true, emit: per_base_csi
    tuple val(meta), path('*.regions.bed.gz')       , optional:true, emit: regions_bed
    tuple val(meta), path('*.regions.bed.gz.csi')   , optional:true, emit: regions_csi
    tuple val(meta), path('*.quantized.bed.gz')     , optional:true, emit: quantized_bed
    tuple val(meta), path('*.quantized.bed.gz.csi') , optional:true, emit: quantized_csi
    tuple val(meta), path('*.thresholds.bed.gz')    , optional:true, emit: thresholds_bed
    tuple val(meta), path('*.thresholds.bed.gz.csi'), optional:true, emit: thresholds_csi
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--fasta ${fasta}" : ""
    def interval = bed ? "--by ${bed}" : ""
    if (bed && args.contains("--by")) {
        error "'--by' can only be specified once when running mosdepth! Either remove input BED file definition or remove '--by' from 'ext.args' definition"
    }
    if (!bed && args.contains("--thresholds")) {
        error "'--thresholds' can only be specified in conjunction with '--by'"
    }

    """
    mosdepth \\
        --threads $task.cpus \\
        $interval \\
        $reference \\
        $args \\
        $prefix \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.global.dist.txt
    touch ${prefix}.region.dist.txt
    touch ${prefix}.summary.txt
    touch ${prefix}.per-base.d4
    echo "" | gzip > ${prefix}.per-base.bed.gz
    touch ${prefix}.per-base.bed.gz.csi
    echo "" | gzip > ${prefix}.regions.bed.gz
    touch ${prefix}.regions.bed.gz.csi
    echo "" | gzip > ${prefix}.quantized.bed.gz
    touch ${prefix}.quantized.bed.gz.csi
    echo "" | gzip > ${prefix}.thresholds.bed.gz
    touch ${prefix}.thresholds.bed.gz.csi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    END_VERSIONS
    """
}

process MOSDEPTH_RAW {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mosdepth:0.3.8--hd299d5a_0' :
        'quay.io/biocontainers/mosdepth:0.3.8--hd299d5a_0'}"

    input:
    tuple val(meta),  path(bam), path(bai)
    tuple val(meta2), path(fasta)
    path(fai)
    path(bed)
    path(bed_index)

    output:
    tuple val(meta), path('*.global.dist.txt')      , emit: global_txt
    tuple val(meta), path('*.summary.txt')          , emit: summary_txt
    tuple val(meta), path('*.region.dist.txt')      , optional:true, emit: regions_txt
    tuple val(meta), path('*.per-base.d4')          , optional:true, emit: per_base_d4
    tuple val(meta), path('*.per-base.bed.gz')      , optional:true, emit: per_base_bed
    tuple val(meta), path('*.per-base.bed.gz.csi')  , optional:true, emit: per_base_csi
    tuple val(meta), path('*.regions.bed.gz')       , optional:true, emit: regions_bed
    tuple val(meta), path('*.regions.bed.gz.csi')   , optional:true, emit: regions_csi
    tuple val(meta), path('*.quantized.bed.gz')     , optional:true, emit: quantized_bed
    tuple val(meta), path('*.quantized.bed.gz.csi') , optional:true, emit: quantized_csi
    tuple val(meta), path('*.thresholds.bed.gz')    , optional:true, emit: thresholds_bed
    tuple val(meta), path('*.thresholds.bed.gz.csi'), optional:true, emit: thresholds_csi
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--fasta ${fasta}" : ""
    def interval = bed ? "--by ${bed}" : ""
    if (bed && args.contains("--by")) {
        error "'--by' can only be specified once when running mosdepth! Either remove input BED file definition or remove '--by' from 'ext.args' definition"
    }
    if (!bed && args.contains("--thresholds")) {
        error "'--thresholds' can only be specified in conjunction with '--by'"
    }

    """
    mosdepth \\
        --threads ${task.cpus} \\
        ${interval} \\
        ${reference} \\
        ${args} \\
        ${prefix}.raw \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.raw.global.dist.txt
    touch ${prefix}.raw.region.dist.txt
    touch ${prefix}.raw.summary.txt
    touch ${prefix}.raw.per-base.d4
    echo "" | gzip > ${prefix}.raw.per-base.bed.gz
    touch ${prefix}.raw.per-base.bed.gz.csi
    echo "" | gzip > ${prefix}.raw.regions.bed.gz
    touch ${prefix}.raw.regions.bed.gz.csi
    echo "" | gzip > ${prefix}.raw.quantized.bed.gz
    touch ${prefix}.raw.quantized.bed.gz.csi
    echo "" | gzip > ${prefix}.raw.thresholds.bed.gz
    touch ${prefix}.raw.thresholds.bed.gz.csi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    END_VERSIONS
    """
}

process MOSDEPTH_CON {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mosdepth:0.3.8--hd299d5a_0' :
        'quay.io/biocontainers/mosdepth:0.3.8--hd299d5a_0'}"

    input:
    tuple val(meta),  path(bam), path(bai)
    tuple val(meta2), path(fasta)
    path(fai)
    path(bed)
    path(bed_index)

    output:
    tuple val(meta), path('*.global.dist.txt')      , emit: global_txt
    tuple val(meta), path('*.summary.txt')          , emit: summary_txt
    tuple val(meta), path('*.region.dist.txt')      , optional:true, emit: regions_txt
    tuple val(meta), path('*.per-base.d4')          , optional:true, emit: per_base_d4
    tuple val(meta), path('*.per-base.bed.gz')      , optional:true, emit: per_base_bed
    tuple val(meta), path('*.per-base.bed.gz.csi')  , optional:true, emit: per_base_csi
    tuple val(meta), path('*.regions.bed.gz')       , optional:true, emit: regions_bed
    tuple val(meta), path('*.regions.bed.gz.csi')   , optional:true, emit: regions_csi
    tuple val(meta), path('*.quantized.bed.gz')     , optional:true, emit: quantized_bed
    tuple val(meta), path('*.quantized.bed.gz.csi') , optional:true, emit: quantized_csi
    tuple val(meta), path('*.thresholds.bed.gz')    , optional:true, emit: thresholds_bed
    tuple val(meta), path('*.thresholds.bed.gz.csi'), optional:true, emit: thresholds_csi
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--fasta ${fasta}" : ""
    def interval = bed ? "--by ${bed}" : ""
    if (bed && args.contains("--by")) {
        error "'--by' can only be specified once when running mosdepth! Either remove input BED file definition or remove '--by' from 'ext.args' definition"
    }
    if (!bed && args.contains("--thresholds")) {
        error "'--thresholds' can only be specified in conjunction with '--by'"
    }

    """
    mosdepth \\
        --threads ${task.cpus} \\
        ${interval} \\
        ${reference} \\
        ${args} \\
        ${prefix}.con \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.con.global.dist.txt
    touch ${prefix}.con.region.dist.txt
    touch ${prefix}.con.summary.txt
    touch ${prefix}.con.per-base.d4
    echo "" | gzip > ${prefix}.con.per-base.bed.gz
    touch ${prefix}.con.per-base.bed.gz.csi
    echo "" | gzip > ${prefix}.con.regions.bed.gz
    touch ${prefix}.con.regions.bed.gz.csi
    echo "" | gzip > ${prefix}.con.quantized.bed.gz
    touch ${prefix}.con.quantized.bed.gz.csi
    echo "" | gzip > ${prefix}.con.thresholds.bed.gz
    touch ${prefix}.con.thresholds.bed.gz.csi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    END_VERSIONS
    """
}

process MOSDEPTH_DUP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mosdepth:0.3.8--hd299d5a_0' :
        'quay.io/biocontainers/mosdepth:0.3.8--hd299d5a_0'}"

    input:
    tuple val(meta),  path(bam), path(bai)
    tuple val(meta2), path(fasta)
    path(fai)
    path(bed)
    path(bed_index)

    output:
    tuple val(meta), path('*.global.dist.txt')      , emit: global_txt
    tuple val(meta), path('*.summary.txt')          , emit: summary_txt
    tuple val(meta), path('*.region.dist.txt')      , optional:true, emit: regions_txt
    tuple val(meta), path('*.per-base.d4')          , optional:true, emit: per_base_d4
    tuple val(meta), path('*.per-base.bed.gz')      , optional:true, emit: per_base_bed
    tuple val(meta), path('*.per-base.bed.gz.csi')  , optional:true, emit: per_base_csi
    tuple val(meta), path('*.regions.bed.gz')       , optional:true, emit: regions_bed
    tuple val(meta), path('*.regions.bed.gz.csi')   , optional:true, emit: regions_csi
    tuple val(meta), path('*.quantized.bed.gz')     , optional:true, emit: quantized_bed
    tuple val(meta), path('*.quantized.bed.gz.csi') , optional:true, emit: quantized_csi
    tuple val(meta), path('*.thresholds.bed.gz')    , optional:true, emit: thresholds_bed
    tuple val(meta), path('*.thresholds.bed.gz.csi'), optional:true, emit: thresholds_csi
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--fasta ${fasta}" : ""
    def interval = bed ? "--by ${bed}" : ""
    if (bed && args.contains("--by")) {
        error "'--by' can only be specified once when running mosdepth! Either remove input BED file definition or remove '--by' from 'ext.args' definition"
    }
    if (!bed && args.contains("--thresholds")) {
        error "'--thresholds' can only be specified in conjunction with '--by'"
    }

    """
    mosdepth \\
        --threads ${task.cpus} \\
        ${interval} \\
        ${reference} \\
        ${args} \\
        ${prefix}.duplex \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.duplex.global.dist.txt
    touch ${prefix}.duplex.region.dist.txt
    touch ${prefix}.duplex.summary.txt
    touch ${prefix}.duplex.per-base.d4
    echo "" | gzip > ${prefix}.duplex.per-base.bed.gz
    touch ${prefix}.duplex.per-base.bed.gz.csi
    echo "" | gzip > ${prefix}.duplex.regions.bed.gz
    touch ${prefix}.duplex.regions.bed.gz.csi
    echo "" | gzip > ${prefix}.duplex.quantized.bed.gz
    touch ${prefix}.duplex.quantized.bed.gz.csi
    echo "" | gzip > ${prefix}.duplex.thresholds.bed.gz
    touch ${prefix}.duplex.thresholds.bed.gz.csi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    END_VERSIONS
    """
}

process MOSDEPTH_SIM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mosdepth:0.3.8--hd299d5a_0' :
        'quay.io/biocontainers/mosdepth:0.3.8--hd299d5a_0'}"

    input:
    tuple val(meta),  path(bam), path(bai)
    tuple val(meta2), path(fasta)
    path(fai)
    path(bed)
    path(bed_index)

    output:
    tuple val(meta), path('*.global.dist.txt')      , emit: global_txt
    tuple val(meta), path('*.summary.txt')          , emit: summary_txt
    tuple val(meta), path('*.region.dist.txt')      , optional:true, emit: regions_txt
    tuple val(meta), path('*.per-base.d4')          , optional:true, emit: per_base_d4
    tuple val(meta), path('*.per-base.bed.gz')      , optional:true, emit: per_base_bed
    tuple val(meta), path('*.per-base.bed.gz.csi')  , optional:true, emit: per_base_csi
    tuple val(meta), path('*.regions.bed.gz')       , optional:true, emit: regions_bed
    tuple val(meta), path('*.regions.bed.gz.csi')   , optional:true, emit: regions_csi
    tuple val(meta), path('*.quantized.bed.gz')     , optional:true, emit: quantized_bed
    tuple val(meta), path('*.quantized.bed.gz.csi') , optional:true, emit: quantized_csi
    tuple val(meta), path('*.thresholds.bed.gz')    , optional:true, emit: thresholds_bed
    tuple val(meta), path('*.thresholds.bed.gz.csi'), optional:true, emit: thresholds_csi
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--fasta ${fasta}" : ""
    def interval = bed ? "--by ${bed}" : ""
    if (bed && args.contains("--by")) {
        error "'--by' can only be specified once when running mosdepth! Either remove input BED file definition or remove '--by' from 'ext.args' definition"
    }
    if (!bed && args.contains("--thresholds")) {
        error "'--thresholds' can only be specified in conjunction with '--by'"
    }

    """
    mosdepth \\
        --threads ${task.cpus} \\
        ${interval} \\
        ${reference} \\
        ${args} \\
        ${prefix}.simplex \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.simplex.global.dist.txt
    touch ${prefix}.simplex.region.dist.txt
    touch ${prefix}.simplex.summary.txt
    touch ${prefix}.simplex.per-base.d4
    echo "" | gzip > ${prefix}.simplex.per-base.bed.gz
    touch ${prefix}.simplex.per-base.bed.gz.csi
    echo "" | gzip > ${prefix}.simplex.regions.bed.gz
    touch ${prefix}.simplex.regions.bed.gz.csi
    echo "" | gzip > ${prefix}.simplex.quantized.bed.gz
    touch ${prefix}.simplex.quantized.bed.gz.csi
    echo "" | gzip > ${prefix}.simplex.thresholds.bed.gz
    touch ${prefix}.simplex.thresholds.bed.gz.csi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    END_VERSIONS
    """
}

process MOSDEPTH_DR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mosdepth:0.3.8--hd299d5a_0' :
        'quay.io/biocontainers/mosdepth:0.3.8--hd299d5a_0'}"

    input:
    tuple val(meta),  path(bam), path(bai)
    tuple val(meta2), path(fasta)
    path(fai)
    path(bed)
    path(bed_index)

    output:
    tuple val(meta), path('*.global.dist.txt')      , emit: global_txt
    tuple val(meta), path('*.summary.txt')          , emit: summary_txt
    tuple val(meta), path('*.region.dist.txt')      , optional:true, emit: regions_txt
    tuple val(meta), path('*.per-base.d4')          , optional:true, emit: per_base_d4
    tuple val(meta), path('*.per-base.bed.gz')      , optional:true, emit: per_base_bed
    tuple val(meta), path('*.per-base.bed.gz.csi')  , optional:true, emit: per_base_csi
    tuple val(meta), path('*.regions.bed.gz')       , optional:true, emit: regions_bed
    tuple val(meta), path('*.regions.bed.gz.csi')   , optional:true, emit: regions_csi
    tuple val(meta), path('*.quantized.bed.gz')     , optional:true, emit: quantized_bed
    tuple val(meta), path('*.quantized.bed.gz.csi') , optional:true, emit: quantized_csi
    tuple val(meta), path('*.thresholds.bed.gz')    , optional:true, emit: thresholds_bed
    tuple val(meta), path('*.thresholds.bed.gz.csi'), optional:true, emit: thresholds_csi
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--fasta ${fasta}" : ""
    def interval = bed ? "--by ${bed}" : ""
    if (bed && args.contains("--by")) {
        error "'--by' can only be specified once when running mosdepth! Either remove input BED file definition or remove '--by' from 'ext.args' definition"
    }
    if (!bed && args.contains("--thresholds")) {
        error "'--thresholds' can only be specified in conjunction with '--by'"
    }

    """
    mosdepth \\
        --threads $task.cpus \\
        $interval \\
        $reference \\
        $args \\
        $prefix \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.global.dist.txt
    touch ${prefix}.region.dist.txt
    touch ${prefix}.summary.txt
    touch ${prefix}.per-base.d4
    echo "" | gzip > ${prefix}.per-base.bed.gz
    touch ${prefix}.per-base.bed.gz.csi
    echo "" | gzip > ${prefix}.regions.bed.gz
    touch ${prefix}.regions.bed.gz.csi
    echo "" | gzip > ${prefix}.quantized.bed.gz
    touch ${prefix}.quantized.bed.gz.csi
    echo "" | gzip > ${prefix}.thresholds.bed.gz
    touch ${prefix}.thresholds.bed.gz.csi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    END_VERSIONS
    """
}
