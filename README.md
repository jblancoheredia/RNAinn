[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/CMOinn/rnainn)

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/RNAinn_logo_dark.svg">
  <img alt="mskcc/cti/rnainn" src="assets/RNAinn_logo_light.svg" width="400">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/RNAinn_metro_dark.svg">
  <img alt="Metro" src="assets/RNAinn_metro_light.svg" width="1500">
</picture>

## Introduction

**MSKCC/CTI/RNAinn** is an nf-core borne, production-ready and comprehensive bioinformatics pipeline for RNA sequencing data processing.

The pipeline processes paired-end RNA-seq data through multiple analysis modules including UMI processing, quality control, alignment, gene expression quantification, fusion detection and variant calling. RNAinn is designed to handle both standard RNA-seq workflows and specialized analyses for cancer genomics research.

The pipeline includes six main analysis tracks:

**Preprocessing**: Quality control, read trimming, and initial metrics collection
- Read QC with ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) and ([`MultiQC`](http://multiqc.info/)) reporting
- Adapter trimming with FastP
- Read counting and sampling with Seqtk

**UMI Processing**: Molecular barcoding analysis with fgbio toolkit
- UMI correction and grouping
- Consensus read calling
- Error rate analysis and quality filtering

**Deduplication**: Molecular barcoding analysis with fgbio toolkit
- UMI correction and grouping
- Consensus read calling
- Error rate analysis and quality filtering

**Gene Expression**: Multiple quantification methods for comprehensive expression analysis
- STAR alignment with two-pass mode
- Transcript quantification with Kallisto, Salmon, StringTie, and RSEM
- Gene-level counting with FeatureCounts
- DESeq2-based quality control metrics

**Fusion Calling**: Multi-tool fusion detection and reporting
- Arriba fusion detection
- STAR-Fusion analysis
- FusionCatcher identification
- Integration of the calls with FusionReport
- FusionInspector validation
- Visualization with FusViz

**Variant Calling**: GATK-based somatic variant detection
- Duplicate marking and base quality score recalibration
- HaplotypeCaller for variant discovery
- Variant filtering and annotation with SNPeff and VEP
- Support for targeted sequencing intervals

The first track is mandatory but highly configurable e.g. the parameter --run_downsamplings set as True, will cap all your samples in the run to the lowest of them, usefull to comparative experiments, or you can set the number of downsampled reads with --ds_totalreads_aim. The rest are optional, you can enable those by setting to True the parameters --run_umiprocessing, --run_fusion_splice, --run_genexpression true, and run_variantdscvry the defaults are set to False so you can compose the run you need.

## Quick Start

Create a samplesheet with your input data:

`samplesheet.csv`:
```csv
sample,fastq_1,fastq_2,strandedness
SAMPLE_1,sample1_R1.fastq.gz,sample1_R2.fastq.gz,reverse
SAMPLE_2,sample2_R1.fastq.gz,sample2_R2.fastq.gz,reverse
```

Run the pipeline:

```bash
nextflow run CMOinn/rnainn \
   -profile <crater/juno/iris/singularity/docker> \
   --input samplesheet.csv \
   --outdir results \
   --genome GRCh38
```

## Usage

### Input Requirements

The pipeline requires a CSV samplesheet with the following columns:
- `sample`: Unique sample identifier
- `fastq_1`: Path to first read file (R1)
- `fastq_2`: Path to second read file (R2) 
- `strandedness`: Library strandedness (`forward`, `reverse`, or `unstranded`)

### Core Parameters

```bash
--input                 # Input samplesheet (required)
--outdir                # Output directory (required)
--genome                # Reference genome [default: GRCh38]
--seq_library           # Sequencing library type [default: Av1]
--seq_center            # Sequencing center [default: MSKCC_IGO]
```

### Analysis Modules

Enable specific analysis modules:

```bash
--run_genexpression     # Gene expression analysis [default: false]
--run_fusion_splice     # Fusion detection [default: false]
--run_variantdscvry     # Variant calling [default: false]
--run_umiprocessing     # UMI processing [default: false]
--run_copynumberalt     # Copy number analysis [default: false]
```

### Tool Selection

Control which tools to run:

```bash
--arriba                # Enable Arriba fusion caller [default: true]
--starfusion            # Enable STAR-Fusion [default: true]
--fusioncatcher         # Enable FusionCatcher [default: true]
--stringtie             # Enable StringTie [default: true]
--portcullis            # Enable Portcullis [default: true]
```

### Advanced Options

```bash
--star_twopass          # Enable STAR two-pass mode [default: true]
--fastp_trim            # Enable FastP trimming [default: true]
--remove_duplicates     # Remove duplicate reads [default: false]
--tools_cutoff          # Minimum tools for fusion filtering [default: 1]
--read_length           # Expected read length [default: 100]
```

### UMI Processing Options

```bash
--f2b_read_structure    # UMI read structure [default: '3M2S+T 3M2S+T']
--group_strategy        # UMI grouping strategy [default: 'edit']
--call_min_reads        # Minimum reads for consensus [default: 1]
--filter_min_reads      # Minimum reads for filtering [default: 2]
```

### Profiles

The pipeline includes several pre-configured profiles:

- `crater`: LSF cluster configuration for MSKCC Crater
- `juno`: LSF cluster configuration for MSKCC Juno
- `iris`: SLURM cluster configuration for MSKCC Iris
- `singularity`: Generic Singularity configuration
- `docker`: Docker container execution
- `conda`: Conda environment management

## Output

The pipeline generates comprehensive outputs organized by analysis module:

```
results/
├── preprocessing/
│   ├── fastqc/
│   ├── multiqc/
│   └── sampling/
├── alignment/
│   ├── star/
│   ├── samtools/
│   └── markduplicates/
├── quantification/
│   ├── kallisto/
│   ├── salmon/
│   ├── rsem/
│   ├── stringtie/
│   └── featurecounts/
├── fusion_analysis/
│   ├── arriba/
│   ├── starfusion/
│   ├── fusioncatcher/
│   ├── fusioninspector/
│   └── fusionreport/
├── variant_calling/
│   ├── gatk_haplotypecaller/
│   ├── variant_filtering/
│   └── variant_annotation/
├── umi_processing/
│   ├── fgbio_consensus/
│   ├── error_correction/
│   └── quality_filtering/
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    └── execution_trace.txt
```

## System Requirements

### Computational Resources

Default resource limits (can be increased via parameters):
- Memory: 64 GB
- CPUs: 12 cores  
- Time: 24 hours

### Reference Data

The pipeline requires pre-built reference datasets including:
- STAR genome indices
- Kallisto/Salmon transcript indices
- GATK reference files (FASTA, VCF, intervals)
- Fusion detection databases (Arriba, STAR-Fusion, FusionCatcher)
- Annotation files (GTF, BED, RefFlat)

Reference paths are configured per-cluster in the profile configurations.

## Configuration

### Custom Configuration

Create a custom configuration file:

```groovy
params {
    // Override default parameters
    max_memory = '128.GB'
    max_cpus = 32
    max_time = '48.h'
    
    // Custom tool parameters
    star_max_intron_size = 1000000
    gatk_hc_call_conf = 30
}

process {
    // Process-specific overrides
    withName: 'STAR_ALIGN' {
        cpus = 16
        memory = '64.GB'
    }
}
```

Run with custom config:
```bash
nextflow run CMOinn/rnainn -c custom.config [other options]
```

### Cluster Profiles

The pipeline includes optimized configurations for MSKCC compute clusters. Each profile sets appropriate:
- Job scheduler settings (LSF/SLURM)
- Container configurations  
- Reference data paths
- Resource allocations

## Troubleshooting

### Common Issues

**Memory errors**: Increase `--max_memory` parameter or use cluster profile
**Timeout errors**: Increase `--max_time` parameter  
**Missing references**: Ensure correct profile is selected for your cluster
**Tool failures**: Check individual tool logs in work directories

### Getting Help

For issues specific to RNAinn:
1. Check the pipeline documentation
2. Review error logs in the work directory
3. Contact the development team

For Nextflow-related issues:
- [Nextflow documentation](https://www.nextflow.io/docs/latest/)
- [nf-core help](https://nf-co.re/help)

## Credits and Citations

RNAinn was developed by the Technology Innovation lab at the Marie-Josée and Henry R. Kravis Center for Molecular Oncology (CMO), Memorial Sloan Kettering Cancer Center (MSKCC).

**Lead Developer**: blancoj@mskcc.org

### Tools and References

This pipeline integrates numerous bioinformatics tools. Key citations include:

- **STAR**: Dobin et al. Bioinformatics 2013
- **Kallisto**: Bray et al. Nature Biotechnology 2016  
- **Salmon**: Patro et al. Nature Methods 2017
- **Arriba**: Uhrig et al. Genome Biology 2021
- **GATK**: McKenna et al. Genome Research 2010
- **fgbio**: Fulcrum Genomics toolkit

A complete list of references can be found in [`CITATIONS.md`](CITATIONS.md).

### Framework

This pipeline uses code and infrastructure developed by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.