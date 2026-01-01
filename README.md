# ATAC-seq Snakemake Pipeline

A comprehensive, reproducible Snakemake workflow for processing and analyzing ATAC-seq (Assay for Transposase-Accessible Chromatin with sequencing) data. This pipeline performs end-to-end analysis from raw FASTQ files to peak calling, annotation, and reproducibility assessment.

## Features

- **Quality Control**: Fastp adapter trimming, alignment quality metrics, TSS enrichment analysis
- **Flexible Alignment**: Supports both Bowtie2 and BWA-MEM2 aligners
- **Prealignment Filtering**: Optional prealignment to filter reads (e.g., mitochondrial DNA)
- **Peak Calling**: MACS2 peak calling with Tn5 shift correction
- **Reproducibility Analysis**: 
  - IDR (Irreproducible Discovery Rate) analysis for replicate consistency
  - **Replicate Correlation**: Pearson and Spearman correlation metrics based on read counts in merged peaks
  - Self-consistency (pseudo-replicates) and true replicate comparisons
  - Optimal peak selection following ENCODE standards
- **Annotation**: UROPA-based peak annotation to genomic features
- **Quantification**: Read counts in consensus peaks, promoter regions, and TSS regions
- **Motif Analysis**: HOMER motif enrichment analysis
- **Visualization**: MultiQC reports, correlation scatter plots, BigWig tracks

## Requirements

- **Snakemake** >= 8.20.1
- **Conda** or **Mamba** (for environment management)
- **Python** 3.8+
- **PEP** (Portable Encapsulated Projects) for sample metadata management

## Installation

1. Clone the repository:
```bash
git clone <repository-url>
cd ATAC-sm
```

2. Install Snakemake (if not already installed):
```bash
conda install -c bioconda snakemake
# or
mamba install -c bioconda snakemake
```

3. Install PEP (if not already installed):
```bash
pip install peppy
```

4. The pipeline uses conda environments defined in `workflow/envs/`. These will be automatically created when you run the pipeline.

## Quick Start

1. **Prepare your configuration files**:
   - Edit `config/config.yaml` with your project settings
   - Edit `pep/project_config.yaml` with your path variables
   - Prepare your sample annotation CSV file (see Configuration section)

2. **Run the pipeline**:
```bash
cd workflow
snakemake --use-conda --cores <number_of_cores>
```

3. **For cluster execution** (SLURM example):
```bash
snakemake --use-conda --profile profiles/default --cores <number_of_cores>
```

## Configuration

### Project Configuration (`config/config.yaml`)

Key configuration sections:

#### Project Settings
```yaml
project:
  name: your_project_name
  genome: hg38  # or mm10, etc.
  out_dir: "{process_dir}/your_output_directory"
  read_type: paired  # or "single"
```

#### Alignment Configuration
```yaml
alignment:
  tool: bwa-mem2  # or "bowtie2"
  sequencing_platform: illumina
  
  # Optional prealignment (e.g., filter mitochondrial reads)
  prealign:
    enabled: true
    indices:
      - name: chrM
        path: "{database_dir}/hg38_rcrds/rcrds"
  
  # BWA-MEM2 specific
  bwa:
    index: "{shared_genome_dir}/bwa-mem2/GRCh38"
    min_score: null
    extra_args: ""
  
  # Bowtie2 specific
  bowtie2:
    index: "{database_dir}/hg38/indices_for_Bowtie2/hg38"
```

#### Peak Calling Configuration
```yaml
peaks:
  macs2_keep_dup: "all"
  macs2_shift: -75      # Tn5 shift correction
  macs2_extsize: 150
  macs2_pval: 1e-4
  disable_tn5_shift: false
```

#### Reproducibility Analysis
```yaml
replicates:
  pval_thresh: 1e-3
  smooth_win: 150
  cap_num_peak: 300000
  peak_combination_method: "idr"  # or "overlap"
  idr_thresh: 0.05
  idr_rank: "signal.value"
```

### PEP Configuration (`pep/project_config.yaml`)

Define path variables that can be referenced in the main config:

```yaml
paths:
  data_dir: "/path/to/raw/data"
  process_dir: "/path/to/output"
  database_dir: "/path/to/reference/genomes"
  shared_genome_dir: "/path/to/shared/genomes"
```

### Sample Annotation File

The sample annotation CSV should contain at minimum:
- `sample_name`: Unique identifier for each sample
- `run`: Run number (integer)
- `R1`: Path to forward read FASTQ file
- `R2`: Path to reverse read FASTQ file (for paired-end)
- `read_type`: "paired" or "single"

Optional columns:
- `replicate_sample_name`: Group identifier for biological replicates (enables reproducibility analysis)
- `passqc`: QC flag (1 = pass, 0 = fail)

Example:
```csv
sample_name,run,R1,R2,read_type,replicate_sample_name
sample1,1,/path/to/sample1_R1.fq.gz,/path/to/sample1_R2.fq.gz,paired,group1
sample1,2,/path/to/sample1_run2_R1.fq.gz,/path/to/sample1_run2_R2.fq.gz,paired,group1
sample2,1,/path/to/sample2_R1.fq.gz,/path/to/sample2_R2.fq.gz,paired,group1
```

## Pipeline Workflow

The pipeline consists of several main phases:

### 1. Preprocessing
- **Adapter Trimming**: Fastp removes adapters and performs quality filtering
- **Prealignment** (optional): Filters reads against specified indices (e.g., mitochondrial DNA)

### 2. Alignment
- Alignment using Bowtie2 or BWA-MEM2
- Duplicate marking with Samblaster
- Quality filtering and sorting

### 3. Peak Calling
- MACS2 peak calling per sample
- Tn5 shift correction for ATAC-seq
- Peak merging across samples

### 4. Reproducibility Analysis (for replicate groups)
- **Pseudo-replicate Generation**: Splits each sample into two pseudo-replicates
- **Self-Consistency**: IDR between pseudo-replicates (N1, N2, ...)
- **True Replicate Comparison**: IDR between biological replicates (Nt)
- **Pooled Pseudo-Replicates**: IDR on pooled data (Np)
- **Replicate Correlation**: Pearson and Spearman correlation based on read counts in group-merged peaks
- **Optimal Peak Selection**: ENCODE-standard selection (max(Nt, Np))

### 5. Annotation & Quantification
- Peak annotation to genomic features (promoters, TSS, etc.)
- Read counting in consensus peaks
- Promoter and TSS region quantification

### 6. Downstream Analysis
- HOMER motif enrichment
- BigWig track generation
- MultiQC report generation

## Contributing

Contributions are welcome! Please open an issue or submit a pull request.
