```markdown name=README.md url=https://github.com/Sasivarman06/Variant-calling-pipeline/blob/main/README.md
# Variant-calling-pipeline

A Snakemake-based Whole Genome Sequencing (WGS) and Whole Exome Sequencing (WES) variant calling pipeline using BWA, GATK, and bcftools.

## Overview

This pipeline automates the complete variant calling workflow for genomic analysis, from raw sequencing reads to variant calls. It leverages Snakemake for workflow orchestration and combines powerful tools for alignment, variant detection, and quality control. 

## Features

- **Automated Workflow**: Snakemake-based pipeline for reproducible analysis
- **Read Alignment**: BWA-based genome mapping
- **Variant Discovery**: GATK for accurate variant calling
- **Variant Analysis**: bcftools for VCF manipulation and filtering
- **WGS/WES Support**: Suitable for both whole genome and exome sequencing projects
- **Configurable**:  Easy-to-modify configuration files for custom analysis parameters

## Pipeline Components

### Tools & Technologies

- **Snakemake**: Workflow management system
- **BWA**: Burrows-Wheeler Aligner for read alignment
- **GATK**:  Genome Analysis Toolkit for variant calling
- **bcftools**: Tools for VCF/BCF file processing
- **Python**: Pipeline scripting and automation

### Directory Structure

```
Variant-calling-pipeline/
├── config/              # Configuration files for pipeline parameters
├── reference/           # Reference genome files
├── tb_data/             # Sample/test data directory
├── workflow/            # Snakemake workflow rules and scripts
└── . gitattributes       # Git attributes configuration
```

## Getting Started

### Prerequisites

- Snakemake
- BWA
- GATK
- bcftools
- Python 3.x

### Installation

1. Clone the repository:
```bash
git clone https://github.com/Sasivarman06/Variant-calling-pipeline.git
cd Variant-calling-pipeline
```

2. Install required dependencies (if using conda):
```bash
conda install -c bioconda snakemake bwa gatk4 bcftools
```

### Usage

1. **Configure Parameters**:  Edit the configuration files in the `config/` directory to specify:
   - Input data paths
   - Reference genome location
   - Pipeline parameters

2. **Prepare Data**: Place your sequencing data in the appropriate directory and reference files in `reference/`

3. **Run the Pipeline**:
```bash
snakemake --use-conda --cores <number_of_cores>
```

## Workflow Steps

The pipeline typically follows these steps:

1. **Alignment**: Map reads to reference genome using BWA
2. **Post-alignment Processing**: Sort and index aligned BAM files
3. **Variant Calling**: Detect variants using GATK HaplotypeCaller or similar
4. **Variant Filtering**: Filter and annotate variants using bcftools
5. **Quality Control**: Generate statistics and quality reports

## Configuration

Key configuration files are located in the `config/` directory. Modify these files to customize:
- Reference genome paths
- Input/output directories
- Variant calling parameters
- Resource allocation (CPU, memory)

## Output

The pipeline generates: 
- Aligned BAM files
- VCF (Variant Call Format) files containing detected variants
- Quality control reports
- Log files for troubleshooting
```
