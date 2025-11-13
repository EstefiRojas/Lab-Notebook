# WARP.md

This file provides guidance to WARP (warp.dev) when working with code in this repository.

## Project Overview

This is a functional genomics research project focused on analyzing various genomic and epigenetic features to distinguish between functional and non-functional RNA elements (mRNAs, lncRNAs, and sncRNAs). The analysis combines intrinsic sequence features, epigenetic marks, and gene functionality metrics to characterize RNA functionality.

## Quick Start Commands

### Docker-based Analysis (Recommended)
```bash
# Build the analysis environment
docker build -t functional-genomics-analysis working_scripts/

# Run the complete analysis pipeline
docker run --rm -v "$(pwd)/results:/results" functional-genomics-analysis
```

### Local R Analysis
```bash
cd working_scripts/

# Restore R environment
R -e "renv::restore()"

# Run complete analysis via R Markdown
R -e "rmarkdown::render('workflow.Rmd')"

# Alternative: run via master script
R -e "source('run_all.R')"
```

### Individual Analysis Components
```bash
# Process epigenetic data
bash scripts/download_epigenetic_data.sh
bash scripts/histone_processing.sh
bash scripts/chrm_acc_processing.sh
bash scripts/methylome_processing.sh

# Compute dinucleotide frequencies
bash scripts/dinucleotide_frequencies.sh

# Join datasets
bash scripts/join_datasets.sh
bash scripts/join_epigenetic_features.sh

# GWAS data processing
bash scripts/get_GWAS_associations.sh
bash scripts/fetch_ensembl_data.sh
```

## Architecture and Project Structure

### Data Processing Pipeline Architecture

The analysis follows a multi-stage pipeline:

1. **Data Loading Phase**: Raw genomic features are loaded from CSV files for functional and negative control datasets
2. **Feature Engineering Phase**: Dinucleotide frequencies, epigenetic marks, and functionality scores are computed
3. **Statistical Analysis Phase**: Robust Z-scores are computed, distance effects analyzed, and dimensionality reduction performed
4. **Visualization Phase**: Heatmaps, violin plots, jitter plots, and PCA visualizations are generated

### Core Component Structure

- **Configuration System**: `working_scripts/scripts/config.R` defines all file paths and analysis parameters
- **Data Loaders**: Modular scripts for loading different feature types (gene functionality, epigenetic, dinucleotide)
- **Statistical Utilities**: Custom implementations of K-S tests, robust z-score calculations, and correlation analysis
- **Visualization Engine**: Standardized plotting functions for violin plots, heatmaps, and scatter plots

### Data Flow Architecture

```
Raw Data → Feature Extraction → Z-score Normalization → Statistical Analysis → Visualization
     ↓              ↓                    ↓                      ↓              ↓
   CSV files    Shell scripts       R functions           R analysis      ggplot2 plots
```

### Gene Type Classification System

The analysis works with three main RNA types, each with functional (+) and negative control (-) variants:
- **mRNA**: Protein-coding transcripts (exon2, exon3)
- **lncRNA**: Long non-coding RNAs (exon1, exon2) 
- **sncRNA**: Short non-coding RNAs

### Feature Categories

1. **Intrinsic Features**: Sequence-based properties (dinucleotide content, structure, conservation)
2. **Epigenetic Features**: Chromatin modifications (histones, chromatin accessibility, methylation)
3. **Functionality Features**: Expression levels, evolutionary conservation, population variation

## Key Analysis Methods

- **Robust Z-score Normalization**: Custom implementation separating functional from non-functional elements
- **Kolmogorov-Smirnov Testing**: Statistical comparison between functional and control groups
- **Distance Effect Analysis**: Spearman correlation analysis of genomic distance effects
- **Dimensionality Reduction**: PCA and MDS for visualizing separation between gene types
- **Multi-scale Visualization**: Violin plots, heatmaps, and scatter plots with standardized aesthetics

## Development Environment

- **R Version**: 4.4.1 (specified in renv.lock and Dockerfile)
- **Package Management**: renv for reproducible R environments
- **Containerization**: Docker with rocker/r-ver base image
- **Key Dependencies**: ggplot2, dplyr, Hmisc, corrplot, viridis, rmarkdown

## External Data Integration

The project integrates with multiple external genomic databases and validation datasets, including the Liang et al. 2024 gRNA essentiality screen for lncRNA functional validation using BLAT alignment and model prediction comparison.

## Testing and Validation

Run the complete analysis to verify data processing and statistical computations:
```bash
# Validate analysis pipeline
cd working_scripts/
R -e "source('run_all.R')"

# Check for expected output directories
ls -la results/
```
