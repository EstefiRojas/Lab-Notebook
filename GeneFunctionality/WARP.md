# WARP.md

This file provides guidance to WARP (warp.dev) when working with code in this repository.

## Project Overview

This project validates lncRNA functionality predictions by comparing model predictions against experimental essentiality data from multiple CRISPR screening studies (Liang et al. 2024, Huang et al. 2024, Montero et al. 2024, Liu et al. 2017, Zhu et al. 2016). The analysis uses BLAT sequence alignment to map gRNA sequences to lncRNA exons, then correlates essentiality status with model-predicted functionality probabilities.

## Key Commands

### Core Analysis Pipeline (Study-specific)

Each study follows a similar pattern. Examples shown for Liang et al.:

```bash
cd scripts/

# Step 1: Convert CSV data to FASTA format
./csv_to_fasta.sh ../data/Liang/LiangMuller_May2025_Table1_Filtered-gRNAs.csv 1-2 4 3 > ../data/Liang/processed/Supplementary_tableS1_gRNAs.fasta

# Step 2: Run BLAT alignments against model predictions
blat -t=dna -q=dna ../data/model_predictions/lncrna_exon1.fasta ../data/Liang/processed/Supplementary_tableS1_gRNAs.fasta -minScore=15 -minIdentity=100 ../data/Liang/processed/tableS1_gRNAs_vs_exon1.psl

# Step 3: Join BLAT matches with model predictions and essentiality data
./join_blat_matches_Liang.sh ../data/Liang/processed/tableS1_gRNAs_vs_exon1.psl ../data/model_predictions/gencode-lncrna-ranking.csv ../data/Liang/LiangMuller_May2025-Table3-Gene-RRA-Ranking.csv > ../data/Liang/processed/annotated_tableS3_gRNAs_vs_exon1_prob.csv

# Step 4: Generate visualizations (run in R)
Rscript essential_prob_distribution_Liang_newData.R
```

### Genome-wide gRNA Mapping Pipeline

For studies requiring genome-wide alignment (Liang, Liu):

```bash
cd scripts/

# Map gRNAs to full genome and intersect with lncRNA annotations
./find_lncRNA_guides.sh ../data/Liang/LiangMuller_May2025_Table1_Filtered-gRNAs.csv ../data/Liang/LiangMuller_May2025-Table3-Gene-RRA-Ranking.csv

# Filter to unique hits
./filter_unique_hits.sh ../results/gRNA_lncRNA_matches.tsv

# Add model probability predictions
./add_probability_field.sh ../results/gRNA_lncRNA_matches_unique_sorted.tsv ../data/model_predictions/gencode-lncrna-ranking.csv

# Select best hit per target gene
./select_best_hit.sh ../results/gRNA_lncRNA_matches_with_prob.tsv
```

### Reference Data Preparation

```bash
cd scripts/

# Extract exon sequences for genes without predictions
./get_exon_sequences.sh ../results/gRNA_lncRNA_matches_unique_ensg_na_prob.tsv

# Analyze GTF annotation statistics
./analyze_gtf_stats.sh ../data/references/gencode.v49.long_noncoding_RNAs.gtf
```

## Architecture

### Data Flow

```
Study CSV → csv_to_fasta → FASTA sequences
                              ↓
                         BLAT alignment → PSL output
                              ↓
                    join_blat_matches_* → Annotated CSV with probabilities
                              ↓
         essential_prob_distribution_*.R → Statistical plots (violin/box/jitter)
```

### Alternative Genome-wide Pipeline

```
Study CSVs → find_lncRNA_guides → genome-wide BLAT + lncRNA intersection
                                      ↓
                              filter_unique_hits
                                      ↓
                             add_probability_field
                                      ↓
                              select_best_hit
                                      ↓
                              R visualization
```

### Script Categories

**Data Conversion**
- `csv_to_fasta.sh`: Convert CSV columns to FASTA format (supports ranges, lists, suffixes)
- `fasta_filter_Montero.sh`: Filter FASTA records by annotation type

**Alignment Processing**
- `find_lncRNA_guides.sh`: Genome-wide BLAT alignment + bedtools intersection with lncRNA GTF
- `find_lncRNA_guides_Liu.sh`: Liu et al.-specific genome alignment
- `join_blat_matches_*.sh`: Study-specific PSL parsers that extract gene IDs, join with model predictions, and annotate essentiality status

**Data Integration**
- `add_probability_field.sh`: Join BLAT results with model predictions by Ensembl ID
- `filter_unique_hits.sh`: Deduplicate gRNA-gene-ENSG combinations
- `select_best_hit.sh`: Keep highest probability lncRNA match per target gene
- `filter_na_prob.sh`: Extract genes without model predictions

**Sequence Extraction**
- `get_exon_sequences.sh`: Extract exon1/exon2 sequences from GTF for missing predictions
- `extract_exons.sh`: General exon extraction utility

**Statistical Analysis (R)**
- `essential_prob_distribution_*.R`: Generate violin/box plots with KS test statistics comparing probability distributions across essentiality groups
- `gwas_prob_distribution.R`: Analyze GWAS beta and p-value correlations with lncRNA probability

### Key Data Files

**Input Data Structure**
- `data/{Study}/`: Raw supplementary tables from publications
- `data/{Study}/processed/`: Generated FASTA, PSL, and annotated CSV files
- `data/model_predictions/`: lncRNA functionality predictions (exon1, exon2, gene-level)
- `data/references/`: GRCh38 genome, GENCODE lncRNA annotations

**Output Structure**
- `results/`: Analysis outputs (TSV matches, plots)
- `tmp/`: Temporary BLAT/BED files

### Study-Specific Parsing Logic

Each `join_blat_matches_*.sh` script implements study-specific identifier parsing:

**Liang et al.**: Extracts `Hum_XLOC_*` identifiers using regex for essentiality lookup
**Huang et al.**: Joins with coPARSE annotation and essential gene mapping
**Montero et al.**: Processes dual gRNA data (gRNA1 + gRNA2)
**Zhu et al.**: Handles sgRNA nomenclature (minimal matches found)

### Analysis Workflow Patterns

1. **Alignment-first approach**: Direct BLAT against model prediction sequences (faster, requires pre-computed predictions)
2. **Genome-wide approach**: BLAT against full genome → bedtools intersect with lncRNA GTF (finds novel genes, slower)

### R Visualization Components

All R scripts follow consistent patterns:
- Load exon1 and exon2 results, merge by gene ID
- Compute `max_prob` across exons
- Keep one record per gene (highest probability)
- Generate violin plots with:
  - Box plots overlaid (except for small groups)
  - KS test statistics between groups
  - Custom legend with sample sizes
  - ggplot2 + ggsignif for statistical annotations

### Prerequisites

**System Tools**
- BLAT (http://hgdownload.cse.ucsc.edu/admin/exe/)
- bedtools (https://bedtools.readthedocs.io/)
- wget, gunzip, awk, sed

**R Packages**
- dplyr, tidyr, stringr
- ggplot2, ggsignif
- viridis (for some plots)

## Common Workflows

### Adding a New Study

1. Create directory: `data/NewStudy/processed/`
2. Convert study data: `csv_to_fasta.sh`
3. Run BLAT: `blat -t=dna -q=dna ...`
4. Create study-specific join script: `join_blat_matches_NewStudy.sh`
5. Create R visualization: `essential_prob_distribution_NewStudy.R`

### Troubleshooting BLAT Matches

```bash
# Check number of alignments
wc -l ../data/Study/processed/file.psl

# Verify unique gene matches
awk 'NR > 5 {print $10}' ../data/Study/processed/file.psl | sort -u | wc -l

# Inspect specific gRNA alignment
grep "gRNA_ID" ../data/Study/processed/file.psl
```

### Analyzing GTF Features

```bash
# Count lncRNA genes
./analyze_gtf_stats.sh ../data/references/gencode.v49.long_noncoding_RNAs.gtf
```

## Data Processing Notes

- PSL files have 5-line header that must be skipped (use `NR > 5` or `FNR > 5`)
- CSV fields may have Windows line endings (`\r`) - scripts use `gsub(/\r$/,"")` to clean
- BLAT match criteria: `-minScore=15 -minIdentity=100` for stringent 100% identity matching
- Composite keys for joining: `GeneID SUBSEP TranscriptID` or similar patterns
- Model predictions use `|` delimiter: `ENSG123|ENST456|exon1` format
