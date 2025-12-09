# Unified lncRNA Essentiality Analysis Workflow

This document outlines the step-by-step process to run the complete unified analysis, from data unification to visualization.


## Step 1: Unify Datasets

Combine gRNA data from the five studies (Huang, Liang, Liu, Montero, Zhu) into a single standardized format.

```bash
cd scripts
./unify_datasets.sh
```

**Output:** `data/unified/unified_gRNAs.csv`

## Step 2: Run Unified Analysis

Perform genome alignment (BLAT), protein off-target checking, and integrate model predictions. This step also applies the Liang study filters (unique hits & best hit selection).

```bash
./run_unified_analysis.sh
```

**Output:** `results/unified_genome_alignments.csv`

## Step 3: Annotate Gene Names

Map the Ensembl IDs (ENSG) to gene names using the GENCODE GTF file.

```bash
./map_lncRNA_genes.sh ../results/unified_genome_alignments.csv ../data/references/gencode.v49.long_noncoding_RNAs.gtf
```

**Output:** `results/annotated_unified_genome_alignments.csv`

## Step 4: Generate Essentiality Matrix

Create a matrix summarizing the essentiality status (Rare, Common, Core) of each lncRNA across the four studies, including functional probability scores.

```bash
./generate_essentiality_matrix.sh
```

**Output:** `results/essentiality_matrix.csv`

## Step 5: Visualize Results

Generate various plots to analyze the consensus and functional probability of essential lncRNAs.

```bash
python3 plot_essentiality_analysis.py
```

**Output:** `results/plots/`
*   `essentiality_heatmap.png`
*   `probability_vs_essentiality_boxplot.png`
*   `probability_vs_essentiality_violin_pooled.png`
*   `study_correlation.png`
*   `upset_plot.png`
*   `consensus_vs_probability.png`
*   `parallel_categories.png`
*   `top_20_ranking.png`
*   `jaccard_similarity.png`
*   `upset_binary_matrix.csv`

## Step 6: Generate UpSet Plot (R)

Generate the refined UpSet plot using the R script, which uses the binary matrix created in the previous step.

```bash
Rscript plot_upset.R
```

**Output:** `results/plots2/upset_plot_R.png`
