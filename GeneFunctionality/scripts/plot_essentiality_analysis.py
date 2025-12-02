import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import os
import warnings

import sys
sys.setrecursionlimit(10000)

# Define paths
INPUT_FILE = "../results/annotated_unified_genome_alignments.csv"
OUTPUT_DIR = "../results/plots2"

# Create output directory
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

print(f"Loading and processing data from {INPUT_FILE}...")
raw_df = pd.read_csv(INPUT_FILE, low_memory=False)

# ---------------------------------------------------------
# Apply Filters (Matching R Script Logic)
# ---------------------------------------------------------
# R: filter(Protein_Off_Target == "NO" & !is.na(Probability_Functional) & Antisense_to_CDS == "NO")
print("Applying filters: Protein_Off_Target='NO', Antisense_to_CDS='NO', Probability_Functional not NA...")
raw_df['Probability_Functional'] = pd.to_numeric(raw_df['Probability_Functional'], errors='coerce')

filtered_df = raw_df[
    (raw_df['Protein_Off_Target'] == 'NO') &
    (raw_df['Antisense_to_CDS'] == 'NO') &
    (raw_df['Probability_Functional'].notna())
].copy()

print(f"Rows after filtering: {len(filtered_df)}")

# ---------------------------------------------------------
# Create Long-Format DataFrame (Matching R Logic for Counts/Distributions)
# ---------------------------------------------------------
# R script groups by (ENSG_ID, Study, Essentiality) and takes max Probability.
# This allows a gene to have multiple classifications in the same study (e.g., Rare AND Non-essential).
print("Creating Long-Format DataFrame (R-style grouping)...")
long_df = filtered_df.sort_values('Probability_Functional', ascending=False).drop_duplicates(['ENSG_ID', 'Study', 'Essentiality'])

# Print Counts to verify against R script
print("\n--- Counts per Study and Essentiality (Long DF) ---")
print(long_df.groupby(['Study', 'Essentiality']).size())
print("---------------------------------------------------\n")

# ---------------------------------------------------------
# Transform to Matrix Format (Gene x Study) for Heatmaps/Comparisons
# ---------------------------------------------------------
# For matrix plots, we need a single classification per Gene-Study pair.
# Strategy: Prioritize highest essentiality (Core > Common > Rare > Non-essential)
# and take the max Probability_Functional for the gene.

# 1. Define Essentiality Rank
ess_rank = {'Core': 3, 'Common': 2, 'Rare': 1, 'Non-essential': 0, '-': -1}
# Handle potential unexpected values safely
filtered_df['Ess_Rank'] = filtered_df['Essentiality'].map(ess_rank).fillna(-1)

# 2. Sort to prioritize higher rank
# We sort by ENSG_ID, Study, Ess_Rank (desc), Probability (desc)
filtered_df = filtered_df.sort_values(by=['ENSG_ID', 'Study', 'Ess_Rank', 'Probability_Functional'], ascending=[True, True, False, False])

# 3. Drop duplicates to get one entry per Gene-Study
# This picks the "most essential" classification for that study
study_entries = filtered_df.drop_duplicates(subset=['ENSG_ID', 'Study'])

# 4. Pivot to wide format
matrix_df = study_entries.pivot(index='ENSG_ID', columns='Study', values='Essentiality')
matrix_df = matrix_df.fillna('-') # Missing studies are '-'

# 5. Get Gene-level info (Probability and Name)
# Max Probability per gene
gene_probs = filtered_df.groupby('ENSG_ID')['Probability_Functional'].max()
# First Gene Name (if available)
gene_names = filtered_df.groupby('ENSG_ID')['Gene_Name'].first()

# 6. Merge into final df
df = matrix_df.join(gene_probs).join(gene_names)
df = df.reset_index()

# 7. Ensure all expected study columns exist
expected_studies = ['Huang', 'Liang', 'Liu', 'Montero']
for col in expected_studies:
    if col not in df.columns:
        df[col] = '-'

# Reorder columns for consistency
cols = ['ENSG_ID', 'Gene_Name', 'Probability_Functional'] + expected_studies
# Ensure columns exist before selecting
cols = [c for c in cols if c in df.columns]
df = df[cols]

print(f"Final Processed Matrix DataFrame: {len(df)} unique genes.")

# ---------------------------------------------------------
# Preprocessing
# ---------------------------------------------------------
# Mapping for numeric conversion
# -: 0, Rare: 1, Common: 2, Core: 3
mapping = {'-': 0, 'Rare': 1, 'Common': 2, 'Core': 3}
reverse_mapping = {0: 'None', 1: 'Rare', 2: 'Common', 3: 'Core'}

# Create numeric columns for the studies
study_cols = ['Huang', 'Liang', 'Liu', 'Montero']
df_numeric = df.copy()

for col in study_cols:
    df_numeric[col] = df_numeric[col].map(mapping).fillna(0)

# Ensure Probability is numeric
df['Probability_Functional'] = pd.to_numeric(df['Probability_Functional'], errors='coerce')

# ---------------------------------------------------------
# Plot 1: Clustered Heatmap
# ---------------------------------------------------------
print("Generating Clustered Heatmap...")
# Prepare data for heatmap
heatmap_data = df_numeric[study_cols]
row_colors = df['Probability_Functional']

# Create a custom colormap
# 0: Grey, 1: Yellow, 2: Orange, 3: Red
cmap = mcolors.ListedColormap(['#f0f0f0', '#fed976', '#fd8d3c', '#e31a1c'])
bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
norm = mcolors.BoundaryNorm(bounds, cmap.N)

# Map probability to colors for annotation bar
# Normalize probability 0-1 to a colormap (e.g., Viridis)
prob_cmap = plt.get_cmap("viridis")
prob_norm = plt.Normalize(vmin=0, vmax=1)
prob_colors = df['Probability_Functional'].map(lambda x: prob_cmap(prob_norm(x)) if pd.notnull(x) else (0.8, 0.8, 0.8, 1.0))

# Plot
g = sns.clustermap(heatmap_data, 
                   row_colors=prob_colors, 
                   cmap=cmap, 
                   norm=norm,
                   figsize=(10, 12),
                   dendrogram_ratio=(0.1, 0.2),
                   cbar_pos=(0.02, 0.8, 0.03, 0.15))

# Custom Legend for Essentiality
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#f0f0f0', edgecolor='k', label='None'),
    Patch(facecolor='#fed976', edgecolor='k', label='Rare'),
    Patch(facecolor='#fd8d3c', edgecolor='k', label='Common'),
    Patch(facecolor='#e31a1c', edgecolor='k', label='Core')
]
g.ax_heatmap.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.2, 1), title="Essentiality")

plt.savefig(f"{OUTPUT_DIR}/essentiality_heatmap.png", dpi=300)
plt.close()

# ---------------------------------------------------------
# Plot 2: Probability vs. Essentiality (Boxplot)
# ---------------------------------------------------------
print("Generating Probability vs. Essentiality Boxplot...")
# Use long_df which mirrors R script's data structure
# Filter out "-" (None) if any (though long_df shouldn't have them from the filter)
plot_data = long_df[long_df['Essentiality'] != '-'].copy()

# Define order
order = ['Rare', 'Common', 'Core']

plt.figure(figsize=(12, 6))
sns.boxplot(data=plot_data, x='Essentiality', y='Probability_Functional', hue='Study', order=order, palette="Set2")
plt.title("Functional Probability Distribution by Essentiality Class")
plt.ylabel("Probability Functional")
plt.xlabel("Essentiality Classification")
plt.legend(title='Study')
plt.savefig(f"{OUTPUT_DIR}/probability_vs_essentiality_boxplot.png", dpi=300)
plt.close()

# Pooled version
plt.figure(figsize=(8, 6))
sns.violinplot(data=plot_data, x='Essentiality', y='Probability_Functional', order=order, hue='Essentiality', palette="OrRd", cut=0)
sns.stripplot(data=plot_data, x='Essentiality', y='Probability_Functional', order=order, color='black', alpha=0.3, size=3)
plt.ylim(0, 1)
plt.title("Functional Probability Distribution (Pooled Studies)")
plt.ylabel("Probability Functional")
plt.xlabel("Essentiality Classification")
plt.savefig(f"{OUTPUT_DIR}/probability_vs_essentiality_violin_pooled.png", dpi=300)
plt.close()

# ---------------------------------------------------------
# Plot 3: Study Correlation
# ---------------------------------------------------------
print("Generating Study Correlation Heatmap...")
corr_matrix = df_numeric[study_cols].corr(method='spearman')

plt.figure(figsize=(8, 6))
sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', vmin=-1, vmax=1, center=0)
plt.title("Spearman Correlation of Essentiality Scores\n(0=None, 1=Rare, 2=Common, 3=Core)")
plt.savefig(f"{OUTPUT_DIR}/study_correlation.png", dpi=300)
plt.close()

print(f"Plots saved to {OUTPUT_DIR}")

# ---------------------------------------------------------
# Plot 4: UpSet Plot
# ---------------------------------------------------------
try:
    from upsetplot import UpSet
    print("Generating UpSet Plot...")
    
    # Prepare data: Convert to boolean (True if Essential, False if None)
    upset_data = df[study_cols].applymap(lambda x: x != '-')
    
    # Create MultiIndex Series
    # We need to group by all columns to get counts for each combination
    upset_series = upset_data.groupby(list(upset_data.columns)).size()
    
    plt.figure(figsize=(10, 6))
    # Suppress FutureWarnings from upsetplot library
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=FutureWarning, module='upsetplot')
        upset = UpSet(upset_series, subset_size='count', show_counts=True, sort_by='cardinality')
        upset.plot()
    plt.savefig(f"{OUTPUT_DIR}/upset_plot.png", dpi=300)
    plt.close()
    print("UpSet plot saved.")

except ImportError:
    print("upsetplot library not found. Skipping UpSet plot.")
except Exception as e:
    print(f"Error generating UpSet plot: {e}")

# ---------------------------------------------------------
# Plot 5: Consensus vs. Probability
# ---------------------------------------------------------
print("Generating Consensus vs. Probability Plot...")
# Calculate Consensus from long_df (R-logic: Is the gene essential in the study?)
# A gene is essential in a study if it has ANY entry that is Rare, Common, or Core.
essential_genes = long_df[long_df['Essentiality'].isin(['Rare', 'Common', 'Core'])]
# Count unique studies per gene that have an essential entry
consensus_counts = essential_genes.groupby('ENSG_ID')['Study'].nunique().to_frame('Consensus_Count')

# Merge consensus count back to the FULL long_df
# This allows us to plot the probability of ALL transcripts (even Non-essential ones) 
# against the gene's consensus score.
long_df_consensus = long_df.merge(consensus_counts, on='ENSG_ID', how='left')
long_df_consensus['Consensus_Count'] = long_df_consensus['Consensus_Count'].fillna(0).astype(int)

plt.figure(figsize=(8, 6))
sns.boxplot(data=long_df_consensus, x='Consensus_Count', y='Probability_Functional', hue='Consensus_Count', palette="viridis")
plt.title("Functional Probability by Consensus Level\n(Plotting all transcript entries)")
plt.xlabel("Number of Studies Identifying Gene as Essential")
plt.ylabel("Probability Functional")
plt.savefig(f"{OUTPUT_DIR}/consensus_vs_probability.png", dpi=300)
plt.close()

# ---------------------------------------------------------
# Plot 6: Parallel Categories (Simplified)
# ---------------------------------------------------------
print("Generating Parallel Categories Plot...")
from pandas.plotting import parallel_coordinates

# Prepare data: Use numeric scores
# We need to handle the color column. Let's bin probability for coloring.
parallel_df = df_numeric[study_cols].copy()
parallel_df['Prob_Bin'] = pd.cut(df['Probability_Functional'].fillna(0), bins=4, labels=['Low', 'Med-Low', 'Med-High', 'High'])

plt.figure(figsize=(12, 6))
parallel_coordinates(parallel_df, 'Prob_Bin', colormap=plt.get_cmap("coolwarm"), alpha=0.2)
plt.title("Parallel Categories: Essentiality Flow across Studies\n(0=None, 1=Rare, 2=Common, 3=Core)")
plt.ylabel("Essentiality Score")
plt.legend(title="Probability Functional")
plt.savefig(f"{OUTPUT_DIR}/parallel_categories.png", dpi=300)
plt.close()

# ---------------------------------------------------------
# Plot 7: Top 20 Ranking
# ---------------------------------------------------------
print("Generating Top 20 Ranking...")
# Composite Score = Sum of Essentiality Scores + Probability
df['Essentiality_Sum'] = df_numeric[study_cols].sum(axis=1)
df['Composite_Score'] = df['Essentiality_Sum'] + df['Probability_Functional'].fillna(0)

top_20 = df.sort_values('Composite_Score', ascending=False).head(20)

plt.figure(figsize=(10, 8))
sns.barplot(data=top_20, x='Composite_Score', y='ENSG_ID', hue='ENSG_ID', palette="magma")
plt.title("Top 20 'Super-Essential' Candidates\n(Score = Sum of Essentiality Levels + Probability)")
plt.xlabel("Composite Score")
plt.ylabel("Gene ID")
plt.savefig(f"{OUTPUT_DIR}/top_20_ranking.png", dpi=300)
plt.close()

# ---------------------------------------------------------
# Plot 8: Jaccard Similarity Heatmap
# ---------------------------------------------------------
print("Generating Jaccard Similarity Heatmap...")

def jaccard_index(set1, set2):
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union if union > 0 else 0

jaccard_matrix = pd.DataFrame(index=study_cols, columns=study_cols, dtype=float)

# Create sets of essential genes for each study
study_sets = {}
for study in study_cols:
    # Get genes where score > 0 (Rare, Common, or Core)
    study_sets[study] = set(df_numeric[df_numeric[study] > 0].index)

for s1 in study_cols:
    for s2 in study_cols:
        jaccard_matrix.loc[s1, s2] = jaccard_index(study_sets[s1], study_sets[s2])

plt.figure(figsize=(8, 6))
sns.heatmap(jaccard_matrix, annot=True, cmap="YlGnBu", vmin=0, vmax=1)
plt.title("Jaccard Similarity of Essential Gene Sets\n(Overlap of any essentiality)")
plt.savefig(f"{OUTPUT_DIR}/jaccard_similarity.png", dpi=300)
plt.close()

print("All advanced plots generated.")

# ---------------------------------------------------------
# Plot 9: Pairwise Confusion Matrices (Using Long Data)
# ---------------------------------------------------------
print("Generating Pairwise Confusion Matrices (R-logic: All Pairs)...")

categories = ['Non-essential', 'Rare', 'Common', 'Core']
n_studies = len(study_cols)

fig, axes = plt.subplots(n_studies, n_studies, figsize=(20, 20))

for i, study_x in enumerate(study_cols):
    for j, study_y in enumerate(study_cols):
        ax = axes[i, j]
        
        if i >= j:
            ax.axis('off')
            continue
            
        # Get subsets for the two studies from long_df
        sub_x = long_df[long_df['Study'] == study_x][['ENSG_ID', 'Essentiality']]
        sub_y = long_df[long_df['Study'] == study_y][['ENSG_ID', 'Essentiality']]
        
        # Merge on ENSG_ID to get all combinations of classifications for shared genes
        merged = sub_x.merge(sub_y, on='ENSG_ID', suffixes=('_x', '_y'))
        
        # Compute confusion matrix
        cm = pd.crosstab(merged['Essentiality_x'], merged['Essentiality_y'])
        
        # Reindex to ensure all categories are present and in order
        cm = cm.reindex(index=categories, columns=categories, fill_value=0)
        
        # Plot heatmap
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', cbar=False, ax=ax)
        
        ax.set_title(f"{study_x} vs {study_y}")
        ax.set_xlabel(study_y)
        ax.set_ylabel(study_x)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/pairwise_confusion_matrices.png", dpi=300)
plt.close()

# ---------------------------------------------------------
# Plot 10: Aggregated Confusion Matrix (Pooled)
# ---------------------------------------------------------
print("Generating Aggregated Confusion Matrix (R-logic)...")

aggregated_cm = pd.DataFrame(0, index=categories, columns=categories)

for i, study_x in enumerate(study_cols):
    for j, study_y in enumerate(study_cols):
        if i == j:
            continue
            
        # Same merge logic as above
        sub_x = long_df[long_df['Study'] == study_x][['ENSG_ID', 'Essentiality']]
        sub_y = long_df[long_df['Study'] == study_y][['ENSG_ID', 'Essentiality']]
        merged = sub_x.merge(sub_y, on='ENSG_ID', suffixes=('_x', '_y'))
        
        cm = pd.crosstab(merged['Essentiality_x'], merged['Essentiality_y'])
        cm = cm.reindex(index=categories, columns=categories, fill_value=0)
        
        aggregated_cm += cm

plt.figure(figsize=(10, 8))
sns.heatmap(aggregated_cm, annot=True, fmt='d', cmap='Blues')
plt.title("Aggregated Confusion Matrix (All Study Pairs Pooled)")
plt.xlabel("Classification (Study B)")
plt.ylabel("Classification (Study A)")
plt.savefig(f"{OUTPUT_DIR}/aggregated_confusion_matrix.png", dpi=300)
plt.close()

# ---------------------------------------------------------
# Plot 11: Categorized Overlap Heatmap (Grouped by Category)
# ---------------------------------------------------------
print("Generating Categorized Overlap Heatmap...")

# Categories of interest (excluding '-' if not requested, but user asked for Non-essential)
# User requested: common, core, rare, non essential
target_categories = ['Core', 'Common', 'Rare', 'Non-essential']

# Create labels: "Category - Study"
# We group by Category first as requested
# Create labels: "Category - Study"
# We group by Category first as requested
labels = []
# We will build a binary DataFrame where columns are "Category-Study" and rows are Genes
# Value is 1 if the gene has that classification in that study (in long_df)

# Get unique genes from long_df to index the binary matrix
unique_genes = long_df['ENSG_ID'].unique()
binary_df = pd.DataFrame(index=unique_genes)

for cat in target_categories:
    for study in study_cols:
        # Skip Core-Liu and Core-Huang as requested (no data)
        if cat == 'Core' and (study == 'Liu' or study == 'Huang'):
            continue
            
        label = f"{cat}\n{study}"
        labels.append(label)
        
        # Find genes that have this (Study, Category) combination in long_df
        # This allows a gene to be marked as both "Rare" and "Non-essential" for the same study
        genes_in_group = set(long_df[
            (long_df['Study'] == study) & 
            (long_df['Essentiality'] == cat)
        ]['ENSG_ID'])
        
        # Create boolean series
        binary_df[label] = binary_df.index.isin(genes_in_group).astype(int)

# binary_df is now populated using long_df logic

# Compute Spearman correlation matrix
# This correlates the binary membership vectors
correlation_matrix = binary_df.corr(method='spearman')

plt.figure(figsize=(16, 14))

# Plot heatmap of correlation coefficients
# Generate a mask for the upper triangle, but keep diagonal (k=1 hides strictly above diagonal)
mask = np.triu(np.ones_like(correlation_matrix, dtype=bool), k=1)

sns.heatmap(correlation_matrix, annot=True, fmt='.2f', cmap='coolwarm', vmin=-1, vmax=1, annot_kws={"size": 8}, mask=mask)

# Add lines to separate the major categories
# Dynamic detection of boundaries since group sizes might vary
for i in range(1, len(labels)):
    prev_cat = labels[i-1].split('\n')[0]
    curr_cat = labels[i].split('\n')[0]
    if prev_cat != curr_cat:
        plt.axhline(i, color='white', linewidth=2)
        plt.axvline(i, color='white', linewidth=2)

plt.title("Categorized Spearman Correlation Matrix\n(Grouped by Essentiality Category)")
plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/categorized_overlap_heatmap.png", dpi=300)
plt.close()

# ---------------------------------------------------------
# Plot 12: Categorized Overlap Heatmap (Counts) - Non-aggregated
# ---------------------------------------------------------
print("Generating Categorized Overlap Heatmap (Counts)...")

# Reuse the labels and binary_cols from Plot 11
# We need to compute the overlap matrix (counts)
n_labels = len(labels)
overlap_matrix = pd.DataFrame(0, index=labels, columns=labels)

# We can compute this efficiently using the binary_df
# The overlap between Col A and Col B is the dot product of their binary vectors
# binary_df is (Genes x Labels)
# Overlap = binary_df.T dot binary_df
overlap_matrix = binary_df.T.dot(binary_df)

plt.figure(figsize=(16, 14))

# Log transform for visualization color scale (add 1 to avoid log(0))
log_overlap_matrix = np.log1p(overlap_matrix)

# Use the same mask as before (Lower triangle + Diagonal)
mask = np.triu(np.ones_like(overlap_matrix, dtype=bool), k=1)

sns.heatmap(log_overlap_matrix, annot=overlap_matrix, fmt='d', cmap='YlGnBu', annot_kws={"size": 8}, mask=mask)

# Add lines to separate the major categories
for i in range(1, len(labels)):
    prev_cat = labels[i-1].split('\n')[0]
    curr_cat = labels[i].split('\n')[0]
    if prev_cat != curr_cat:
        plt.axhline(i, color='white', linewidth=2)
        plt.axvline(i, color='white', linewidth=2)

plt.title("Categorized Overlap Matrix (Counts)\n(Grouped by Essentiality Category, Log Color Scale)")
plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/categorized_overlap_heatmap_counts.png", dpi=300)
plt.close()

