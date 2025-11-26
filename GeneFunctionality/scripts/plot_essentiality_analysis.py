import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import os

# Define paths
INPUT_FILE = "../results/essentiality_matrix.csv"
OUTPUT_DIR = "../results/plots"

# Create output directory
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

print(f"Loading data from {INPUT_FILE}...")
df = pd.read_csv(INPUT_FILE)

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
# Melt the dataframe
melted_df = df.melt(id_vars=['ENSG_ID', 'Probability_Functional'], 
                    value_vars=study_cols, 
                    var_name='Study', 
                    value_name='Essentiality')

# Filter out "-" (None)
melted_df = melted_df[melted_df['Essentiality'] != '-']

# Define order
order = ['Rare', 'Common', 'Core']

plt.figure(figsize=(12, 6))
sns.boxplot(data=melted_df, x='Essentiality', y='Probability_Functional', hue='Study', order=order, palette="Set2")
plt.title("Functional Probability Distribution by Essentiality Class")
plt.ylabel("Probability Functional")
plt.xlabel("Essentiality Classification")
plt.legend(title='Study')
plt.savefig(f"{OUTPUT_DIR}/probability_vs_essentiality_boxplot.png", dpi=300)
plt.close()

# Pooled version
plt.figure(figsize=(8, 6))
sns.violinplot(data=melted_df, x='Essentiality', y='Probability_Functional', order=order, palette="OrRd", cut=0)
sns.stripplot(data=melted_df, x='Essentiality', y='Probability_Functional', order=order, color='black', alpha=0.3, size=3)
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
    upset_series = upset_data.value_counts()
    
    plt.figure(figsize=(10, 6))
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
# Count how many studies consider the gene essential (score > 0)
df['Consensus_Count'] = (df_numeric[study_cols] > 0).sum(axis=1)

plt.figure(figsize=(8, 6))
sns.boxplot(data=df, x='Consensus_Count', y='Probability_Functional', palette="viridis")
plt.title("Functional Probability by Consensus Level")
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
sns.barplot(data=top_20, x='Composite_Score', y='ENSG_ID', palette="magma")
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

