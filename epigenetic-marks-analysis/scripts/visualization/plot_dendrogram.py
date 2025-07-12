import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage


# Variables
#select_datasets = ['protein-coding-exon2','protein-coding-exon3']
#select_datasets = ['lncrna-exon1','lncrna-exon2']
select_datasets = ['short-ncrna']
#select_datasets = ['protein-coding-exon2','protein-coding-exon3','lncrna-exon1','lncrna-exon2','short-ncrna']

# 1. Load features
data = pd.read_csv("../data/features/gene_functionality_features_latest1000all_new_din.csv", sep=",")

# Features of interest
#select_features = ["Random number",
#                   "GC content",
#                   "CpG","GA","GG","TA","TC",
#                   "low_complexity_density","phyloP max_241w","phyloP max_100w",
#                   "RPKM_tissue","RPKM_primary cell","Copy number","Repeat free","RNAcode","Max covariance",
#                  "MFE","RNAalifold","Interaction_ave","gnomAD_SNP_density","gnomAD_MAF",
#                   "H3K27ac","H3K36me3","H3K79me2","chromatin_acc","methylome"]

select_features = ["Random number",
                     "GC content",
                     "AA","AC","AG","AT","CA","CC","CpG","CT","GA","GC","GG","GT","TA","TC","TG","TT",
                     "low_complexity_density","phyloP max_241w","phyloP max_100w","GERP_91_mammals_max","GERP_63_amniotes_max",
                     "RPKM_tissue","RPKM_primary cell","Copy number","Repeat free","Fickett_score","RNAcode","Max covariance",
                     "MFE","accessibility","RNAalifold","Interaction_ave","gnomAD_SNP_density","gnomAD_MAF",
                     "H3K27ac","H3K36me3","H3K79me2","chromatin_acc","methylome"]

# 1.1 Preprocess data
#data = data.dropna() # Drop rows with missing values
data = data[((data['Dataset'] == select_datasets[0]))] #| 
             #(data['Dataset'] == select_datasets[1]))]


# 2. Features (X)
X = data.drop(["Functional","Dataset"], axis=1)
X = X[select_features]
#X_abs = X.abs()
df = pd.DataFrame(X)

# Calculate the correlation matrix
correlation_matrix = df.corr()
corr_matrix_abs = correlation_matrix.abs()

# Convert the correlation matrix to a distance matrix
distance_matrix = 1 - corr_matrix_abs

# Perform hierarchical clustering
Z = linkage(distance_matrix, 'ward')

# Plot the dendrogram
plt.figure(figsize=(12, 8))
dendrogram(Z, labels=correlation_matrix.columns, leaf_rotation=90, leaf_font_size=14)
plt.title('Dendrogram Based on Feature Correlation - sincRNA', fontsize=18)
plt.xlabel('Features', fontsize=14)
plt.ylabel('Distance', fontsize=14)
plt.show()