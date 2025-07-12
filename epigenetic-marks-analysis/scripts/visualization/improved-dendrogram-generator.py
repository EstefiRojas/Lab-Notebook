import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage

# Constants
DATA_PATH = "../data/features/gene_functionality_features_latest1000all_new_din.csv"
FEATURES = [
    "Random number",
    "GC content",
    "AA","AC","AG","AT","CA","CC","CpG","CT","GA","GC","GG","GT","TA","TC","TG","TT",
    "low_complexity_density","Conservation","phyloP max_100w","GERP_91_mammals_max","GERP_63_amniotes_max",
    "Expression","RPKM_primary cell","Copy number","Repeat free","Fickett_score","RNAcode","RNA Covariance",
    "MFE","accessibility","RNAalifold","RNA-RNA interactions","SNP density","gnomAD_MAF",
    "H3K27ac","H3K36me3","H3K79me2","chromatin_acc","methylome"
]

def load_data(file_path):
    """
    Load data from CSV and select specified features.
    """
    return pd.read_csv(file_path, sep=",")[FEATURES]
    #return data[FEATURES]

def calculate_correlation_matrix(df):
    """
    Calculate the correlation matrix for the given dataframe.
    """
    return df.corr().abs()

def calculate_distance_matrix(correlation_matrix):
    """
    Convert correlation matrix to distance matrix.
    """
    return 1 - correlation_matrix

def perform_hierarchical_clustering(distance_matrix):
    """
    Perform hierarchical clustering on the distance matrix.
    """
    linkage_matrix = linkage(distance_matrix, 'ward')
    plot_dendrogram(linkage_matrix, distance_matrix.columns)
    return distance_matrix

def plot_dendrogram(linkage_matrix, labels):
    """
    Plot dendrogram based on the linkage matrix.
    """
    #labels = linkage_matrix.columns
    plt.figure(figsize=(12, 8))
    dendrogram(linkage_matrix, labels=labels, leaf_rotation=90, leaf_font_size=8)
    plt.title('Dendrogram Based on Feature Correlation for all gene types', fontsize=18)
    plt.xlabel('Features', fontsize=14)
    plt.ylabel('Distance', fontsize=14)
    plt.tight_layout()
    plt.show()

def main():
    
    #df = load_data(DATA_PATH, FEATURES)
    
    
    #corr_matrix = calculate_correlation_matrix(df)
    #dist_matrix = calculate_distance_matrix(corr_matrix)
    
    # Perform hierarchical clustering
    #linkage_matrix = perform_hierarchical_clustering(dist_matrix)
    
    # Plot dendrogram
    #plot_dendrogram(linkage_matrix, corr_matrix.columns)

    df = (pd.DataFrame()
            .pipe(lambda x: load_data(DATA_PATH)) # Load and preprocess data
            .pipe(calculate_correlation_matrix) # Calculate correlation matrix
            .pipe(calculate_distance_matrix) # Calculate distance matrix
            .pipe(perform_hierarchical_clustering) # Calculate correlation and distance matrices
           )

if __name__ == "__main__":
    main()