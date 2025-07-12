import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


def remove_outliers_IQR(df, col):
    Q1 = df[col].quantile(0.25)
    Q3 = df[col].quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR
    return df[(df[col] >= lower_bound) & (df[col] <= upper_bound)][col]

def mad(x, nan_policy='omit'):
    return stats.median_abs_deviation(x, nan_policy=nan_policy)


def get_robust_zscores(functional_df, negative_df, selected_features):
    zscores_dict = {}
    methods_dict = {}
    
    for col in selected_features:
        positive_col = remove_outliers_IQR(functional_df, col)
        negative_col = remove_outliers_IQR(negative_df, col)
        
        if len(positive_col) == 0:
            positive_col = functional_df[col]
        if len(negative_col) == 0:
            negative_col = negative_df[col]
        
        mad_value = mad(negative_col, nan_policy='omit')
        median_value = np.median(negative_col)
        
        mean_value = np.mean(negative_col)
        absolute_deviations = np.abs(negative_col - mean_value)
        meanAD_value = np.mean(absolute_deviations)
        
        if mad_value != 0:  # MAD method
            zscores_dict[col] = (positive_col - median_value) / mad_value
            methods_dict[col] = "MAD"
        elif meanAD_value != 0:  # meanAD estimator
            zscores_dict[col] = (positive_col - median_value) / (1.2533 * meanAD_value)
            methods_dict[col] = "meanAD"
        else:  # regular z-score
            zscores_dict[col] = (positive_col - np.mean(negative_col)) / np.std(negative_col, ddof=1)
            methods_dict[col] = "regular z-score"
    
    # Create a DataFrame from the zscores dictionary
    zscores_df = pd.DataFrame(zscores_dict)
    
    # Add a 'method' column to store the method used for each feature
    zscores_df['method'] = pd.Series(methods_dict)
    
    return zscores_df


def calculate_delta_zscores(zscores_df1, zscores_df2):
    """
    Calculate the delta z-scores between two z-score DataFrames.
    
    Parameters:
    zscores_df1 (pd.DataFrame): First z-score DataFrame
    zscores_df2 (pd.DataFrame): Second z-score DataFrame
    
    Returns:
    pd.DataFrame: Delta z-scores
    """
    # Ensure both DataFrames have the same columns
    common_columns = zscores_df1.columns.intersection(zscores_df2.columns)
    common_columns = [col for col in common_columns if col != 'method']
    
    if len(common_columns) == 0:
        raise ValueError("No common features found between the two DataFrames")
    
    # Calculate delta z-scores
    delta_zscores = zscores_df2[common_columns] - zscores_df1[common_columns]
    
    # Add method columns
    delta_zscores['method_1'] = zscores_df1['method']
    delta_zscores['method_2'] = zscores_df2['method']
    
    return delta_zscores


def calculate_binned_delta_zscores(df, bin_feature, bin_edges, labels, selected_features):
    """
    Calculate delta z-scores based on binning a specific feature.
    
    Parameters:
    df (pd.DataFrame): Input DataFrame containing all features
    bin_feature (str): Name of the feature to be binned (e.g., 'DistanceGene')
    bin_edges (list): List of bin edges for the binning process
    labels (list): List of names to give the different bin slots
    selected_features (list): List of features to calculate z-scores for
    
    Returns:
    pd.DataFrame: Delta z-scores for each bin compared to the first bin
    """
    # Create a new column for distance bins
    df['bin'] = pd.cut(df[bin_feature], bins=bins, labels=labels, include_lowest=True)
    #df['bin'] = pd.cut(df[bin_feature], bins=bin_edges, labels=False)
    df['bin'].to_csv('../data/delta_zscores/bin_cut.csv')
    
    # Function to calculate z-scores for a group
    def group_zscores(group):
        functional_df = group[group['Functional'] == 1]
        negative_df = group[group['Functional'] == 0]
        return get_robust_zscores(functional_df, negative_df, selected_features)
    
    # Calculate z-scores for each bin
    bin_zscores = df.groupby('bin').apply(group_zscores)
    bin_zscores.to_csv('../data/delta_zscores/bin_zscores.csv')
    #pd.write_csv(bin_zscores,'../data/delta_zscores/bin_zscores.csv')
    
    # Calculate delta z-scores (comparing each bin to the first bin)
    #first_bin_zscores = bin_zscores.iloc[0]
    first_bin_zscores = group_zscores(df) # Calculate positives z-scores against all negatives.
    delta_zscores = bin_zscores.subtract(first_bin_zscores)
    
    # Reshape the result for easier interpretation
    delta_zscores_reshaped = delta_zscores.reset_index()
    delta_zscores_reshaped = delta_zscores_reshaped.melt(id_vars=['bin', 'level_1'],
                                                         var_name='feature',
                                                         value_name='delta_zscore')
    delta_zscores_reshaped = delta_zscores_reshaped.pivot(index=['bin', 'feature'],
                                                          columns='level_1',
                                                          values='delta_zscore')
    delta_zscores_reshaped.columns.name = None
    delta_zscores_reshaped = delta_zscores_reshaped.reset_index()
    
    return delta_zscores_reshaped

# 1. Load features
df = pd.read_csv("../data/features/gene_functionality_features_latest1000all_distance.csv", sep=",")

# Features of interest
select_features = ["Functional","Distance","Random number",
                   "GC content"]#,
                   #"CpG","GA","GG","TA","TC",
                   #"low_complexity_density","phyloP max_241w","phyloP max_100w",
                   #RPKM_tissue","RPKM_primary cell","Copy number","Repeat free","RNAcode","Max covariance",
                   #"MFE","RNAalifold","Interaction_ave","gnomAD_SNP_density","gnomAD_MAF",
                   #"H3K27ac","H3K36me3","H3K79me2","chromatin_acc","methylome"]

df = df[select_features]
df = df.dropna()

# 2. Define bins and labels
bins = [-1, 1, 5000, 50000, 500000, 5000000]
labels = ['Functional', '1-5k', '5k-50k', '50k-500k', '500k-5M']
df['Distance'] = pd.to_numeric(df['Distance'])
#df['distance_bin'] = pd.cut(df['Distance'], bins=bins, labels=labels, include_lowest=True)

# Calculate delta z-scores
delta_zscores = calculate_binned_delta_zscores(df, 'Distance', bins, labels, select_features)

# List of feature columns
features = delta_zscores.drop(columns=['Functional','Distance'], axis=1).columns.tolist()

### VIOLIN PLOTS ###
# Create a violin plot for each feature grouped by distance bins
for feature in select_features:
  plt.figure(figsize=(10, 6))
  sns.violinplot(x='bin', y=feature, data=delta_zscores)

  plt.title(f'Violin Plot of {feature} by Distance to Gene')
  plt.xlabel('Distance to Functional')
  plt.ylabel(feature)

  # Save the plot to a file
  plt.savefig(f'../results/delta-zscores/all-{feature}_violin_plot.png')  # Change the filename as needed
  plt.close()  # Close the figure to avoid display
