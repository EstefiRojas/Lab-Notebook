import pandas as pd
import numpy as np
from scipy.stats import ks_2samp
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from sklearn.preprocessing import RobustScaler
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer

select_datasets1 = ['protein-coding-exon2','protein-coding-exon3','protein-exon2-negative-control','protein-exon3-negative-control']
select_datasets2 = ['lncrna-exon1','lncrna-exon2','lncrna-exon1-negative-control','lcnrna-exon2-negative-control']
select_datasets3 = ['short-ncrna','short-ncrna','short-ncrna-negative-control','short-ncrna-negative-control']

title1 = 'Protein Coding'
title2 = 'LncRNA'
title3 = 'Short ncRNA'

# Data preparation
data = pd.read_csv('../data/features/gene_functionality_features_latest1000all.csv')
#print(data)

# Features of interest
select_features = ["Random_number",
                   "GC_percentage",
                   "CpG","GA","GG","TA","TC",
                   "lowComplexity_density","phyloP_max_241w","phyloP_max_100w",
                   "RPKM_tissue","RPKM_primary.cell","copy_number","Dfam_sum","RNAcode_score","Max_covariance",
                   "MFE","RNAalifold_score","Interaction_ave","gnomAD_SNP_density","gnomAD_MAF_avg",
                   "H3K27ac_AvgSignal","H3K36me3_AvgSignal","H3K79me2_AvgSignal","chrm_acc_AvgSignal","methylome"]

def get_data(select_datasets):
    # Separate the data into functional and control groups
    functional_data = data[(data['Functional'] == 1) & 
                           ((data['Dataset'] == select_datasets[0]) | 
                           (data['Dataset'] == select_datasets[1]))]
    functional_data = functional_data.drop(["Dataset","Functional"], axis=1)
    functional_data = functional_data.dropna() # Drop rows with missing values
    functional_data = functional_data[select_features]
    #print(functional_data)

    control_data = data[(data['Functional'] == 0) & 
                        ((data['Dataset'] == select_datasets[2]) | 
                        (data['Dataset'] == select_datasets[3]))]
    control_data = control_data.drop(["Dataset","Functional"], axis=1)
    control_data = control_data.dropna()
    control_data = control_data[select_features]
    #print(control_data)
    return functional_data, control_data


functional_data1, control_data1 = get_data(select_datasets1)
functional_data2, control_data2 = get_data(select_datasets2)
functional_data3, control_data3 = get_data(select_datasets3)


def compute_ks_test(functional_data, control_data):
    # List to store KS test results
    ks_results = []
    # Perform KS test for each feature
    features = functional_data.columns[:]  # Exclude non-numeric columns
    for feature in features:
        functional_values = functional_data[feature]
        control_values = control_data[feature]
        ks_statistic, p_value = ks_2samp(functional_values, control_values)
        ks_results.append({'Feature': feature, 'KS Statistic': ks_statistic, 'P-value': p_value})

    # Convert results to DataFrame
    return pd.DataFrame(ks_results)


ks_results_df1 = compute_ks_test(functional_data1, control_data1)
ks_results_df2 = compute_ks_test(functional_data2, control_data2)
ks_results_df3 = compute_ks_test(functional_data3, control_data3)


def compute_significant_column(pval, dsat):
    if(pval < 0.05):
        return dsat
    else:
        return None


# Add a column to indicate significance (e.g., p-value < 0.05)
ks_results_df1['D-statistic'] = ks_results_df1['P-value'] < 0.05
ks_results_df2['D-statistic'] = ks_results_df2['P-value'] < 0.05
ks_results_df3['D-statistic'] = ks_results_df3['P-value'] < 0.05

# Apply the function to create a new 'Area' column
ks_results_df1['D-statistic'] = ks_results_df1.apply(lambda row: compute_significant_column(row['P-value'], row['KS Statistic']), axis=1)
ks_results_df2['D-statistic'] = ks_results_df2.apply(lambda row: compute_significant_column(row['P-value'], row['KS Statistic']), axis=1)
ks_results_df3['D-statistic'] = ks_results_df3.apply(lambda row: compute_significant_column(row['P-value'], row['KS Statistic']), axis=1)
#print(ks_results_df)

# Pivot the results DataFrame for heatmap plotting (reshape to long format for visualization)
#ks_heatmap_data = ks_results_df.pivot(index='Feature', columns='Significant', values='KS Statistic')
ks_heatmap_data1 = ks_results_df1[['Feature', 'D-statistic']].set_index('Feature')
ks_heatmap_data2 = ks_results_df2[['Feature', 'D-statistic']].set_index('Feature')
ks_heatmap_data3 = ks_results_df3[['Feature', 'D-statistic']].set_index('Feature')


# Create a custom colormap
colors = ["pink", "red", "orange", "gold", "green"]
n_bins = 1000  # Discretizes the interpolation into bins
cmap_name = 'custom_heatmap'
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)


# Function to map values to colors
def value_to_color(val, vmin=0, vmax=0.9):
    norm = plt.Normalize(vmin, vmax)
    return cm(norm(val))

# Function to annotate the heatmap with colored text
def annotate_heatmap(ax, data, vmin=0, vmax=0.9, **kwargs):
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data.iloc[i, j]
            color = value_to_color(val, vmin, vmax)
            # Set white background
            ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=True, color='white', ec='white'))
            # Add text annotation
            ax.text(j + 0.5, i + 0.5, f'{val:.2f}', ha='center', va='center', color=color)


# Create a figure with three subplots side by side
fig, axs = plt.subplots(1, 3, figsize=(36, 8), sharey=True)

# Create a heatmap of the KS statistics
#plt.figure(figsize=(12, 8))
sns.heatmap(ks_heatmap_data1, annot=False, cmap=cm, vmin=0, vmax=0.9, ax=axs[0], cbar=False)
annotate_heatmap(axs[0], ks_heatmap_data1)
axs[0].set_title(title1)
axs[0].set_xlabel('')
axs[0].set_ylabel('Genomic Features')

sns.heatmap(ks_heatmap_data2, annot=False, cmap=cm, vmin=0, vmax=0.9, ax=axs[1], cbar=False)
annotate_heatmap(axs[1], ks_heatmap_data2)
axs[1].set_title(title2)
axs[1].set_xlabel('')
axs[1].set_ylabel('')

sns.heatmap(ks_heatmap_data3, annot=False, cmap=cm, vmin=0, vmax=0.9, ax=axs[2])
annotate_heatmap(axs[2], ks_heatmap_data3)
axs[2].set_title(title3)
axs[2].set_xlabel('')
axs[2].set_ylabel('')

# Show plots
plt.show()
