
import pandas as pd
import sys
import os

# Determine absolute path relative to this script
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(script_dir)
input_file = os.path.join(project_root, "results", "annotated_unified_genome_alignments.csv")
output_file = input_file

print(f"Target File: {input_file}")

if not os.path.exists(input_file):
    print(f"Error: {input_file} not found.")
    sys.exit(1)

print(f"Reading data...")
try:
    df = pd.read_csv(input_file, low_memory=False)
except Exception as e:
    print(f"Error reading csv: {e}")
    sys.exit(1)

print(f"Initial row count: {len(df)}")

# Define hierarchy
# Higher number = Higher priority
rank_map = {
    'Core': 4,
    'Common': 3,
    'Rare': 2,
    'Non-essential': 1
}

# Function to map rank, handling unknown values as 0 (lowest)
def get_rank(val):
    return rank_map.get(val, 0)

if 'Essentiality' not in df.columns:
    print("Error: Column 'Essentiality' not found in dataframe.")
    sys.exit(1)

if 'Study' not in df.columns or 'ENSG_ID' not in df.columns:
     print("Error: Required columns 'Study' or 'ENSG_ID' not found.")
     sys.exit(1)

df['Rank'] = df['Essentiality'].apply(get_rank)

# Find max rank per (Study, ENSG_ID)
print("Calculating max rank per (Study, ENSG_ID)...")
# We use transform to broadcast the max rank back to the original index
df['Max_Rank'] = df.groupby(['Study', 'ENSG_ID'])['Rank'].transform('max')

# Filter: Keep rows where the row's Rank equals the Max_Rank for that group
# This preserves all gRNAs that match the highest found category
df_filtered = df[df['Rank'] == df['Max_Rank']].copy()

final_count = len(df_filtered)
dropped_count = len(df) - final_count

print(f"Filtered row count: {final_count}")
print(f"Dropped {dropped_count} rows derived from lower-priority classifications.")

# Clean up temporary columns
df_filtered = df_filtered.drop(columns=['Rank', 'Max_Rank'])

print(f"Saving to {output_file}...")
df_filtered.to_csv(output_file, index=False)
print("Filtering complete.")
