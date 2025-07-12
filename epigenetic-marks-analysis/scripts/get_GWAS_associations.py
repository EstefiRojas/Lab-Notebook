#!/usr/bin/env python3
"""
This script retrieves the minimum p-value of SNPs present at lncRNAs from 
GWAS summary-statistics API using corresponding ENSEMBL coordinates as filter.
It also filters by p-value < 1e-5.

Dependencies: 
- requests
- pandas
- tqdm

Author: Estefania Rojas (Python version)

Input: 
1. Csv file: containing the dataset of sequences to be analysed. This 
   file should have the exon 1 and 2 coordinates at columns 10-21.
   Ex. ../data/lncRNAs_ensembl_exon_coords.csv
2. Name of the output csv file: minimum p-value for each coordinate in input file 1.
   Ex. lncrna-gwas-pval-feature
"""

import sys
import os
from datetime import datetime
import pandas as pd
import requests
from tqdm import tqdm
import json

def convert_chromosome(chrm):
    """Convert X to 23 and Y to 24 while leaving numeric values unchanged."""
    if pd.isna(chrm):  # Handle NA/NaN values
        return None
    chrm = str(chrm).upper()
    if chrm == 'X':
        return '23'
    elif chrm == 'Y':
        return '24'
    return chrm

def is_numeric(value):
    """Check if a value is numeric."""
    if pd.isna(value):  # Handle NA/NaN values
        return False
    try:
        float(value)
        return True
    except (ValueError, TypeError):
        return False

def is_valid_chromosome(chrm):
    """Validate chromosome format."""
    if pd.isna(chrm):  # Handle NA/NaN values
        return False
    chrm = str(chrm).upper().replace('CHR', '')
    valid_chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y', 'M', 'MT']
    return chrm in valid_chromosomes

def process_line(line):
    """Process a line of data from the CSV file."""
    # Ensure all fields are present
    fields = line + [''] * (21 - len(line))  # Pad with empty strings if needed
    
    # Extract fields
    coordinates = {
        'tl_exon1': {'chrm': fields[9], 'start': fields[10], 'end': fields[11]},
        'tl_exon2': {'chrm': fields[12], 'start': fields[13], 'end': fields[14]},
        'tl': {'chrm': fields[15], 'start': fields[16], 'end': fields[17]},
        'tr': {'chrm': fields[18], 'start': fields[19], 'end': fields[20]}
    }
    
    # Validate coordinates
    for region, data in coordinates.items():
        if data['chrm'] and not is_valid_chromosome(data['chrm']):
            print(f"Warning: Invalid chromosome format in {region}: {data['chrm']}", file=sys.stderr)
        
        for pos_type in ['start', 'end']:
            if data[pos_type] and not is_numeric(data[pos_type]):
                print(f"Warning: Invalid position in {region} {pos_type}: {data[pos_type]}", file=sys.stderr)
    
    return coordinates

def get_unique_gwas_variants(coordinates):
    """Get unique GWAS variants for the given coordinates."""
    unique_variants = set()
    base_url = "https://gwas.mrcieu.ac.uk/api/chromosomes/"
    
    for region, data in coordinates.items():
        # Skip if any required data is missing or NA
        if pd.isna(data['chrm']) or pd.isna(data['start']) or pd.isna(data['end']):
            continue
            
        try:
            chrm = convert_chromosome(str(data['chrm']).replace('chr', ''))
            if chrm is None:  # Skip if chromosome conversion failed
                continue
                
            response = requests.get(f"{base_url}{chrm}/{data['start']}/{data['end']}")
            if response.status_code == 200:
                variants = response.json()
                if variants:
                    unique_variants.update(variants)
        except requests.exceptions.RequestException as e:
            print(f"Error fetching variants for {region}: {e}", file=sys.stderr)
        except (AttributeError, ValueError) as e:
            print(f"Error processing data for {region}: {e}", file=sys.stderr)
            continue
    
    return list(unique_variants)

def main():
    # Check command line arguments
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} regions_to_extract.csv outname.csv")
        sys.exit(1)

    regions_file = sys.argv[1]
    output_name = sys.argv[2]

    # Print current date
    current_date = datetime.now().strftime("%A %B %d, %Y")
    print(current_date)

    # Read input file
    try:
        df = pd.read_csv(regions_file)
    except Exception as e:
        print(f"Error reading input file: {e}", file=sys.stderr)
        sys.exit(1)

    # Process each line and get GWAS variants
    results = []
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Processing regions"):
        coordinates = process_line(row.tolist())
        variants = get_unique_gwas_variants(coordinates)
        
        # Process variants and find minimum p-value
        min_pval = float('inf')
        for variant in variants:
            try:
                pval = float(variant.get('p', float('inf')))
                if pval < min_pval and pval < 1e-5:
                    min_pval = pval
            except (ValueError, TypeError):
                continue
        
        results.append({
            'coordinates': coordinates,
            'min_pval': min_pval if min_pval != float('inf') else None
        })

    # Create output DataFrame and save to CSV
    output_df = pd.DataFrame(results)
    output_df.to_csv(output_name, index=False)
    print(f"Results saved to {output_name}")

if __name__ == "__main__":
    main()
