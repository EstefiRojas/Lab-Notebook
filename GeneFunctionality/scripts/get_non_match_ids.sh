# Create the list of all targets that were analyzed
cut -f3 tmp/gRNA_metadata_map.tsv | sort -u > tmp/all_target_ids.txt

# Create the list of successful targets from the final unique report
tail -n +2 results/gRNA_lncRNA_matches_unique_sorted.tsv | awk -F'\t' '$9 != "NA"' | cut -f3 | sort -u > tmp/successful_target_ids.txt

# Find the non-matching targets and save the report
grep -v -F -f tmp/successful_target_ids.txt tmp/all_target_ids.txt > results/non_matching_Target_Genes.txt

# Clean up temporary files
rm tmp/all_target_ids.txt tmp/successful_target_ids.txt

echo "Report of non-matching Target_Gene_IDs created at results/non_matching_Target_Genes.txt"