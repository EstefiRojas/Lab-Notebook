#####################
# Mandatory scripts #
#####################
# Commands needed to obtain epigenetic analysis features
cd scripts/

# Dowloads Marks from GENCODE
./download_epigenetic_data.sh ../listOfCodes.txt # Download and create appropiate folders in data/

# Process histone mark features
# All resulting feature files stored in "data/histone_feature/<histone_name>/<histone_name>_<sequence_type>_matrix.csv"
./run_histone.sh H3k4me1
./run_histone.sh H3k4me3
./run_histone.sh H3k9me3
./run_histone.sh H3k27ac
./run_histone.sh H3k27me3
./run_histone.sh H3k36me3

# Process Chromatin Accessibility mark feature
# All resulting feature files stored in "data/chrm_acc_feature/crhm_acc_<sequence_type>_matrix.csv"
./run_chrm_acc.sh chrm

# Process Methylome marks feature
# All resulting feature files stored in "data/methylome_feature/<sequence_type>_matrix.csv"
./run_methylome.sh
##################



#########################
# Non mandatory scripts #
#########################
# Run Spearman correlation analysis from R file
r histone_spearman_correlation_analysis.R


# Commands used to obtain a HeatMap of epigenetic features
./join_epigenetic_features.sh ../data/chrm_acc_feature/chrm_acc_lncrna_matrix.csv list_of_files_lncrna.txt lncrna_epigenetic_features_matrix.csv
./join_epigenetic_features.sh ../data/chrm_acc_feature/chrm_acc_short_ncrna_matrix.csv list_of_files_short_ncrna.txt short_ncrna_epigenetic_features_matrix.csv
./join_epigenetic_features.sh ../data/chrm_acc_feature/chrm_acc_protein_matrix.csv list_of_files_protein.txt protein_epigenetic_features_matrix.csv

r heat_map.R


# Commands used to format gene annotation data from different databases into UpSet Plot format (Just chr22)
# Preprocess databases:
# Reformat ncbi database
gtf2bed < genomic.gtf > ncbi.bed

# Extract positions from chromosome 22
awk -F '\t' '$1=="chr22" {print $1"\t"$2"\t"$3"\t"$10}' OFS="\t", gencode.v44.annotation.bed > gencode.v44.chr22.annotation.bed
awk -F '\t' '$1=="chr22" {print $1"\t"$2"\t"$3"\t"$10}' OFS="\t", gencode.v45.annotation.bed > gencode.v45.chr22.annotation.bed
awk -F '\t' '$1=="chr22" {print $1"\t"$2"\t"$3"\t"$12}' OFS=, hgnc.bed > hgnc_chr22.bed
awk -F '\t' '$1=="NC_000022.11" {print $1"\t"$2"\t"$3"\t"$8}' OFS="\t", ncbi.bed > NCBI_chr22.bed
awk -F '\t' '$1=="chr22" {print $1"\t"$2"\t"$3"\t"$14}' OFS="\t", rna\ central_homo_sapiens.GRCh38.bed > rnacentral_chr22.bed
awk -F '\t' '$1=="chr22" {print $1"\t"$2"\t"$3}' OFS="\t", uniprotswissprot.bed > uniprot_chr22.bed

# Add a column with name of database
awk -F '\t' '{print $0"\tGencodeV44"}' OFS="\t", chr22/gencode.v44.chr22.genetype.bed > chr22/gencodev44_chr22_database.bed
awk -F '\t' '{print $0"\tGencodeV45"}' OFS="\t", chr22/gencode.v45.chr22.genetype.bed > chr22/gencodev45_chr22_database.bed
awk -F '\t' '{print $0"\tHGNC"}' OFS="\t", chr22/hgnc_chr22.bed > chr22/hgnc_chr22_database.bed
awk -F '\t' '{print $0"\tNCBI"}' OFS="\t", chr22/NCBI_chr22.bed > chr22/NCBI_chr22_database.bed
awk -F '\t' '{print $0"\tRNACentral"}' OFS="\t", chr22/rnacentral_chr22.bed > chr22/rnacentral_chr22_database.bed
awk -F '\t' '{print $0"\tUniprot"}' OFS="\t", chr22/uniprot_chr22.bed > chr22/uniprot_chr22_database.bed

# Ignore genetype column
cut -d $'\t' -f 1-3,5 chr22/hgnc_chr22_database.bed > chr22/hgnc.bed
cut -d $'\t' -f 1-3,5 chr22/NCBI_chr22_database.bed > chr22/ncbi.bed
cut -d $'\t' -f 1-3,5 chr22/uniprot_chr22_database.bed > chr22/uniprot.bed
cut -d $'\t' -f 1-3,5 chr22/rnacentral_chr22_database.bed > chr22/rnacentral.bed
cut -d $'\t' -f 1-3,5 chr22/gencodev44_chr22_database.bed > chr22/gencodev44.bed
cut -d $'\t' -f 1-3,5 chr22/gencodev45_chr22_database.bed > chr22/gencodev45.bed

# Format databases
./upset_formater.sh ../data/UpSet/chr22/gencodev44.bed UpSetDatabaseList.txt ../data/UpSet/chr22/chr22/gencodev44.bed
./upset_formater.sh ../data/UpSet/chr22/gencodev45.bed UpSetDatabaseList.txt ../data/UpSet/chr22/chr22/gencodev45.bed
./upset_formater.sh ../data/UpSet/chr22/hgnc.bed UpSetDatabaseList.txt ../data/UpSet/chr22/chr22/hgnc.bed
./upset_formater.sh ../data/UpSet/chr22/ncbi.bed UpSetDatabaseList.txt ../data/UpSet/chr22/chr22/ncbi.bed
./upset_formater.sh ../data/UpSet/chr22/rnacentral.bed UpSetDatabaseList.txt ../data/UpSet/chr22/chr22/rnacentral.bed
./upset_formater.sh ../data/UpSet/chr22/uniprot.bed UpSetDatabaseList.txt ../data/UpSet/chr22/chr22/uniprot.bed

r UpSetPlot.R
#############