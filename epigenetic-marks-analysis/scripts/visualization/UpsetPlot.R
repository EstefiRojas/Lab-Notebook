# Install necessary packages (if not already installed)
if (!requireNamespace("UpSetR", quietly = TRUE)) install.packages("UpSetR")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

# Load required libraries
library(dplyr)
library(UpSetR)

# Define file paths
file_paths <- list(
  gencodev44 = "../data/UpSet/chr22/result/gencodev44.bed",
  gencodev45 = "../data/UpSet/chr22/result/gencodev45.bed",
  hgnc = "../data/UpSet/chr22/result/hgnc.bed",
  ncbi = "../data/UpSet/chr22/result/ncbi.bed",
  rnacentral = "../data/UpSet/chr22/result/rnacentral.bed",
  uniprot = "../data/UpSet/chr22/result/uniprot.bed"
)

# Function to read BED files
read_bed_file <- function(file_path) {
  read.csv(file_path, header = TRUE, sep = "\t")
}

# Read all BED files into a list of data frames
data_frames <- lapply(file_paths, read_bed_file)

# Combine all data frames into one
all_annotations <- bind_rows(data_frames)

# Provide summary of the combined data
summary(all_annotations)

# Convert 'Start' and 'End' columns to numeric
all_annotations <- all_annotations %>%
  mutate(
    Start = as.numeric(Start),
    End = as.numeric(End)
  )

# Check the class of the 'Start' column
print(class(all_annotations$Start))

# Rename the 'NCBI' column to 'RefSeq' if it exists
if ("NCBI" %in% colnames(all_annotations)) {
  names(all_annotations)[names(all_annotations) == 'NCBI'] <- 'RefSeq'
}

# Generate UpSet plot
upset(all_annotations, sets = c("GencodeV44", "GencodeV45", "HGNC", "RefSeq", "RNACentral", "UniProt"), 
      order.by = "freq")
