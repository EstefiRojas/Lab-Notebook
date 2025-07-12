load_epigenetic_zscores <- function(subset = FALSE) {
  source("scripts/config.R")
  library(dplyr)
  library(tidyr)

  if (subset) {
    # Load subset zscores files
    funcProtZscores <- data.frame(read.csv(MRNA_ZSCORES_SUBSET_FILE, header=TRUE))
    funcLncrnaZscores <- data.frame(read.csv(LNCRNA_ZSCORES_SUBSET_FILE, header=TRUE))
    funcSncrnaZscores <- data.frame(read.csv(SNCRNA_ZSCORES_SUBSET_FILE, header=TRUE))

    protNCZscores <- data.frame(read.csv(MRNA_NC_ZSCORES_SUBSET_FILE, header=TRUE))
    lncrnaNCZscores <- data.frame(read.csv(LNCRNA_NC_ZSCORES_SUBSET_FILE, header=TRUE))
    sncrnaNCZscores <- data.frame(read.csv(SNCRNA_NC_ZSCORES_SUBSET_FILE, header=TRUE))
  } else {
    # Load all zscores files
    funcProtZscores <- data.frame(read.csv(MRNA_ZSCORES_FILE, header=TRUE))
    funcLncrnaZscores <- data.frame(read.csv(LNCRNA_ZSCORES_FILE, header=TRUE))
    funcSncrnaZscores <- data.frame(read.csv(SNCRNA_ZSCORES_FILE, header=TRUE))

    protNCZscores <- data.frame(read.csv(MRNA_NC_ZSCORES_FILE, header=TRUE))
    lncrnaNCZscores <- data.frame(read.csv(LNCRNA_NC_ZSCORES_FILE, header=TRUE))
    sncrnaNCZscores <- data.frame(read.csv(SNCRNA_NC_ZSCORES_FILE, header=TRUE))
  }

  # Join in a single dataframe
  zscores <- rbind(data.frame(funcProtZscores), 
                       data.frame(funcLncrnaZscores), 
                       data.frame(funcSncrnaZscores),
                       
                       data.frame(protNCZscores), 
                       data.frame(lncrnaNCZscores), 
                       data.frame(sncrnaNCZscores))

  return(zscores)
}