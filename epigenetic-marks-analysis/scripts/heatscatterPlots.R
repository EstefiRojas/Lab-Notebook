# Load Libraries
library("dplyr")
library(LSD)

# Load z-scores data #

# Prepare z-score data #
# PROTEIN: #
prot_pos_zscores <- protein_functional_z_scores %>%
  mutate(distance=1)

prot_neg_zscores <- protein_negative_z_scores %>%
  mutate(distance = rbind(
    protExon2NCData %>% dplyr::select(DistanceGene),
    protExon3NCData %>% dplyr::select(DistanceGene)
  )$DistanceGene)

prot_pos_zscores$group <- "mrna(+)"
prot_neg_zscores$group <- "mrna(-)"
prot_zscores_all <- rbind(prot_pos_zscores,prot_neg_zscores)
prot_zscores_all$group <- factor(prot_zscores_all$group, levels = c("mrna(+)","mrna(-)"))


# LNCRNA: #
lncrna_pos_zscores <- lncrna_functional_z_scores %>%
  mutate(distance=1)

lncrna_neg_zscores <- lncrna_negative_z_scores %>%
  mutate(distance = rbind(
    lncrnaExon1NCData %>% dplyr::select(DistanceGene),
    lncrnaExon2NCData %>% dplyr::select(DistanceGene)
  )$DistanceGene)

lncrna_pos_zscores$group <- "lncrna(+)"
lncrna_neg_zscores$group <- "lncrna(-)"
lncrna_zscores_all <- rbind(lncrna_pos_zscores,lncrna_neg_zscores)
lncrna_zscores_all$group <- factor(lncrna_zscores_all$group, levels = c("lncrna(+)","lncrna(-)"))


# SNCRNA: #
sncrna_pos_zscores <- sncrna_functional_z_scores %>%
  mutate(distance=1)

sncrna_neg_zscores <- sncrna_negative_z_scores %>%
  mutate(distance = (sncrnaNCData %>% dplyr::select(DistanceGene)
  )$DistanceGene)

sncrna_pos_zscores$group <- "sncrna(+)"
sncrna_neg_zscores$group <- "sncrna(-)"
sncrna_zscores_all <- rbind(sncrna_pos_zscores,sncrna_neg_zscores)
sncrna_zscores_all$group <- factor(sncrna_zscores_all$group, levels = c("sncrna(+)","sncrna(-)"))


##############################
### USE heatscatter from LSD R PACKAGE.
#reset plot
dev.off()
layout(1)

# For high-res PNG
png("effectOnDistanceJoined5k_5M.png", width=4920, height=3240, res=300)
nf <- layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))
#GC%
par(mar = c(7, 5, 7, 2))  # Increase top margin (third number)
heatscatter(prot_neg_zscores$distance, 
            prot_neg_zscores$GC_percentage,
            daltonize = TRUE, 
            cor = FALSE,
            log = "x",
            xlim = c(7500, 5100000),
            xlab = "Distance to Gene",
            ylab = "Robust Z-score",
            main = "",
            xaxt = 'n',
            cex.axis = 2,
            cex.lab = 2.5)
axis(1, at = c(1000, 10000, 100000, 1000000, 5000000), 
     labels = c("1kb", "10kb", "100kb", "1Mb", "5Mb"),
     cex.axis = 2)
title(main = "mRNA", cex.main = 3)

# Filter out NA values first
valid_data <- !is.na(prot_neg_zscores$distance) & !is.na(prot_neg_zscores$GC_percentage)

# Sort the valid data
sorted_idx <- order(prot_neg_zscores$distance[valid_data])
sorted_distance <- prot_neg_zscores$distance[valid_data][sorted_idx]
sorted_gc <- prot_neg_zscores$GC_percentage[valid_data][sorted_idx]

# Fit LOESS on sorted, valid data
loess_fit <- loess(sorted_gc ~ sorted_distance, span = 0.1)
pred <- predict(loess_fit, se = TRUE)

# Add LOESS curve and confidence intervals
lines(sorted_distance, pred$fit, col = "red", lwd = 3)
lines(sorted_distance, pred$fit + 1.96 * pred$se.fit, col = "grey", lty = 3)
lines(sorted_distance, pred$fit - 1.96 * pred$se.fit, col = "grey", lty = 3)

par(mar = c(7, 5, 7, 2))  # Increase top margin (third number)
heatscatter(sncrna_neg_zscores$distance, 
            sncrna_neg_zscores$GC_percentage,
            daltonize = TRUE, 
            cor = FALSE,
            log = "x",
            xlim = c(7500, 5100000),
            xlab = "Distance to Gene",
            ylab = "Robust Z-score",
            main = "",
            xaxt = 'n',
            cex.axis = 2,
            cex.lab = 2.5)
axis(1, at = c(1000, 10000, 100000, 1000000, 5000000), 
     labels = c("1kb", "10kb", "100kb", "1Mb", "5Mb"),
     cex.axis = 2)
title(main = "sncRNA", cex.main = 3)

# Filter out NA values first
valid_data <- !is.na(sncrna_neg_zscores$distance) & !is.na(sncrna_neg_zscores$GC_percentage)

# Sort the valid data
sorted_idx <- order(sncrna_neg_zscores$distance[valid_data])
sorted_distance <- sncrna_neg_zscores$distance[valid_data][sorted_idx]
sorted_gc <- sncrna_neg_zscores$GC_percentage[valid_data][sorted_idx]

# Fit LOESS on sorted, valid data
loess_fit <- loess(sorted_gc ~ sorted_distance, span = 0.1)
pred <- predict(loess_fit, se = TRUE)

# Add LOESS curve and confidence intervals
lines(sorted_distance, pred$fit, col = "red", lwd = 3)
lines(sorted_distance, pred$fit + 1.96 * pred$se.fit, col = "grey", lty = 3)
lines(sorted_distance, pred$fit - 1.96 * pred$se.fit, col = "grey", lty = 3)



heatscatter(lncrna_neg_zscores$distance, 
            lncrna_neg_zscores$GC_percentage,
            daltonize = TRUE, 
            cor = FALSE,
            log = "x",
            xlim = c(7500, 5100000),
            xlab = "Distance to Gene",
            ylab = "Robust Z-score",
            main = "",
            xaxt = 'n',
            cex.axis = 2,
            cex.lab = 2.5)
axis(1, at = c(1000, 10000, 100000, 1000000, 5000000), 
     labels = c("1kb", "10kb", "100kb", "1Mb", "5Mb"),
     cex.axis = 2)
title(main = "lncRNA", cex.main = 3)

# Filter out NA values first
valid_data <- !is.na(lncrna_neg_zscores$distance) & !is.na(lncrna_neg_zscores$GC_percentage)

# Sort the valid data
sorted_idx <- order(lncrna_neg_zscores$distance[valid_data])
sorted_distance <- lncrna_neg_zscores$distance[valid_data][sorted_idx]
sorted_gc <- lncrna_neg_zscores$GC_percentage[valid_data][sorted_idx]

# Fit LOESS on sorted, valid data
loess_fit <- loess(sorted_gc ~ sorted_distance, span = 0.1)
pred <- predict(loess_fit, se = TRUE)

# Add LOESS curve and confidence intervals
lines(sorted_distance, pred$fit, col = "red", lwd = 3)
lines(sorted_distance, pred$fit + 1.96 * pred$se.fit, col = "grey", lty = 3)
lines(sorted_distance, pred$fit - 1.96 * pred$se.fit, col = "grey", lty = 3)


# SNP density
heatscatter(prot_neg_zscores$distance, 
            prot_neg_zscores$SNP_density,
            daltonize = TRUE, 
            cor = FALSE,
            log = "x",
            xlim = c(7500, 5100000),
            ylim = c(-3, 3),
            xlab = "Distance to Gene",
            ylab = "Robust Z-score",
            main = "",
            xaxt = 'n',
            cex.axis = 2,
            cex.lab = 2.5)
axis(1, at = c(1000, 10000, 100000, 1000000, 5000000), 
     labels = c("1kb", "10kb", "100kb", "1Mb", "5Mb"),
     cex.axis = 2)
title(main = "mRNA", cex.main = 3)

# Filter out NA values first
valid_data <- !is.na(prot_neg_zscores$distance) & !is.na(prot_neg_zscores$SNP_density)

# Sort the valid data
sorted_idx <- order(prot_neg_zscores$distance[valid_data])
sorted_distance <- prot_neg_zscores$distance[valid_data][sorted_idx]
sorted_gc <- prot_neg_zscores$SNP_density[valid_data][sorted_idx]

# Fit LOESS on sorted, valid data
loess_fit <- loess(sorted_gc ~ sorted_distance, span = 0.1)
pred <- predict(loess_fit, se = TRUE)

# Add LOESS curve and confidence intervals
lines(sorted_distance, pred$fit, col = "red", lwd = 3)
lines(sorted_distance, pred$fit + 1.96 * pred$se.fit, col = "grey", lty = 3)
lines(sorted_distance, pred$fit - 1.96 * pred$se.fit, col = "grey", lty = 3)

heatscatter(sncrna_neg_zscores$distance, 
            sncrna_neg_zscores$SNP_density,
            daltonize = TRUE, 
            cor = FALSE,
            log = "x",
            xlim = c(7500, 5100000),
            ylim = c(-3, 3),
            xlab = "Distance to Gene",
            ylab = "Robust Z-score",
            main = "",
            xaxt = 'n',
            cex.axis = 2,
            cex.lab = 2.5)
axis(1, at = c(1000, 10000, 100000, 1000000, 5000000), 
     labels = c("1kb", "10kb", "100kb", "1Mb", "5Mb"),
     cex.axis = 2)
title(main = "sncRNA", cex.main = 3)

# Filter out NA values first
valid_data <- !is.na(sncrna_neg_zscores$distance) & !is.na(sncrna_neg_zscores$SNP_density)

# Sort the valid data
sorted_idx <- order(sncrna_neg_zscores$distance[valid_data])
sorted_distance <- sncrna_neg_zscores$distance[valid_data][sorted_idx]
sorted_gc <- sncrna_neg_zscores$SNP_density[valid_data][sorted_idx]

# Fit LOESS on sorted, valid data
loess_fit <- loess(sorted_gc ~ sorted_distance, span = 0.1)
pred <- predict(loess_fit, se = TRUE)

# Add LOESS curve and confidence intervals
lines(sorted_distance, pred$fit, col = "red", lwd = 3)
lines(sorted_distance, pred$fit + 1.96 * pred$se.fit, col = "grey", lty = 3)
lines(sorted_distance, pred$fit - 1.96 * pred$se.fit, col = "grey", lty = 3)

heatscatter(lncrna_neg_zscores$distance, 
            lncrna_neg_zscores$SNP_density,
            daltonize = TRUE, 
            cor = FALSE,
            log = "x",
            xlim = c(7500, 5100000),
            ylim = c(-3, 3),
            xlab = "Distance to Gene",
            ylab = "Robust Z-score",
            main = "",
            xaxt = 'n',
            cex.axis = 2,
            cex.lab = 2.5)
axis(1, at = c(1000, 10000, 100000, 1000000, 5000000), 
     labels = c("1kb", "10kb", "100kb", "1Mb", "5Mb"),
     cex.axis = 2)
title(main = "lncRNA", cex.main = 3)
# Filter out NA values first
valid_data <- !is.na(lncrna_neg_zscores$distance) & !is.na(lncrna_neg_zscores$SNP_density)

# Sort the valid data
sorted_idx <- order(lncrna_neg_zscores$distance[valid_data])
sorted_distance <- lncrna_neg_zscores$distance[valid_data][sorted_idx]
sorted_gc <- lncrna_neg_zscores$SNP_density[valid_data][sorted_idx]

# Fit LOESS on sorted, valid data
loess_fit <- loess(sorted_gc ~ sorted_distance, span = 0.1)
pred <- predict(loess_fit, se = TRUE)

# Add LOESS curve and confidence intervals
lines(sorted_distance, pred$fit, col = "red", lwd = 3)
lines(sorted_distance, pred$fit + 1.96 * pred$se.fit, col = "grey", lty = 3)
lines(sorted_distance, pred$fit - 1.96 * pred$se.fit, col = "grey", lty = 3)
dev.off()

