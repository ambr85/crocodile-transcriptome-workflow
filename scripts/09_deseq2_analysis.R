###########################################################
#Program:    Differential Expression DESeq2
#Author:       Marisol Buenfil-Rojas
#Date:       2025 (FINAL)
###########################################################

# Set working directory ← Replace with your actual path
directory <- "/your_directory/"
setwd(directory)
getwd()

# Load required libraries
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(zoo)
library(ggplot2)
library(dplyr)
library(ggrepel)

# ----------------------
# 1. Read HTSeq count files
# ----------------------
sampleFiles <- c("REFERENCE1count.out", "REFERENCE2count.out", "REFERENCE3count.out", 
                 "SITE1SAMPLE1count.out", "SITE1SAMPLE2count.out", "SITE1SAMPLE3count.out", 
                 "SITE2SAMPLE1count.out", "SITE2SAMPLE2count.out", "SITE2SAMPLE3count.out")

sampleCondition <- c("REF", "REF", "REF",
                     "S1S1", "S1S2", "S1S3", 
                     "S2S1", "S2S2", "S2S3")

sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design = ~ condition)

# Set factor levels (control first)
dds$condition <- factor(dds$condition, levels = c("REF", "S1", "S2"))

# ----------------------
# 2. Filter low-expressed genes
# ----------------------
dds <- dds[rowSums(counts(dds)) >= 10, ]

# ----------------------
# 3. Run DESeq2
# ----------------------
dds <- DESeq(dds)

# ----------------------
# 4. rlog normalization
# ----------------------
rld <- rlog(dds, blind = TRUE)
rld_mat <- assay(rld)

rld_mat_filtered <- rld_mat[!grepl("^LOC", rownames(rld_mat)), ]
sum(grepl("^LOC", rownames(rld_mat)))

# ----------------------
# 5. Extract DE results
# ----------------------
resC1 <- results(dds, contrast = c("condition", "REF", "S1"))
resC2 <- results(dds, contrast = c("condition", "REF", "S2"))

# Save results
write.csv(as.data.frame(resC1), file = "DESeq2_Results_REF_vs_S1.csv")
write.csv(as.data.frame(resC2), file = "DESeq2_Results_REF_vs_S2.csv")

# For mixOmics:
# Convert to data.frame and add the names of genes as first column
expr_df <- as.data.frame(expr_matrix)
expr_df$Gene <- rownames(expr_df)

# Reorganize to set GENE as first column
expr_df <- expr_df[, c("Gene", setdiff(colnames(expr_df), "Gene"))]
write.csv(expr_df, file = "NormalizedCounts_for_mixOmics.csv")

# ----------------------
# 6. Select significant genes and plot heatmap
# ----------------------
#Explore Concentrations individually
# Threshold
alpha <- 0.05
lfc_threshold <- 1  # You can remove this line if you don’t want to filter by fold change

# --- REF vs S1 ---
resC1 <- results(dds, contrast = c("condition", "REF", "S1"))
resC1 <- resC1[order(resC1$padj), ]  
summary(resC1)

# Filter by padj or pvalue, and optional log2FC
sigC1_genes <- rownames(resC1)[
  !is.na(resC1$pvalue) & 
    resC1$pvalue <= alpha &
    abs(resC1$log2FoldChange) >= lfc_threshold
  # Replace with: resC1$pvalue <= alpha  ← if padj too strict
]

# Save results
write.csv(as.data.frame(resC1), file = "DESeq2_Results_REF_vs_S1.csv") #Site1_vs_REF
write.csv(as.data.frame(resC1[sigC1_genes, ]), file = "Significant_Genes_REFvsS1.csv") #Significant genes REFvsSite1

# --- REF vs S2 ---
resC2 <- results(dds, contrast = c("condition", "REF", "S2"))
resC2 <- resC2[order(resC2$padj), ]
summary(resC2)

# Filter by padj or pvalue, and optional log2FC
sigC2_genes <- rownames(resC2)[
  !is.na(resC2$pvalue) &
    resC2$pvalue <= alpha &
    abs(resC2$log2FoldChange) >= lfc_threshold
  # Replace with: resC2$pvalue <= alpha  ← for exploratory plots if padj too strict
]

write.csv(as.data.frame(resC2), file = "DESeq2_Results_REF_vs_S2.csv")
write.csv(as.data.frame(resC2[sigC2_genes, ]), file = "Significant_Genes_REFvsS2.csv")

# Optional: heatmap per set (editar)
col.pan <- colorpanel(50, "blue", "white", "red")

# Create heatmap input from rlog data
rld <- rlog(dds, blind = TRUE)
rld_mat <- assay(rld)

# Select only the samples REF & Sites (cada uno)
selected_samples_C1 <- colData(dds)$condition %in% c("REF", "S1")
selected_samples_C2 <- colData(dds)$condition %in% c("REF", "S2")

# Subset of the rlog matrix with the genes and the samples you want to check
heat_data_C1 <- rld_mat[sigC1_genes, selected_samples_C1]
heat_data_C2 <- rld_mat[sigC2_genes, selected_samples_C2]

# Etiquetas de condición para esas muestras
sampleLabels_C1 <- colData(dds)$condition[selected_samples_C1]
sampleLabels_C2 <- colData(dds)$condition[selected_samples_C2]

# Subset and plot if there are enough genes
if (length(sigC1_genes) >= 2) {
  heatmap.2(as.matrix(heat_data_C1),
            trace = "none", density = "none", scale = "row",
            col = col.pan, cexRow = 0.7, cexCol = 0.9,
            labCol = sampleLabels_C1,
            main = "REF vs S1: DEGs")
}

if (length(sigC2_genes) >= 2) {
  heatmap.2(as.matrix(heat_data_C2),
            trace = "none", density = "none", scale = "row",
            col = col.pan, cexRow = 0.7, cexCol = 0.9,
            labCol = sampleLabels_C2,
            main = "REF vs S2: DEGs")
}

#Explore the whole dataset (you can change padj to pvalue if too strict) (THIS IS IN CASE YOU DON'T DO THE SEPARATE HEATMAPS)
alpha <- 0.05
sigC1 <- rownames(resC1)[resC1$padj <= alpha & !is.na(resC1$padj)]
sigC2 <- rownames(resC2)[resC2$padj <= alpha & !is.na(resC2$padj)]

# Union or intersection accordingly
union_genes <- union(sigC1_genes, sigC2_genes)  # Genes that change at least in one condition
#intersect_genes <- intersect(sigC1, sigC2)     # Genes that change in both conditions
# To check using the genes from the individual heatmaps, replace for sigC1_genes and sigC2_genes, but if you do it 
# directly, use sigC1 and sigC2.

# Subset para heatmap
heat_data <- rld_mat[union_genes, ]
#heat_data <- rld_mat[intersect_genes, ]

# Handle -Inf and NA values
heat_data[heat_data == -Inf] <- NA
heat_data <- na.aggregate(heat_data, FUN = mean)

# Heatmap
col.pan <- colorpanel(50, "blue", "white", "red") #Blue (downregulated), White (no change), Red (upregulated)
sampleLabels <- colData(dds)$condition

heatmap.2(as.matrix(heat_data),
          trace = "none",
          density = "none",
          scale = "row",
          col = col.pan,
          cexRow = 0.7,
          cexCol = 0.9,
          labCol = sampleLabels,
          main = "Heatmap of DEGs in Sites")

write.csv(as.data.frame(rld_mat[union_genes, ]), file = "Rlog_Transformed_Union_DEGs.csv")
# write.csv(as.data.frame(rld_mat[intersect_genes, ]), file = "Rlog_Transformed_Intersect_DEGs.csv")
write.csv(as.data.frame(rld_mat), file = "Rlog_AllGenes.csv")

# ----------------------
# 7. Original MA plot (All dataset)
# ----------------------
plotMA(dds, ylim = c(-5, 5), main = "MAplot")
dev.copy(png, "deseq2_MAplotREFvsSites.png")
dev.off()

# ----------------------
# 8. rlog-transformed MA plot (manual)
# ----------------------
# Adjust condition names as needed — using REF vs Site1 here
conditionA <- colData(rld)$condition == "S1" #Site here
conditionB <- colData(rld)$condition == "REF" #REF

log2FC <- rowMeans(rld_mat[, conditionA]) - rowMeans(rld_mat[, conditionB])
baseMean <- rowMeans(rld_mat)

df1 <- data.frame(BaseMean = baseMean, log2FC = log2FC)

ggplot(df1, aes(x = BaseMean, y = log2FC)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red") +
  scale_x_log10() +
  labs(title = "MA Plot - rlog Normalized REF vs S1", x = "Mean Expression", y = "Log2 Fold Change") +
  theme_minimal()

# Adjust condition names as needed — using REF vs Site2 here
conditionA <- colData(rld)$condition == "S2"
conditionB <- colData(rld)$condition == "REF"

log2FC <- rowMeans(rld_mat[, conditionA]) - rowMeans(rld_mat[, conditionB])
baseMean <- rowMeans(rld_mat)

df2 <- data.frame(BaseMean = baseMean, log2FC = log2FC)

ggplot(df2, aes(x = BaseMean, y = log2FC)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red") +
  scale_x_log10() +
  labs(title = "MA Plot - rlog Normalized REF vs S2", x = "Mean Expression", y = "Log2 Fold Change") +
  theme_minimal()

# ----------------------
# 9. PCA and labeling
# ----------------------
z <- plotPCA(rld, intgroup = c("condition"))
z + geom_text(size = 2, aes(label = name), nudge_y = -0.5)

# ----------------------
# 10. Volcano plots
# ----------------------
# Significance criteria
alpha <- 0.05
lfc_cutoff <- 1

# === Volcano plot for resC1 (padj)===
resC1 <- resC1[!is.na(resC1$pvalue) & !is.na(resC1$padj), ]

colsC1 <- densCols(resC1$log2FoldChange, -log10(resC1$pvalue))

plot(resC1$log2FoldChange, -log10(resC1$padj),
     col = colsC1,
     panel.first = grid(),
     main = "Volcano Plot - REF vs S1", #SiteX vs REF
     xlab = "log2(Fold Change)",
     ylab = "-log10(p-adj)",
     pch = 20,
     cex = 0.6)

abline(v = c(-lfc_threshold, lfc_threshold), col = "brown", lty = 2)
abline(h = -log10(alpha), col = "brown", lty = 2)

gnC1 <- abs(resC1$log2FoldChange) > lfc_threshold & resC1$padj < alpha

text(resC1$log2FoldChange[gnC1],
     -log10(resC1$padj)[gnC1],
     labels = rownames(resC1)[gnC1],
     cex = 0.6, pos = 3, offset = 0.2, col = "black")

# === Volcano plot for resC2 (padj)===
resC2 <- resC2[!is.na(resC2$pvalue) & !is.na(resC2$padj), ]

colsC2 <- densCols(resC2$log2FoldChange, -log10(resC2$pvalue))

plot(resC2$log2FoldChange, -log10(resC2$padj),
     col = colsC2,
     panel.first = grid(),
     main = "Volcano Plot - REF vs S2",
     xlab = "log2(Fold Change)",
     ylab = "-log10(p-adj)",
     pch = 20,
     cex = 0.6)

abline(v = c(-lfc_threshold, lfc_threshold), col = "brown", lty = 2)
abline(h = -log10(alpha), col = "brown", lty = 2)

gnC2 <- abs(resC2$log2FoldChange) > lfc_threshold & resC2$padj < alpha

text(resC2$log2FoldChange[gnC2],
     -log10(resC2$padj)[gnC2],
     labels = rownames(resC2)[gnC2],
     cex = 0.6, pos = 3, offset = 0.2, col = "black")

###Using p-value instead of padj (if data are too strictly curated by padj filter)
# Parameters
alpha <- 0.05  # p-value threshold
lfc_threshold <- 1  # log2 Fold Change threshold

# Volcano plot for Site1

# Add labels
resC1_df <- as.data.frame(resC1) %>%
  mutate(Gene = rownames(resC1),
         Significant = ifelse(!is.na(pvalue) & pvalue < alpha & abs(log2FoldChange) > lfc_threshold, "Yes", "No"))

# Save significant genes
sigC1_subset <- resC1_df %>% filter(Significant == "Yes")
write.csv(sigC1_subset, file = "Volcano_SignificantGenes_C1.csv", row.names = FALSE)

# Plot
ggplot(resC1_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Significant)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(data = sigC1_subset, aes(label = Gene),
                  size = 3, max.overlaps = 15, box.padding = 0.4) +
  scale_color_manual(values = c("No" = "gray", "Yes" = "firebrick")) +
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), color = "darkblue", linetype = "dashed") +
  geom_hline(yintercept = -log10(alpha), color = "darkblue", linetype = "dashed") +
  labs(title = "Volcano Plot: REF vs S1", #Site1 vs REF
       x = "log2(Fold Change)",
       y = "-log10(p-value)") +
  theme_minimal()

# Volcano plot for Site2 (EDIT)

resC2_df <- as.data.frame(resC2) %>%
  mutate(Gene = rownames(resC2),
         Significant = ifelse(!is.na(pvalue) & pvalue < alpha & abs(log2FoldChange) > lfc_threshold, "Yes", "No"))

sigC2_subset <- resC2_df %>% filter(Significant == "Yes")
write.csv(sigC2_subset, file = "Volcano_SignificantGenes_C2.csv", row.names = FALSE)

ggplot(resC2_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Significant)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(data = sigC1_subset, aes(label = Gene),
                  size = 3, max.overlaps = 15, box.padding = 0.4) +
  scale_color_manual(values = c("No" = "gray", "Yes" = "firebrick")) +
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), color = "darkblue", linetype = "dashed") +
  geom_hline(yintercept = -log10(alpha), color = "darkblue", linetype = "dashed") +
  labs(title = "Volcano Plot: REF vs S2",
       x = "log2(Fold Change)",
       y = "-log10(p-value)") +
  theme_minimal()

###############################
## In case FDR is too stringent, p-value can be used. In the article state: "Due to the small number of DEGs passing FDR (False Discovery Rate, padj value) thresholds, 
##  we explored genes with p-value < 0.1 and |log2FC| > 1 as potential candidates."
