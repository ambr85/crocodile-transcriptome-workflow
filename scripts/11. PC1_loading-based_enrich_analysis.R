################################################################################
# PC1 loading-based enrichment analysis
# Purpose: Interpret the biological meaning of PC1
#          Clarify "What does PC1 represent? Structure? Stress? Metabolism?"
# This script is for: 
#   1. Extract PC1 loadings from PCA results
################################################################################

# Required libraries
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)  # or appropriate organism database
library(ggplot2)
library(dplyr)
library(tidyr)
library(enrichplot)

#Load data
directory <- "/YOUR_DIRECTORY"
setwd(directory)
getwd()

# ---- 1. Extract PC1 loadings ----
# Get normalized data from DESeq2 object
# (Assuming rlog/vst transformation is already done)
dds <- readRDS("dds_full.rds") #Your DESeq2 normalized data and transformed
rld <- rlog(dds, blind=FALSE)  # or vst(dds)

# Perform PCA (skip if already performed)
pca_data <- plotPCA(rld, intgroup="condition", returnData=TRUE)
z <- plotPCA(rld, intgroup="condition")
z + geom_text(size=2, aes(label=name), nudge_y = -0.5)
pca_result <- prcomp(t(assay(rld)))

# Extract PC1 loadings
pc1_loadings <- pca_result$rotation[, 1]
names(pc1_loadings) <- rownames(pca_result$rotation)

# Sort by loading (descending order)
pc1_loadings_sorted <- sort(pc1_loadings, decreasing = TRUE)

# Convert to data frame
pc1_loadings_df <- data.frame(
  gene = names(pc1_loadings_sorted),
  loading = pc1_loadings_sorted,
  rank = 1:length(pc1_loadings_sorted)
)

# ---- 2. Select top/bottom loading genes ----
# Top positive loading genes (e.g., top 300 genes)
n_top <- 300
top_positive_genes <- head(pc1_loadings_df$gene, n_top)

# Top negative loading genes (e.g., bottom 300 genes)
top_negative_genes <- tail(pc1_loadings_df$gene, n_top)

# dir.create("./Step2/results", showWarnings = FALSE, recursive = TRUE)
# file.exists("./Step2/results")  # TRUE if exists

save(pc1_loadings_df,  file = "./Step2/results/pc1_loadings_df.dat")
save(pc1_loadings_sorted, file = "./Step2/results/pc1_loadings_sorted.dat")
save(top_positive_genes,  file = "./Step2/results/top_positive_genes.dat")
save(top_negative_genes,  file = "./Step2/results/top_negative_genes.dat")

write.csv(pc1_loadings_df, 
          "./Step2/results/PC1_loadings_all_genes.csv", 
          row.names = FALSE)

################################################################################
# Fig.2A: PC1 loading-based enrichment analysis
#This script is for:
#  1.5 gene ID conversion for  
################################################################################

# Required libraries
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)  # or appropriate organism database
library(ggplot2)
library(dplyr)
library(tidyr)
library(enrichplot)

# ---- 3. Prepare gene ID conversion ----
# If Alligator annotation is not available, conversion to human ortholog may be needed
# Example: conversion from ENSEMBL ID to ENTREZ ID

# load("./Step2/pc1_loadings_df.dat")

# If gene IDs are in ENSEMBL format
gene_ids <- pc1_loadings_df$gene
sym <- names(pc1_loadings_sorted)

# For simplicity, assume ENTREZ IDs or symbols are already available
# universe (background gene set) = all 18,106 detected genes
universe_genes <- pc1_loadings_df$gene

# Convert to ENTREZID (adjust according to organism database)
# library(biomaRt)
# mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# gene_mapping <- getBM(
#   attributes=c('ensembl_gene_id', 'entrezgene_id', 'external_gene_name'),
#   filters='ensembl_gene_id',
#   values=gene_ids_clean,
#   mart=mart
# )

# Map SYMBOL -> ENTREZ (only once)
sym2entrez <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys     = gene_ids,  # universe_genes in symbol
  keytype  = "SYMBOL",
  column   = "ENTREZID",
  multiVals = "first"
)

sym2entrez2 <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys     = sym,
  keytype  = "SYMBOL",
  column   = "ENTREZID",
  multiVals = "first"
)

# Replace the CONTENT (does not change the name of the objects)
universe_genes   <- unname(sym2entrez[universe_genes])
top_positive_genes<- unname(sym2entrez[top_positive_genes])
top_negative_genes<- unname(sym2entrez[top_negative_genes])

## Remove those that were not mapped
keep <- !is.na(sym2entrez2) & nzchar(sym2entrez2)
vals <- pc1_loadings_sorted[keep]
names(vals) <- sym2entrez2[keep]    # ← replace SYMBOL for ENTREZ (only ENTREZ)

# Filter NAs (recommended for ORA/GSEA)
universe_genes    <- unique(universe_genes[!is.na(universe_genes)])
top_positive_genes <- unique(top_positive_genes[!is.na(top_positive_genes)])
top_negative_genes <- unique(top_negative_genes[!is.na(top_negative_genes)])

## Collapse duplicates (if one ENTREZ appears multiple times: take the major |loading|)
tmp <- split(as.numeric(vals), names(vals))
pc1_entrez <- sapply(tmp, function(x) x[which.max(abs(x))])
pc1_entrez <- sort(pc1_entrez, decreasing = TRUE)

save(universe_genes, file = "./Step2/universe_genes.dat")
save(top_positive_genes, file = "./Step2/top_positive_genes.dat")
save(top_negative_genes, file = "./Step2/top_negative_genes.dat")

################################################################################
# Fig.2A: PC1 loading-based enrichment analysis
#This script is for:
#  2. Perform ORA (Over-Representation Analysis) using top/bottom loading genes
################################################################################

# Required libraries
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)  # or appropriate organism database
library(ggplot2)
library(dplyr)
library(tidyr)
library(enrichplot)

# ---- 0. load data -----------------------------
load("./Step2/top_positive_genes.dat")
load("./Step2/top_negative_genes.dat")
load("./Step2/universe_genes.dat")

#universe_genes <- pc1_loadings_df$gene

# ---- 4. ORA (Over-Representation Analysis) ----
# 4-1. GO enrichment for positive loading genes
ego_positive <- enrichGO(
  gene          = top_positive_genes,
  universe      = universe_genes,
  OrgDb         = org.Hs.eg.db,  # Change to appropriate organism DB
  keyType       = "ENTREZID",      # or "ENTREZID", "SYMBOL"
  ont           = "BP",            # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# 4-2. GO enrichment for negative loading genes
ego_negative <- enrichGO(
  gene          = top_negative_genes,
  universe      = universe_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# 4-3. KEGG enrichment (positive loading)
ekegg_positive <- enrichKEGG(
  gene         = top_positive_genes,
  universe     = universe_genes,
  organism     = "hsa",  # human; change to appropriate code for alligator
  keyType      = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# 4-4. KEGG enrichment (negative loading)
ekegg_negative <- enrichKEGG(
  gene         = top_negative_genes,
  universe     = universe_genes,
  organism     = "hsa",
  keyType      = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# ---- save data -----------------------------
save(ego_positive, file = "./Step2/results/ego_positive.dat")
save(ego_negative, file = "./Step2/results/ego_negative.dat")
save(ekegg_positive, file = "./Step2/results/ekegg_positive.dat")
save(ekegg_negative, file = "./Step2/results/ekegg_negative.dat")

# ORA results
write.csv(as.data.frame(ego_positive), 
          "./Step2/results/PC1_ORA_GO_positive.csv", 
          row.names = FALSE)
write.csv(as.data.frame(ego_negative), 
          "./Step2/results//PC1_ORA_GO_negative.csv", 
          row.names = FALSE)
write.csv(as.data.frame(ekegg_positive), 
          "./Step2/results//PC1_ORA_KEGG_positive.csv", 
          row.names = FALSE)
write.csv(as.data.frame(ekegg_negative), 
          "./Step2/results//PC1_ORA_KEGG_negative.csv", 
          row.names = FALSE)

################################################################################
# Fig.2A: PC1 loading-based enrichment analysis
#This script is for:
#  3. Rank all genes by PC1 loading and perform GSEA
################################################################################

# Required libraries
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)  # or appropriate organism database
library(ggplot2)
library(dplyr)
library(tidyr)
library(enrichplot)

# ---- 0. load data -----------------------------
load("./Step2/pc1_loadings_sorted.dat")

# ---- 5. GSEA (Gene Set Enrichment Analysis) ----
# Rank all genes by PC1 loading values
gene_list <- pc1_loadings_sorted
gene_list2 <- pc1_entrez

names(gene_list) <- names(pc1_loadings_sorted)  # Already named
names(gene_list2) <- names(pc1_entrez)  # Already named

# 5-1. GO GSEA
gsea_go <- gseGO(
  geneList      = gene_list,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  eps           = 0
)

# 5-2. KEGG GSEA
gsea_kegg <- gseKEGG(
  geneList      = gene_list2,
  organism      = "hsa",
  keyType       = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  eps           = 0
)

# ---- save data -----------------------------
save(gsea_go, file = "./Step2/results/gsea_go.dat")
save(gsea_kegg, file = "./Step2/results/gsea_kegg.dat")

# GSEA results
write.csv(as.data.frame(gsea_go), 
          "./Step2/results/PC1_GSEA_GO.csv", 
          row.names = FALSE)
write.csv(as.data.frame(gsea_kegg), 
          "./Step2/results/PC1_GSEA_KEGG.csv", 
          row.names = FALSE)

################################################################################
# Fig.2A: Visualization of PC1 loading-based enrichment results
# Purpose: Visualize enrichment results from PC1 analysis
# Methods: 
#   - Barplot/Dotplot for ORA results
#   - Enrichment map
#   - GSEA enrichment plot for top pathways
################################################################################

# Required libraries (if not already loaded)
library(ggplot2)
library(enrichplot)
library(cowplot)
library(RColorBrewer)

# ----0. load data ---------------------
load("./Step2/results/ego_positive.dat")
load("./Step2/results/ego_negative.dat")

# ---- 6. Barplot of ORA results (positive and negative) ----
# GO enrichment for positive loading genes (top 10 terms)
p1 <- barplot(ego_positive, 
              showCategory = 10,
              x = "GeneRatio", 
              title = "GO BP enriched in PC1+ loading genes\n(Farm-associated)") +
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold"))

# GO enrichment for negative loading genes (top 10 terms)
p2 <- barplot(ego_negative, 
              showCategory = 10,
              x = "GeneRatio", 
              title = "GO BP enriched in PC1- loading genes\n(Wild-associated)") +
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold"))

# ---- 7.plot Barplot of ORA results (positive and negative) ----
# Option 1: Arrange ORA barplots
fig2a_v1 <- plot_grid(p1, p2, 
                      ncol = 1, 
                      labels = c("A", "B"),
                      label_size = 14)

ggsave("./Step2/figures/Fig2A_PC1_enrichment_v1.pdf", 
       fig2a_v1, 
       width = 10, height = 12, 
       dpi = 300)
                     
################################################################################
# Fig.2A: Visualization of PC1 loading-based enrichment results
# Purpose: Visualize enrichment results from PC1 analysis
# Methods: 
#   - Barplot/Dotplot for ORA results
#   - Enrichment map
#   - GSEA enrichment plot for top pathways
################################################################################

# Required libraries (if not already loaded)
library(ggplot2)
library(enrichplot)
library(cowplot)
library(RColorBrewer)

# ----0. load data ---------------------
load("./Step2/results/ego_positive.dat")
load("./Step2/results/ego_negative.dat")
load("./Step2/results/gsea_go.dat")

# ---- 8. Dotplot (more informative visualization) ----
p3 <- dotplot(ego_positive, 
              showCategory = 15,
              x = "GeneRatio", 
              title = "PC1+ GO enrichment") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 9))

p4 <- dotplot(ego_negative, 
              showCategory = 15,
              x = "GeneRatio", 
              title = "PC1- GO enrichment") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 9))

# 8-1. GSEA dotplot (top and bottom terms)
p7 <- dotplot(gsea_go, 
              showCategory = 20,
              x = "GeneRatio", 
              split = ".sign",
              title = "GO GSEA based on PC1 loading") +
  facet_grid(. ~ .sign) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 8))

# ---- 9.plot Dotplot (more informative visualization) ----
# Option 2: Dotplot + GSEA
fig2a_v2 <- plot_grid(p3, p4, p7,
                      ncol = 1,
                      labels = c("A", "B", "C"),
                      label_size = 9,
                      rel_heights = c(1, 1, 1.2))

ggsave("./Step2/figures/Fig2A_PC1_enrichment_v2.pdf", 
       fig2a_v2, 
       width = 12, height = 24, 
       dpi = 300)

################################################################################
# Fig.2A: Visualization of PC1 loading-based enrichment results
# Purpose: Visualize enrichment results from PC1 analysis
# Methods: 
#   - Barplot/Dotplot for ORA results
#   - Enrichment map
#   - GSEA enrichment plot for top pathways
################################################################################

# Required libraries (if not already loaded)
library(ggplot2)
library(enrichplot)
library(cowplot)
library(RColorBrewer)

# ----0. load data ---------------------
load("./Step2/results/ekegg_positive.dat")
load("./Step2/results/ekegg_negative.dat")
load("./Step2/results/gsea_go.dat")

# ---- 10. KEGG enrichment barplot ----
if(nrow(as.data.frame(ekegg_positive)) > 0) {
  p5 <- barplot(ekegg_positive, 
                x = "GeneRatio", 
                showCategory = 10,
                title = "KEGG pathways (PC1+)") +
    theme_classic()
} else {
  p5 <- ggplot() + 
    annotate("text", x = 1, y = 1, 
             label = "No significant KEGG pathways") +
    theme_void()
}

if(nrow(as.data.frame(ekegg_negative)) > 0) {
  p6 <- barplot(ekegg_negative, 
                x = "GeneRatio", 
                showCategory = 10,
                title = "KEGG pathways (PC1-)") +
    theme_classic()
} else {
  p6 <- ggplot() + 
    annotate("text", x = 1, y = 1, 
             label = "No significant KEGG pathways") +
    theme_void()
}

# Option 3: 4-panel (ORA + KEGG)
fig2a_v3 <- plot_grid(p5, p6,
                      ncol = 2,
                      labels = c("A", "B"),
                      label_size = 14)

ggsave("./Step2/figures/Fig2A_PC1_enrichment_v3.pdf", 
       fig2a_v3, 
       width = 14, height = 10, 
       dpi = 300)

################################################################################
# Fig.2A: Visualization of PC1 loading-based enrichment results
# Purpose: Visualize enrichment results from PC1 analysis
# Methods: 
#   - Barplot/Dotplot for ORA results
#   - Enrichment map
#   - GSEA enrichment plot for top pathways
################################################################################

# Required libraries (if not already loaded)
library(ggplot2)
library(enrichplot)
library(cowplot)
library(RColorBrewer)

# ----0. load data ---------------------
load("./Step2/results/gsea_go.dat")

# 11. Enrichment plot for specific pathways
# Example: ECM-related terms like "extracellular matrix organization"
# Check term IDs
gsea_go_df <- as.data.frame(gsea_go)
head(gsea_go_df$Description, 20)

# Select interesting terms (e.g., most significant positive/negative enrichment)
top_term_pos <- gsea_go_df %>% 
  filter(NES > 0) %>% 
  slice_min(order_by = pvalue, n = 1) %>% 
  pull(ID)

top_term_neg <- gsea_go_df %>% 
  filter(NES < 0) %>% 
  slice_min(order_by = pvalue, n = 1) %>% 
  pull(ID)

# Enrichment plot
if(length(top_term_pos) > 0) {
  p8 <- gseaplot2(gsea_go, 
                  geneSetID = top_term_pos,
                  title = gsea_go_df[gsea_go_df$ID == top_term_pos, "Description"],
                  pvalue_table = TRUE)
}

if(length(top_term_neg) > 0) {
  p9 <- gseaplot2(gsea_go, 
                  geneSetID = top_term_neg,
                  title = gsea_go_df[gsea_go_df$ID == top_term_neg, "Description"],
                  pvalue_table = TRUE)
}

# Save individual GSEA plots
if(exists("p8")) {
  ggsave("./Step2/figures/Fig2A_GSEA_top_positive.pdf", 
         p8, 
         width = 8, height = 6, 
         dpi = 300)
}

if(exists("p9")) {
  ggsave("./Step2/figures/Fig2A_GSEA_top_negative.pdf", 
         p9, 
         width = 9, height = 6, 
         dpi = 300)
}                     
