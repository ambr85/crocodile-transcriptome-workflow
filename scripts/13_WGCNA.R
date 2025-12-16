################################################################################
# Fig.3: WGCNA - Weighted Gene Co-expression Network Analysis
# Purpose: Separate "structural" and "stress" genes into distinct modules
#          Confirm that structural and stress response genes operate as 
#          separate modules in the co-expression network
#This script is for:
#   1. Select highly variable genes from rlog/vst expression data
################################################################################

# Required libraries
library(devtools)
library(WGCNA)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)  # or appropriate organism database
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Enable multi-threading for WGCNA (adjust based on your system)
enableWGCNAThreads(nThreads = 4)

# ---- 0. load data ----
#in_f <- "set path for your expression data"
#rld <- get(load(in_f))

#dds <- readRDS("dds_full.rds") #already loaded for Step 2 and in main directory

# ---- 1. Prepare expression data for WGCNA ----
# Get normalized expression data (rlog or vst transformed)
# Assuming 'rld' is your rlog/vst transformed DESeq2 object
#rld <- rlog(dds, blind=FALSE) #Not necessary if you have loaded everything in the same project and directory

expr_mat <- assay(rld)  # genes x samples

# Transpose for WGCNA (samples x genes)
datExpr <- t(expr_mat)

# ---- 2. Select highly variable genes ----
# Calculate variance for each gene
gene_vars <- apply(expr_mat, 1, var)

# Select top variable genes (e.g., top 5000-10000)
n_genes <- 8000
top_var_genes <- names(sort(gene_vars, decreasing = TRUE)[1:n_genes])

# Subset expression matrix
datExpr_filtered <- datExpr[, top_var_genes]

# Check for genes and samples with too many missing values
gsg <- goodSamplesGenes(datExpr_filtered, verbose = 3)
if (!gsg$allOK) {
  # Remove problematic genes/samples if necessary
  if (sum(!gsg$goodGenes) > 0)
    datExpr_filtered <- datExpr_filtered[, gsg$goodGenes]
  if (sum(!gsg$goodSamples) > 0)
    datExpr_filtered <- datExpr_filtered[gsg$goodSamples, ]
}

# ----- save data ---------------
save(rld, file = "./Step4/results/rld.dat")
save(datExpr_filtered, file = "./Step4/results/datExpr_filtered.dat")

# ---- 3. Sample clustering to detect outliers ----
sampleTree <- hclust(dist(datExpr_filtered), method = "average")

# Plot sample dendrogram
pdf("./Step4/figures/WGCNA_sample_clustering.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", 
     sub = "", xlab = "", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

# Optional: Remove outlier samples if needed
# Visually inspect and set cutoff
cutHeight <- 95  # example
clust <- cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 2)
keepSamples <- (clust == 1)
datExpr_filtered <- datExpr_filtered[keepSamples, ]

# ---- 4. Choose soft-thresholding power (beta) ----
# Test different powers to find scale-free topology
powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))

sft <- pickSoftThreshold(datExpr_filtered, 
                         powerVector = powers,
                         networkType = "signed",
                         verbose = 5)

# Save soft threshold selection plot
pdf("./Step4/figures/WGCNA_soft_threshold_selection.pdf", width = 12, height = 5)
par(mfrow = c(1, 2))
cex1 <- 0.9

# Scale-free topology fit index
plot(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = paste("Scale independence"))
text(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.80, col = "red")  # Threshold for R^2

# Mean connectivity
plot(sft$fitIndices[, 1], 
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity",
     type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], 
     sft$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")
dev.off()

# Select power (choose lowest power with R^2 > 0.80)
softPower <- sft$powerEstimate
if (is.na(softPower)) {
  softPower <- 6  # default if automatic selection fails
  warning("Automatic power selection failed. Using default power = 6")
}
print(paste("Selected soft threshold power:", softPower))

################################################################################
# Fig.3: WGCNA - Weighted Gene Co-expression Network Analysis
# Purpose: Separate "structural" and "stress" genes into distinct modules
#          Confirm that structural and stress response genes operate as 
#          separate modules in the co-expression network
#This script is for:
#   2. Estimate co-expression modules using WGCNA
################################################################################

# Required libraries
library(WGCNA)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)  # or appropriate organism database
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# --- 0. load data ------------------------------
load("./Step4/results/datExpr_filtered.dat")
softPower = 16 #select softPower

# ---- 5. Construct network and detect modules (one-step approach) ----
# This uses automatic module detection with dynamic tree cut
net <- blockwiseModules(
  datExpr_filtered,
  power = softPower,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,           # Minimum module size
  reassignThreshold = 0,
  mergeCutHeight = 0.25,        # Merge similar modules
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "./Step4/results/WGCNA_TOM",
  verbose = 3,
  maxBlockSize = 10000          # Adjust based on RAM
)

# Convert numeric labels to colors
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)

# Number of modules
nModules <- length(unique(moduleColors)) - 1  # Exclude grey (unassigned)
print(paste("Number of modules detected:", nModules))

# Module sizes
table(moduleColors)

lab_num   <- as.integer(net$colors)                 # 0,1,2,...
col_name  <- WGCNA::labels2colors(lab_num)          # grey, blue, ...
map_df    <- unique(data.frame(label=lab_num, color=col_name))
map_df    <- map_df[order(map_df$label), ]

lut_label2color <- setNames(map_df$color, as.character(map_df$label))  # "1"->"brown"
lut_color2label <- setNames(as.character(map_df$label), map_df$color)  # "brown"->"1"

# ---- save -------------------------------------
save(net, file = "./Step4/results/net.dat")
save(moduleColors, file = "./Step4/results/moduleColors.dat")

################################################################################
# Fig.3: WGCNA - Weighted Gene Co-expression Network Analysis
# Purpose: Separate "structural" and "stress" genes into distinct modules
#          Confirm that structural and stress response genes operate as 
#          separate modules in the co-expression network
#This script is for:
#    3. Calculate module eigengenes for each sample
################################################################################

# --- 0. load data ------------------------------
load("./Step4/results/rld.dat")
load("./Step4/results/datExpr_filtered.dat")
load("./Step4/results/net.dat")
load("./Step4/results/moduleColors.dat")

# Create binary condition (farm/wild) and site (before in condition in the original metadata)

## 1) Define 'site' (your real labels) from your original column
colData(rld)$site <- factor(as.character(colData(rld)$condition),
                            levels = c("CC","SL","CH","CBR","ML","IL","NH"))

## 2) Define binary 'condition' Farm/Wild from 'site'
farmed_sites <- c("CC")
colData(rld)$condition <- factor(
  ifelse(colData(rld)$site %in% farmed_sites, "Farm", "Wild"),
  levels = c("Farm","Wild")
)

# ---- 6. Calculate module eigengenes (MEs) ----
MEs <- net$MEs
MEs <- orderMEs(MEs)  # Order by hierarchical clustering

# Rename for clarity
colnames(MEs) <- gsub("^ME", "", colnames(MEs))

modules_num_all  <- colnames(MEs)                  # includes "0" if it exists
modules_num_keep <- setdiff(modules_num_all, "0")  # leave out grey

me_colors <- lut_label2color[colnames(MEs)]  # vector relating names-numbers: "1"->"brown", "0"->"grey"

# Add sample information
sample_info <- data.frame(
  sample = rownames(datExpr_filtered),
  condition = colData(rld)[rownames(datExpr_filtered), "condition"],  # Farm/Wild
  site = colData(rld)[rownames(datExpr_filtered), "site"]            # Site info
)

MEs_with_info <- cbind(sample_info, MEs)

# ---- 7. Correlate module eigengenes with traits ----
# Create trait data (Farm = 1, Wild = 0)
traitData <- data.frame(
  sample = sample_info$sample,
  Farm_vs_Wild = as.numeric(sample_info$condition == "Farm")
)

# Add PC1 scores if available
if (exists("pca_result")) {
  pc1_scores <- pca_result$x[rownames(datExpr_filtered), 1]
  traitData$PC1 <- pc1_scores
}

# Match order
traitData <- traitData[match(rownames(MEs), traitData$sample), ]

# Calculate correlations
moduleTraitCor <- cor(MEs, traitData[, -1], use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr_filtered))

mod_rows <- modules_num_keep  # or modules_num_all if you want to see grey
moduleTraitCor    <- moduleTraitCor[mod_rows, , drop=FALSE]
moduleTraitPvalue <- moduleTraitPvalue[mod_rows, , drop=FALSE]

yLabs <- paste0("ME", mod_rows, " (", lut_label2color[mod_rows], ")")

################################################################################
# Fig.3: WGCNA - Weighted Gene Co-expression Network Analysis
# Purpose: Separate "structural" and "stress" genes into distinct modules
#          Confirm that structural and stress response genes operate as 
#          separate modules in the co-expression network
#This script is for:
#   4. Annotate each module with GO/KEGG enrichment
################################################################################

# --- 0. load data ------------------------------
load("./Step4/results/rld.dat")
load("./Step4/results/datExpr_filtered.dat")
load("./Step4/results/net.dat")
load("./Step4/results/moduleColors.dat")

# Helper: SYMBOL -> ENTREZID 
sym2entrez <- function(sym_vec) {
  m <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys     = unique(sym_vec),
    keytype  = "SYMBOL",
    column   = "ENTREZID",
    multiVals = "first"
  )
  out <- unname(m[sym_vec])   # respect input order
  out <- out[!is.na(out)]     # removed unmapped (e.g. LOC*)
  unique(out)                 # optional: remove duplicates
}

# Before the loop: universe in ENTREZ
universe_entrez <- sym2entrez(colnames(datExpr_filtered))

# ---- 8. Module annotation with GO/KEGG enrichment ----
# Create list to store enrichment results
module_enrichment <- list()

#unique_modules <- unique(moduleColors)
moduleLabels_num  <- as.integer(moduleLabels)  # 0,1,2,... by gene

for (m in modules_num_keep) {
  # Get genes in this module
  module_genes <- colnames(datExpr_filtered)[moduleLabels_num == as.integer(m)] #now labeled by number
  # A ENTREZ
  module_genes_entrez <- sym2entrez(module_genes)

  # GO enrichment
  ego <- tryCatch({
    enrichGO(
      gene = module_genes_entrez,
      universe = universe_entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )
  }, error = function(e) NULL)
  
  # KEGG enrichment
  ekegg <- tryCatch({
    enrichKEGG(
      gene = module_genes_entrez,
      universe = universe_entrez,
      organism = "hsa",
      keyType = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
  }, error = function(e) NULL)
  
  module_enrichment[[m]] <- list(
    GO = ego,
    KEGG = ekegg,
    n_genes = length(module_genes_entrez)
  )
  
  # Save results
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    write.csv(as.data.frame(ego), 
              paste0("./Step4/results/WGCNA_module_", m, "_GO.csv"),
              row.names = FALSE)
  }
  
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    write.csv(as.data.frame(ekegg), 
              paste0("./Step4/results/WGCNA_module_", m, "_KEGG.csv"),
              row.names = FALSE)
  }
}

################################################################################
# Fig.3: WGCNA - Weighted Gene Co-expression Network Analysis
# Purpose: Separate "structural" and "stress" genes into distinct modules
#          Confirm that structural and stress response genes operate as 
#          separate modules in the co-expression network
#This script is for:
#   5. Evaluate relationship between modules and Farm vs Wild / PC1 scores
################################################################################

# ---- 9. Create summary table of module annotations ----
module_summary <- data.frame(
  module = modules_num_keep,
  n_genes = sapply(modules_num_keep, function(m) 
    sum(moduleLabels_num == as.integer(m))),
  stringsAsFactors = FALSE
)

# Add top enriched terms
module_summary$top_GO_terms <- sapply(modules_num_keep, function(m) {
  ego <- module_enrichment[[m]]$GO
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    terms <- head(as.data.frame(ego)$Description, 3)
    paste(terms, collapse = "; ")
  } else {
    "No significant enrichment"
  }
})

# Add correlation with Farm vs Wild
module_summary$cor_Farm_vs_Wild <- moduleTraitCor[as.character(modules_num_keep), "Farm_vs_Wild"]
module_summary$pval_Farm_vs_Wild <- moduleTraitPvalue[as.character(modules_num_keep), "Farm_vs_Wild"]

# Add correlation with PC1 if available
if ("PC1" %in% colnames(moduleTraitCor)) {
  module_summary$cor_PC1 <- moduleTraitCor[as.character(modules_num_keep), "PC1"]
  module_summary$pval_PC1 <- moduleTraitPvalue[as.character(modules_num_keep), "PC1"]
}

# ---- 10. Save results ----
write.csv(module_summary, 
          "./Step4/results/WGCNA_module_summary.csv", 
          row.names = FALSE)

write.csv(MEs_with_info, 
          "./Step4/results/WGCNA_module_eigengenes.csv", 
          row.names = FALSE)

gene_assign <- data.frame(
  gene         = colnames(datExpr_filtered),
  module_num   = moduleLabels_num,
  module_color = moduleColors,
  stringsAsFactors = FALSE
)
write.csv(gene_assign, "./Step4/results/WGCNA_gene_module_assignment.csv", row.names = FALSE)

# Save workspace
save(net, MEs, moduleColors, moduleLabels_num, 
     module_enrichment, module_summary,
     file = "./Step4/results/WGCNA_results.RData")

print("Step 4 WGCNA analysis completed. Results saved in 'results/' directory.")

  ################################################################################
# Fig.3: Visualization of WGCNA results
# Purpose: Visualize co-expression modules and their relationships
# Methods: 
#   - Module dendrogram with color assignments
#   - Module-trait relationship heatmap
#   - Module eigengene boxplots (Farm vs Wild)
#   - Module eigengene correlation with PC1
#   - Network visualization for key modules
################################################################################

# Required libraries (if not already loaded)
library(WGCNA)
library(ggplot2)
library(pheatmap)
library(cowplot)
library(reshape2)
library(ComplexHeatmap)
library(circlize)

# Load WGCNA results if needed
# load("./Step4/results/WGCNA_results.RData")

# ---- 1. Module eigengene boxplots (Farm vs Wild comparison) ----
# Prepare data for plotting
MEs_long <- MEs_with_info %>%
  pivot_longer(
    cols = -c(sample, condition, site),
    names_to = "module",
    values_to = "eigengene"
  )

## Make sure 'module' in MEs_long is character for inspection
MEs_long$module <- sub("^ME", "", as.character(MEs_long$module))

# Color of the module according with the LUT you already created 
MEs_long$module_color <- unname(lut_label2color[MEs_long$module])

# Label for the facet: "color (number)"
MEs_long$module_label <- paste0(MEs_long$module_color, " (", MEs_long$module, ")")

## Create column 'module_color' robust:
##    - if 'module' is numeric (0,1,2,...) -> mape with labels2colors()
##    - if is text: remove "ME" if exists; if only number is left -> map; if color is left -> use
#MEs_long <- MEs_long %>%
#  mutate(tmp = sub("^ME", "", module),
#         module_color = case_when(
#           grepl("^[0-9]+$", tmp) ~ labels2colors(as.numeric(tmp)),
#           TRUE                   ~ tmp
#         )) %>%
#  select(-tmp)

# Select modules significantly correlated with Farm vs Wild
sig_modules <- module_summary %>%
  dplyr::filter(pval_Farm_vs_Wild < 0.05) %>%
  dplyr::pull(module)

sig_modules <- as.character(sig_modules)

MEs_sig <- dplyr::filter(MEs_long, module %in% sig_modules)

# Create boxplots for significant modules
p1 <- ggplot(MEs_sig, aes(x = condition, y = eigengene, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  facet_wrap(~ module, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Farm" = "#E69F00", "Wild" = "#56B4E9")) +
  labs(
    title = "Module eigengenes: Farm vs Wild comparison",
    subtitle = "Showing modules significantly correlated with Farm/Wild status",
    x = "Condition",
    y = "Module eigengene"
  ) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom"
  )

ggsave("./Step4/figures/Fig3A_module_eigengenes_farm_vs_wild.pdf",
       p1, width = 12, height = 10, dpi = 300)

# ---- 2. Plot dendrogram with module colors ----
pdf("./Step4/figures/Fig3A_WGCNA_dendrogram.pdf", width = 12, height = 6)

cols_block1 <- WGCNA::labels2colors(net$colors)[net$blockGenes[[1]]]

plotDendroAndColors(
  net$dendrograms[[1]],
  cols_block1,
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)
dev.off()

# ---- 3. Save figures ----
ggsave("./Step4/figures/Fig3A_module_eigengenes_farm_vs_wild.pdf",
       p1, width = 12, height = 10, dpi = 300)

# ---- 4. Site-specific patterns (boxplots by site) ----
p2 <- ggplot(MEs_sig, aes(x = site, y = eigengene, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  facet_wrap(~ module, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Farm" = "#E69F00", "Wild" = "#56B4E9")) +
  labs(
    title = "Module eigengenes by sampling site",
    x = "Site",
    y = "Module eigengene"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom"
  )

ggsave("./Step4/figures/Fig3B_module_eigengenes_by_site.pdf",
       p2, width = 14, height = 10, dpi = 300)

# ---- 5. Module-trait relationship heatmap ----
# Prepare data for heatmap
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

module_ids   <- rownames(moduleTraitCor)              # "1","6","15",...,"0"
module_cols  <- labels2colors(as.numeric(module_ids)) # color WGCNA of each id
ylabs        <- paste0(module_ids, " (", module_cols, ")")

pdf("./Step4/figures/Fig3B_module_trait_heatmap.pdf", width = 8, height = 10)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(moduleTraitCor),
  yLabels = ylabs,
  ySymbols = ylabs,
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = paste("Module-trait relationships")
)
dev.off()

# ---- 6. Module eigengene correlation with PC1 ----
if ("PC1" %in% colnames(traitData)) {
  # Create scatter plots for key modules
  key_modules <- head(sig_modules, 6)
  key_modules <- intersect(key_modules, colnames(MEs_with_info))
  
  MEs_pc1 <- MEs_with_info %>%
    select(sample, condition, all_of(key_modules)) %>%
    mutate(PC1 = traitData$PC1[match(sample, traitData$sample)]) %>%
    pivot_longer(
      cols = all_of(key_modules),
      names_to = "module",
      values_to = "eigengene"
    )
  
  p3 <- ggplot(MEs_pc1, aes(x = PC1, y = eigengene, color = condition)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
    facet_wrap(~ module, scales = "free_y", ncol = 3) +
    scale_color_manual(values = c("Farm" = "#E69F00", "Wild" = "#56B4E9")) +
    labs(
      title = "Module eigengenes vs PC1 scores",
      x = "PC1 score",
      y = "Module eigengene"
    ) +
    theme_classic() +
    theme(
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "bottom"
    )
}

ggsave("./Step4/figures/Fig3C_module_eigengenes_vs_PC1.pdf",
       p3, width = 12, height = 8, dpi = 300)

# ---- 7. Heatmap of module eigengenes ----
# Prepare annotation
annotation_col <- data.frame(
  Condition = sample_info$condition,
  Site = sample_info$site,
  row.names = sample_info$sample
)

# Color schemes
ann_colors <- list(
  Condition = c(Farm = "#E69F00", Wild = "#56B4E9"),
  Site = setNames(
    rainbow(length(unique(sample_info$site))),
    unique(sample_info$site)
  )
)

mod_ids  <- colnames(MEs)                                      # "1","6","15",...
num2col  <- setNames(labels2colors(as.numeric(mod_ids)), mod_ids)  # number -> color
row_lab  <- paste0(mod_ids, " (", num2col[mod_ids], ")")       # label: 1 (turquoise), etc.

# Create heatmap
p4 <- pheatmap(
  t(MEs),
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  show_rownames = TRUE,
  labels_row = row_lab,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Module eigengenes across all samples",
  fontsize = 10,
  fontsize_row = 8,
  filename = "figures/Fig3D_module_eigengene_heatmap.pdf",
  width = 10,
  height = 8
)

# ---- 8. Module enrichment summary plot ----
# Extract top GO terms for each module
module_go_summary <- lapply(sig_modules, function(mod) {
  ego <- module_enrichment[[mod]]$GO
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    df <- as.data.frame(ego) %>%
      head(5) %>%
      mutate(module = mod) %>%
      select(module, Description, p.adjust, Count)
    return(df)
  } else {
    return(NULL)
  }
})

module_go_df <- bind_rows(module_go_summary)

mods     <- unique(as.character(module_go_df$module))                # IDs "1","6",...
mod_cols <- setNames(labels2colors(as.numeric(mods)), mods)          # number -> color WGCNA

if (nrow(module_go_df) > 0) {
  p5 <- ggplot(module_go_df, 
               aes(x = reorder(Description, -log10(p.adjust)), 
                   y = -log10(p.adjust),
                   fill = module)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~ module, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = mod_cols) +
    labs(
      title = "Top enriched GO terms by module",
      x = "GO term",
      y = "-log10(adjusted p-value)"
    ) +
    theme_classic() +
    theme(
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )
}

if (exists("p5")) {
  ggsave("./Step4/figures/Fig3E_module_GO_enrichment.pdf",
         p5, width = 10, height = 12, dpi = 300)
}

################################################################################
# Fig.3: Visualization of WGCNA results
# This script is to visualize modules that we considered priority based
# on the category Farm vs Wild.
################################################################################

# Required libraries (if not already loaded)
library(WGCNA)
library(ggplot2)
library(pheatmap)
library(cowplot)
library(reshape2)
library(ComplexHeatmap)
library(circlize)

# Load WGCNA results if needed
# load("./Step4/results/WGCNA_results.RData")

# ---- 1. Prioritizing modules Farm/Wild ----

mods_all <- sort(unique(moduleLabels))
num2col  <- setNames(labels2colors(mods_all), mods_all)              # "11" -> "greenyellow"

# ---- 2. Module eigengene boxplots (Farm vs Wild comparison) ----
# Prepare data for plotting

lut_label2color <- setNames(
  labels2colors(sort(unique(moduleLabels))),
  as.character(sort(unique(moduleLabels)))
)

MEs_long <- MEs_with_info %>%
  pivot_longer(
    cols = -c(sample, condition, site),
    names_to = "module",
    values_to = "eigengene"
  )

MEs_long$module <- sub("^ME", "", as.character(MEs_long$module))
MEs_long$module_color <- unname(lut_label2color[MEs_long$module])
MEs_long$module_label <- paste0(MEs_long$module, " (", MEs_long$module_color, ")")

modules_priority <- module_summary %>%
  dplyr::filter(pval_Farm_vs_Wild < 0.05) %>%
  dplyr::arrange(pval_Farm_vs_Wild,
                 dplyr::desc(abs(cor_Farm_vs_Wild)),
                 dplyr::desc(n_genes))

sig_modules <- head(as.character(modules_priority$module), 9)

MEs_sig <- dplyr::filter(MEs_long, module %in% sig_modules)

p1 <- ggplot(MEs_sig, aes(x = condition, y = eigengene, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  facet_wrap(~ module_label, scales = "free_y", ncol = 3) +  # <-- cambiar aquí
  scale_fill_manual(values = c("Farm" = "#E69F00", "Wild" = "#56B4E9")) +
  labs(
    title = "Module eigengenes: Farm vs Wild comparison",
    subtitle = "Showing modules significantly correlated with Farm/Wild status",
    x = "Condition",
    y = "Module eigengene"
  ) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom"
  )

ggsave("./Step4/figures/Fig3A_module_eigengenes_farm_vs_wild(Prioritize).pdf",
       p1, width = 12, height = 10, dpi = 300)

# ---- 3. Heatmap módulo+rasgo ----

# Already generated:
# MEs_long$module <- sub("^ME", "", as.character(MEs_long$module))
# MEs_long$module_color <- unname(lut_label2color[MEs_long$module])
# MEs_long$module_label <- paste0(MEs_long$module_color, " (", MEs_long$module, ")")

## (if not generated) build the LUT number->color, only once
if (!exists("lut_label2color")) {
  map_df <- unique(data.frame(label = moduleLabels,
                              color = moduleColors,
                              stringsAsFactors = FALSE))
  lut_label2color <- setNames(map_df$color, as.character(map_df$label))
}

## significant modules by Farm/Wild (PRIORITY)
sig_modules <- module_summary %>%
  dplyr::filter(pval_Farm_vs_Wild < 0.05) %>%
  dplyr::arrange(pval_Farm_vs_Wild,
                 dplyr::desc(abs(cor_Farm_vs_Wild)),
                 dplyr::desc(n_genes)) %>%
  dplyr::pull(module) %>%
  as.character()

## (optional) push those that are not associated with PC1
if ("pval_PC1" %in% colnames(module_summary)) {
  ord_df <- module_summary %>%
    dplyr::mutate(overlaps_PC1 = pval_PC1 < 0.05) %>%
    dplyr::filter(module %in% sig_modules) %>%
    dplyr::arrange(overlaps_PC1,          # FALSE primero
                   pval_Farm_vs_Wild,
                   dplyr::desc(abs(cor_Farm_vs_Wild)),
                   dplyr::desc(n_genes))
  sig_modules <- as.character(ord_df$module)
}

## 2) top N (for facet) and make sure they are in your MEs
sig_modules_top <- head(sig_modules, 9)
# Attention: your columns of MEs are "ME1","ME2",... and in MEs_long we already removed "ME" -> "1","2",...
sig_modules_top <- intersect(sig_modules_top, unique(MEs_long$module))

## 3) Filter and build the labels
MEs_sig <- dplyr::filter(MEs_long, module %in% sig_modules_top)

MEs_sig$module_label <- paste0(
  unname(lut_label2color[MEs_sig$module]), " (", MEs_sig$module, ")"
)

MEs_sig$module_label <- factor(
  MEs_sig$module_label,
  levels = paste0(unname(lut_label2color[sig_modules_top]), " (", sig_modules_top, ")")
)

# ---- 4. Site-specific patterns (boxplots by site) ----
p2 <- ggplot(MEs_sig,
             aes(x = factor(site, levels = unique(sample_info$site)),  # <-- CHANGE
                 y = eigengene, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  facet_wrap(~ module_label, scales = "free_y", ncol = 3) +            # <-- CHANGE
  scale_fill_manual(values = c("Farm" = "#E69F00", "Wild" = "#56B4E9")) +
  labs(
    title = "Module eigengenes by sampling site",
    x = "Site",
    y = "Module eigengene"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom"
  )

ggsave("./Step4/figures/Fig3B_module_eigengenes_by_site_(Prioritize).pdf",
       p2, width = 14, height = 10, dpi = 300)

# ---- 5. Farm/Wild effect summary for key modules ----
# Use the modules prioritized in the previous step
key_modules <- as.character(sig_modules_top)

df_fw <- module_summary %>%
  dplyr::filter(module %in% key_modules) %>%
  dplyr::transmute(
    module      = as.character(module),
    cor_fw      = cor_Farm_vs_Wild,
    p_fw        = pval_Farm_vs_Wild,
    color       = unname(lut_label2color[as.character(module)]),    # e.g. "salmon"
    facet_label = paste0(color, " (", module, ")")                  # "salmon (11)"
  )

# order by |cor| (major effect first)
df_fw$facet_label <- factor(
  df_fw$facet_label,
  levels = df_fw$facet_label[order(-abs(df_fw$cor_fw))]
)

# optional: stars of significance
stars <- function(p) ifelse(p < 0.001, "***",
                            ifelse(p < 0.01, "**",
                                   ifelse(p < 0.05, "*", "")))
df_fw$signif <- stars(df_fw$p_fw)

# Plot: barplot of correlations in modules ~ Farm/Wild with p-value
p3 <- ggplot(df_fw, aes(x = facet_label, y = cor_fw, fill = cor_fw)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_col(alpha = 0.85) +
  geom_text(aes(label = paste0("p=", signif(p_fw, 2), " ", signif)),
            hjust = ifelse(df_fw$cor_fw >= 0, -0.05, 1.05),
            size = 3) +
  coord_flip() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, guide = "none") +
  labs(
    title = "Farm vs Wild effect by module (prioritized)",
    x = "Module (color [number])",
    y = "Correlation with Farm_vs_Wild"
  ) +
  theme_classic() +
  theme(
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.title.x = element_text(margin = margin(t = 8))
  )

ggsave("./Step4/figures/Fig3C_module_FarmWild_effect_summary.pdf",
       p3, width = 8, height = 6, dpi = 300)

# ---- 6. GO/KEGG enrichment por módulo ----
# use the already defined:
# module_enrichment (from your loop)
# sig_modules        (all the significant FW)
# sig_modules_top    (e.g. top 9 for facet)
# lut_label2color    (number -> color), num2col equivalent if you prefer

# --- Compact GO only for modules FW ---
module_go_summary_fw <- lapply(sig_modules, function(mod) {
  ego <- module_enrichment[[mod]]$GO
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    as.data.frame(ego) |>
      head(5) |>
      dplyr::mutate(
        module = as.character(mod),
        module_color = unname(lut_label2color[as.character(mod)]),
        module_label = paste0(module_color, " (", module, ")")
      ) |>
      dplyr::select(module, module_color, module_label, Description, p.adjust, Count)
  } else {
    NULL
  }
})

module_go_df <- dplyr::bind_rows(module_go_summary_fw)

# order facets based on FW priority (sig_modules_top first)
if (nrow(module_go_df) > 0) {
  module_go_df$module <- factor(module_go_df$module, levels = as.character(sig_modules))
  module_go_df$module_label <- factor(
    module_go_df$module_label,
    levels = paste0(unname(lut_label2color[as.character(sig_modules)]), " (", sig_modules, ")")
  )
  # colors by number of module
  mods     <- levels(module_go_df$module)
  mod_cols <- setNames(labels2colors(as.numeric(mods)), mods)
  
  p5 <- ggplot(module_go_df, aes(x = reorder(Description, -log10(p.adjust)),
                                 y = -log10(p.adjust),
                                 fill = module)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~ module_label, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = mod_cols, guide = "none") +
    labs(title = "Top enriched GO terms (FW-prioritized modules)",
         x = "GO term", y = "-log10(adjusted p-value)") +
    theme_classic() +
    theme(strip.background = element_rect(fill = "white", color = "black"),
          strip.text = element_text(face = "bold"))
  
  ggsave("./Step4/figures/Fig3E_module_GO_enrichment_FW.pdf",
         p5, width = 10, height = 12, dpi = 300)
  
  # (optional) save the combined filtered table
  readr::write_csv(module_go_df, "./Step4/results/WGCNA_GO_enrichment_FW_only.csv")
}

module_kegg_summary_fw <- lapply(sig_modules, function(mod) {
  ekegg <- module_enrichment[[mod]]$KEGG
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    as.data.frame(ekegg) |>
      head(5) |>
      dplyr::mutate(
        module = as.character(mod),
        module_color = unname(lut_label2color[as.character(mod)]),
        module_label = paste0(module_color, " (", module, ")")
      ) |>
      dplyr::select(module, module_color, module_label, Description, p.adjust, Count)
  } else {
    NULL
  }
})

module_kegg_df <- dplyr::bind_rows(module_kegg_summary_fw)

if (nrow(module_kegg_df) > 0) {
  module_kegg_df$module <- factor(module_kegg_df$module, levels = as.character(sig_modules))
  module_kegg_df$module_label <- factor(
    module_kegg_df$module_label,
    levels = paste0(unname(lut_label2color[as.character(sig_modules)]), " (", sig_modules, ")")
  )
  mods     <- levels(module_kegg_df$module)
  mod_cols <- setNames(labels2colors(as.numeric(mods)), mods)
  
  p5k <- ggplot(module_kegg_df, aes(x = reorder(Description, -log10(p.adjust)),
                                    y = -log10(p.adjust),
                                    fill = module)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~ module_label, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = mod_cols, guide = "none") +
    labs(title = "Top enriched KEGG pathways (FW-prioritized modules)",
         x = "KEGG pathway", y = "-log10(adjusted p-value)") +
    theme_classic() +
    theme(strip.background = element_rect(fill = "white", color = "black"),
          strip.text = element_text(face = "bold"))
  
  ggsave("./Step4/figures/Fig3E_module_KEGG_enrichment_FW.pdf",
         p5k, width = 10, height = 12, dpi = 300)
  
  readr::write_csv(module_kegg_df, "./Step4/results/WGCNA_KEGG_enrichment_FW_only.csv")
}
