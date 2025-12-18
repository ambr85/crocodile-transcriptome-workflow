################################################################################
# ssGSEA - Single Sample Gene Set Enrichment Analysis
# Methods: 
#   1. Get Raw ssGSEA scores
#   2. Heatmaps of all gene sets
#   3. Perform module-trait correlations
#   4. Generate statistical tests for biological interpretation
################################################################################

# ========= Download gene sets from MSigDB ====================
# cat("Files generated:\n")
# cat("1. ssGSEA_scores.csv - Raw ssGSEA scores\n")
# cat("2. MSigDB_selected_genesets_human.gmt - Downloaded gene sets\n")
# cat("3. Fig4A_ssGSEA_boxplot_by_site.pdf - Boxplots by site\n")
# cat("4. Fig4B_ssGSEA_heatmap.pdf - Heatmap of all gene sets\n")
# cat("5. Fig4C_Module_ssGSEA_correlation.pdf - Module-trait correlation\n")
# cat("6. Fig4D_Module_ssGSEA_scatter_example.pdf - Example scatter plot\n")
# cat("7. ssGSEA_category_statistics.csv - Statistical test results\n")
# cat("8. Module_ssGSEA_correlation.csv - Correlation values\n")
# cat("9. Strong_Module_ssGSEA_correlations.csv - Significant correlations\n\n")
# 
# cat("Next steps:\n")
# cat("- Review statistical results in ssGSEA_category_statistics.csv\n")
# cat("- Identify which categories show strongest differences\n")
# cat("- Check Module-ssGSEA correlations for biological interpretation\n")
# ==================================================

library(msigdbr)

# ---------- 0. prepare ------------
# check the species in the database
msigdbr_species()

# get Human gene set（for ortholog conversion）
msigdb_human <- msigdbr(species = "Homo sapiens")

# define gene set
target_genesets <- list(
  
  # 1. Oxidative stress
  "GO_RESPONSE_TO_OXIDATIVE_STRESS" = c("C5", "GO:BP"),
  "GO_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS" = c("C5", "GO:BP"),
  "GO_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS" = c("C5", "GO:BP"),
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION" = c("H", NA),
  
  # 2. Xenobiotic / Chemical detox
  "KEGG_METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450" = c("C2", "CP:KEGG"),
  "KEGG_DRUG_METABOLISM_CYTOCHROME_P450" = c("C2", "CP:KEGG"),
  "KEGG_GLUTATHIONE_METABOLISM" = c("C2", "CP:KEGG"),
  "GO_XENOBIOTIC_METABOLIC_PROCESS" = c("C5", "GO:BP"),
  "GO_RESPONSE_TO_TOXIC_SUBSTANCE" = c("C5", "GO:BP"),
  
  # 3. Immune / Inflammatory
  "HALLMARK_INFLAMMATORY_RESPONSE" = c("H", NA),
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB" = c("H", NA),
  "HALLMARK_INTERFERON_ALPHA_RESPONSE" = c("H", NA),
  "HALLMARK_INTERFERON_GAMMA_RESPONSE" = c("H", NA),
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY" = c("C2", "CP:KEGG"),
  
  # 4. Epidermal / Barrier
  "GO_KERATINOCYTE_DIFFERENTIATION" = c("C5", "GO:BP"),
  "GO_EPIDERMAL_CELL_DIFFERENTIATION" = c("C5", "GO:BP"),
  "GO_SKIN_DEVELOPMENT" = c("C5", "GO:BP"),
  "GO_CORNIFIED_ENVELOPE" = c("C5", "GO:CC"),
  
  # 5. Mitochondrial OXPHOS
  "KEGG_OXIDATIVE_PHOSPHORYLATION" = c("C2", "CP:KEGG"),
  "KEGG_CITRATE_CYCLE_TCA_CYCLE" = c("C2", "CP:KEGG"),
  "GO_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT" = c("C5", "GO:BP"),
  
  # 6. Translation / Ribosome
  "KEGG_RIBOSOME" = c("C2", "CP:KEGG"),
  "GO_TRANSLATION" = c("C5", "GO:BP"),
  "GO_RIBOSOME_BIOGENESIS" = c("C5", "GO:BP")
)

# ----------- 1. extraction of gene set ---------------------
extract_geneset <- function(geneset_name, collection_subcollection) {
  collection <- collection_subcollection[1]
  subcollection <- collection_subcollection[2]
  
  if (is.na(subcollection)) {
    # for Hallmark
    genes <- msigdb_human %>%
      filter(gs_collection == collection,
             gs_name == geneset_name) %>%
      pull(gene_symbol) %>%
      unique()
  } else {
    # other
    genes <- msigdb_human %>%
      filter(gs_collection == collection,
             gs_subcollection == subcollection,
             gs_name == geneset_name) %>%
      pull(gene_symbol) %>%
      unique()
  }
  
  if (length(genes) == 0) {
    warning(paste("No genes found for:", geneset_name))
  } else {
    cat(geneset_name, ":", length(genes), "genes\n")
  }
  
  return(genes)
}

# get whole gene sets 
genesets_human <- lapply(names(target_genesets), function(gs_name) {
  extract_geneset(gs_name, target_genesets[[gs_name]])
})
names(genesets_human) <- names(target_genesets)

# remove 0
genesets_human <- genesets_human[sapply(genesets_human, length) > 0]

cat("\n=== Downloaded", length(genesets_human), "gene sets from MSigDB ===\n")
print(sapply(genesets_human, length))

# ------------- 2. save the data -----------------------------
write_gmt <- function(genesets, filename) {
  gmt_lines <- lapply(names(genesets), function(gs_name) {
    genes <- genesets[[gs_name]]
    paste(c(gs_name, "description", genes), collapse = "\t")
  })
  writeLines(unlist(gmt_lines), filename)
}

write_gmt(genesets_human, "./YOUR_DIRECTORY/MSigDB_selected_genesets_human.gmt")

# ======== Conversion Crocodile → Alligator → Human ortholog ===========
# your data（Organism gene ID）
# expression_matrix: gene（row）× sample（colum）
# rownames = Organism gene IDs
# Example = Alligator mississippiensis
# ==============================================================================
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("biomaRt")

library(biomaRt)

# ------------- 3. convert ID -----------------------------
# check species
# listDatasets(useEnsembl("ensembl")) %>% filter(grepl("croc|allig", description, ignore.case=TRUE))

# for Alligator mississippiensis #there is no alligator in ensembl. Modify accordingly.
tryCatch({
  # Alligator database
  ensembl_alligator <- useEnsembl("ensembl", 
                                  dataset = "amississippiensis_gene_ensembl")
  
  # Human database
  ensembl_human <- useEnsembl("ensembl", 
                              dataset = "hsapiens_gene_ensembl")
  
  # Alligator gene IDs
  crocodile_genes <- rownames(expression_matrix)
  
  # Alligator → Human ortholog
  orthologs <- getLDS(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "ensembl_gene_id",
    values = crocodile_genes,
    mart = ensembl_alligator,
    attributesL = c("ensembl_gene_id", "hgnc_symbol"),
    martL = ensembl_human,
    uniqueRows = TRUE
  )
  
  colnames(orthologs) <- c("Crocodile_ID", "Crocodile_Symbol", 
                           "Human_Ensembl", "Human_Symbol")
  
  cat("Orthologs found:", nrow(orthologs), "\n")
  
}, error = function(e) {

  # --- option 2: Maching by gene symbol ---

  orthologs <- data.frame(
    Crocodile_Symbol = rownames(expression_matrix),
    Human_Symbol = rownames(expression_matrix),  
    stringsAsFactors = FALSE
  )
  
  cat("Using direct gene symbol matching\n")
})

# ------------- 4. handling one-to-many mappings ----------------- 
# When a single Crocodile gene maps to multiple human orthologs,  
# select the first one.  

orthologs_unique <- orthologs %>%
  group_by(Crocodile_Symbol) %>%
  slice(1) %>%  
  ungroup()

# ------------- 5. convert row and column names -----------------
expression_matrix_human <- expression_matrix[
  rownames(expression_matrix) %in% orthologs_unique$Crocodile_Symbol, 
]

# convert row names to Human symbol
rownames_mapping <- setNames(
  orthologs_unique$Human_Symbol,
  orthologs_unique$Crocodile_Symbol
)

rownames(expression_matrix_human) <- rownames_mapping[rownames(expression_matrix_human)]

save(expression_matrix_human, file = "./YOUR_DIRECTORY/expression_matrix_human.dat")

###########################################################################
### ALTERNATIVE TO CONVERSION (USED FOR OUR ANALYSIS) WITHOUT USING biomaRt
###########################################################################
# Our current matrix with SYMBOL genes, already used for WGCNA
expr_sym <- t(datExpr_filtered)

# Clean the names of the sample: remove suffix "countREV.out"
colnames(expr_sym) <- sub("countREV\\.out$", "", colnames(expr_sym))

# Collapse duplicated genes by symbol (mean); change to sum/median if you prefer
expr_sym <- as.data.frame(expr_sym)
expr_sym$Gene <- rownames(expr_sym)
expression_matrix_human <- expr_sym |>
  dplyr::group_by(Gene) |>
  dplyr::summarise(dplyr::across(dplyr::everything(), mean), .groups = "drop") |>
  tibble::column_to_rownames("Gene") |>
  as.matrix()

# Save with the SAME name/route expected for your process
save(expression_matrix_human, file = "./YOUR_DIRECTORY/expression_matrix_human.dat")

# ============== MAKE EXPRESSION PROFILE ====================
# your data（Crocodile gene ID）
# expression_matrix: gene（row）× sample（colum）
# rownames = Crocodile gene IDs
# ==============================================================================

library(dplyr)
library(tibble)

# ------------- 0. Load data ------------------
load("./YOUR_DIRECTORY/expression_matrix_human.dat")

# ------------- 3. Make expression matrix -----------------
expression_matrix_human <- as.matrix(
  expression_matrix_human %>%
    as.data.frame() %>%
    mutate(Gene = rownames(.)) %>%
    group_by(Gene) %>%
    summarise(across(everything(), mean)) %>%
    column_to_rownames("Gene")
)

cat("Genes after ortholog conversion:", nrow(expression_matrix_human), "\n")

# ------------- 4. filtering expression matrix -----------------
#Filter the gene sets to include only genes that are present in the Crocodile expression data.
available_genes <- rownames(expression_matrix_human)

genesets_filtered <- lapply(genesets_human, function(genes) {
  intersect(genes, available_genes)
})

# Remove gene sets that contain fewer than five genes.
genesets_filtered <- genesets_filtered[sapply(genesets_filtered, length) >= 5]

cat("\n=== Gene sets after filtering ===\n")
print(sapply(genesets_filtered, length))

# Check the final number of gene sets.
cat("\nFinal number of gene sets:", length(genesets_filtered), "\n")

save(genesets_filtered, file = "./YOUR_DIRECTORY/genesets_filtered.dat")

# ========== Running ssGSEA ====================
#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("GSVA")

library(GSVA)
library(BiocParallel)
set.seed(123)

packageVersion("GSVA")

# ------------- 0. Load data ------------------
load("./YOUR_DIRECTORY/expression_matrix_human.dat")
load("./YOUR_DIRECTORY/genesets_filtered.dat")

# ------------- 5. run ssGSEA ------------------
x <- expression_matrix_human

param <- ssgseaParam(x,  # matrix: genes x samples (SYMBOLS humans)
  geneSets      = genesets_filtered,        # list of genesets (SYMBOLS humans)
  minSize       = 5,
  maxSize       = 500
)

# Parallel (adjust to your PC)
bp <- SnowParam(workers = 4)  # use in Windows MulticoreParam(..., type="SOCK")

# --- Execute ssGSEA ---
ssgsea_scores <- gsva(param, BPPARAM = bp, verbose = TRUE)

cat("ssGSEA completed!\n")
cat("Output dimensions:", dim(ssgsea_scores), "\n")
cat("Gene sets:", nrow(ssgsea_scores), "\n")
cat("Samples:", ncol(ssgsea_scores), "\n")
cat("\n=== ssGSEA score summary ===\n")
print(summary(as.vector(ssgsea_scores)))

# ------------- 6. Save data ------------------
write.csv(ssgsea_scores, "./YOUR_DIRECTORY/ssGSEA_scores.csv")

# ====== Group the gene sets by category (for visualization) ==========

# ------------- 0. Load data ------------------
load("./YOUR_DIRECTORY/genesets_filtered.dat")
 
# ------------- 7. Classify the gene set into category ------------------
geneset_categories <- data.frame(
  Geneset = names(genesets_filtered),
  Category = case_when(
    grepl("OXIDATIVE|ROS|REACTIVE_OXYGEN", names(genesets_filtered)) ~ "Oxidative Stress",
    grepl("XENOBIOTIC|DRUG_METABOLISM|GLUTATHIONE|TOXIC", names(genesets_filtered)) ~ "Xenobiotic/Detox",
    grepl("INFLAMMATORY|TNFA|INTERFERON|TOLL_LIKE|CYTOKINE", names(genesets_filtered)) ~ "Immune/Inflammatory",
    grepl("KERATINOCYTE|EPIDERMAL|SKIN|CORNIFIED", names(genesets_filtered)) ~ "Epidermal/Barrier",
    grepl("OXIDATIVE_PHOSPHORYLATION|CITRATE_CYCLE|ATP_SYNTHESIS", names(genesets_filtered)) ~ "Mitochondrial OXPHOS",
    grepl("RIBOSOME|TRANSLATION|RRNA", names(genesets_filtered)) ~ "Translation/Ribosome",
    TRUE ~ "Other"
  ),
  stringsAsFactors = FALSE
)

# ------------- 8. Calculate score of category ------------------
ssgsea_by_category <- sapply(unique(geneset_categories$Category), function(cat) {
  genesets_in_cat <- geneset_categories$Geneset[geneset_categories$Category == cat]
  if (length(genesets_in_cat) == 1) {
    ssgsea_scores[genesets_in_cat, ]
  } else {
    colMeans(ssgsea_scores[genesets_in_cat, , drop=FALSE])
  }
})

ssgsea_by_category <- t(ssgsea_by_category)

cat("\n=== Category-level scores calculated ===\n")
print(dim(ssgsea_by_category))

# ------------- 9. Save data ------------------
save(ssgsea_by_category, file = "./YOUR_DIRECTORY/ssgsea_by_category.dat")

# ======== Visualization by sites (ssGSEA score boxplot by sites) ========
library(ggplot2)
library(tidyr)
library(dplyr)
library(rstatix)

# ------------- 0. Load data ------------------
load("./YOUR_DIRECTORY/ssgsea_scores")

# Prepare the metadata (adjust as needed for your dataset).
# metadata: 
# - Sample: sample name
# - Site: site name（Farm, Wild_Site1, Wild_Site2, ...）
# - Group: Farm or Wild

# Use the data already created in: sample_info
meta_df <- sample_info |>
  dplyr::select(sample, condition, site) |>
  dplyr::rename(Sample = sample) |>
  as.data.frame()

# Make sure the names of samples and ssGSEA are the same (no 'countREV.out')
meta_df$Sample <- sub("countREV\\.out$", "", meta_df$Sample)
colnames(ssgsea_scores) <- sub("countREV\\.out$", "", colnames(ssgsea_scores))

# If you already have a vector named as: c("REACTOME_X"="Xenobiotic", ...)
# convert to tibble with the columns Geneset and Category
gs_cat_df <- tibble::tibble(
  Geneset  = names(geneset_categories),
  Category = as.character(unname(geneset_categories))
)

# ------------- 10. Convert ssGSEA score for long format ------------------
ssgsea_long <- ssgsea_scores %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Geneset") %>%
  tidyr::pivot_longer(-Geneset, names_to = "Sample", values_to = "ssGSEA_score") %>%
  dplyr::left_join(meta_df, by = "Sample") %>%
  dplyr::left_join(gs_cat_df, by = "Geneset")

# ------------- 11. Create the category-level data ------------------
ssgsea_category_long <- ssgsea_by_category %>%
  as.data.frame() %>%
  rownames_to_column("Category") %>%
  pivot_longer(-Category, names_to = "Sample", values_to = "Category_score") %>%
  left_join(meta_df, by = "Sample")

# ------------- 12. Make boxplot (by category) ------------------
fig4a <- ggplot(ssgsea_category_long, 
                aes(x = site, y = Category_score, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  facet_wrap(~ Category, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Farm" = "red", "Wild" = "lightblue")) +
  labs(title = "Category-level ssGSEA scores by site",
       x = "Site",
       y = "ssGSEA score") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "gray90"),
        legend.position = "top")

ggsave("./YOUR_DIRECTORY/Fig4A_ssGSEA_boxplot_by_site.pdf", fig4a, width = 12, height = 10)

# ------------- 13. Statistical test (Farm vs Wild) ------------------
stat_results <- ssgsea_category_long %>%
  group_by(Category) %>%
  wilcox_test(Category_score ~ condition) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

write.csv(stat_results, "./YOUR_DIRECTORY/ssGSEA_category_statistics.csv", row.names = FALSE)

cat("\n=== Statistical tests (Farm vs Wild) ===\n")
print(stat_results)

# ======= Visualization by sites (ssGSEA score heatmap) =========
library(ggplot2)
library(tidyr)
library(dplyr)
library(rstatix)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# ------------- 0. Load data ------------------
load("./YOUR_DIRECTORY/ssgsea_scores")
# Prepare the metadata (adjust as needed for your dataset).
# metadata: 
# - Sample: sample name
# - Site: site name（Farm, Wild_Site1, Wild_Site2, ...）
# - Group: Farm or Wild

# ------------- 14. Prepare annotation of Heatmap (geneset x sample) ------------------
# sampleannotation
col_anno <- HeatmapAnnotation(
  Site = meta_df$site,
  Group = meta_df$condition,
  col = list(
    Group = c("Farm" = "red", "Wild" = "lightblue"),
    Site = setNames(
      colorRampPalette(brewer.pal(8, "Set2"))(length(unique(meta_df$site))),
      unique(meta_df$site)
    )
  ),
  annotation_name_side = "left"
)

# Vector of categories aligned to the rows of the heatmap
cats <- gs_cat_df$Category[ match(rownames(ssgsea_scores), gs_cat_df$Geneset) ]
cats[is.na(cats)] <- "Uncategorized"   # avoids NA 

# Palette named by level -> color
levs <- unique(cats)
n <- length(levs)

pal <- if (n <= 3) {
  RColorBrewer::brewer.pal(3, "Dark2")[seq_len(n)]
} else if (n <= 8) {
  RColorBrewer::brewer.pal(n, "Dark2")
} else {
  grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(n)
}
names(pal) <- levs

# 3) Anotación de filas
row_anno <- ComplexHeatmap::rowAnnotation(
  Category = cats,
  col = list(Category = pal)
)

# ------------- 15. Make Heatmap ------------------
fig4b <- Heatmap(
  ssgsea_scores,
  name = "ssGSEA\nscore",
  top_annotation = col_anno,
  left_annotation = row_anno,
  
  # Clustering
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",

  col = colorRamp2(
    c(min(ssgsea_scores), 0, max(ssgsea_scores)),
    c("blue", "white", "red")
  ),
  
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  row_title = "Gene Sets",
  column_title = "Samples",
  
  column_split = meta_df$site,
  column_gap = unit(2, "mm"),
  
  # Legend
  heatmap_legend_param = list(
    title = "ssGSEA score",
    direction = "vertical",
    title_position = "leftcenter-rot"
  )
)

pdf("./YOUR_DIRECTORY/Fig4B_ssGSEA_heatmap.pdf", width = 14, height = 10)
draw(fig4b)
dev.off()

# ========= Correlation with WGCNA module eigengenes ===========
library(WGCNA)

# ------------- 0. Load data ------------------
# WGCNA results (assuming already computed)
# MEs: module eigengenes (samples × modules)
# moduleColors: module assignments for each gene

# ------------- 16. Preparations ------------------
# 0) Helper to normalize IDs of samples
.fix <- function(v) sub("countREV\\.out$", "", as.character(v))

# 1) Align sample IDs in both objects
#    - MEs: normally has rownames = sample IDs
#    - ssgsea_traits: transposed by samples x categories
if (is.null(rownames(MEs))) {
  # if MEs are in MEs_with_info: set rownames = sample and extract only eigengenes
  stop("rownames(MEs) es NULL. Use rownames(MEs) <- MEs_with_info$sample and leave only columns EIGEN.")
}

rownames(MEs)            <- .fix(rownames(MEs))
rownames(ssgsea_traits)  <- .fix(rownames(ssgsea_traits))

# 2) Subsets with common samples
common_samples <- intersect(rownames(MEs), rownames(ssgsea_traits))
if (length(common_samples) == 0) {
  stop("There are no common samples after standarize IDs. Check head(rownames(MEs)) and head(rownames(ssgsea_traits)).")
}

MEs_subset            <- MEs[common_samples, , drop = FALSE]
ssgsea_traits_subset  <- ssgsea_traits[common_samples, , drop = FALSE]

# 3) Make sure they are numeric
MEs_subset           <- data.matrix(MEs_subset)
ssgsea_traits_subset <- data.matrix(ssgsea_traits_subset)

# Use ssGSEA scores as traits
# Use category-level scores (easier to interpret)
#ssgsea_traits <- t(ssgsea_by_category)
#colnames(ssgsea_traits) <- paste0("ssGSEA_", colnames(ssgsea_traits))

# Align the sample order
#common_samples <- intersect(rownames(MEs), rownames(ssgsea_traits))
#MEs_subset <- MEs[common_samples, ]
#ssgsea_traits_subset <- ssgsea_traits[common_samples, ]

# ------------- 17. Calculate Module-ssGSEA correlation ------------------
moduleTraitCor_ssGSEA <- cor(MEs_subset, ssgsea_traits_subset, use = "pairwise.complete.obs")
moduleTraitPvalue_ssGSEA  <- WGCNA::corPvalueStudent(moduleTraitCor_ssGSEA, nSamples = nrow(MEs_subset))
#moduleTraitCor_ssGSEA <- cor(MEs_subset, ssgsea_traits_subset, use = "p")
nSamples <- nrow(MEs_subset)
#moduleTraitPvalue_ssGSEA <- corPvalueStudent(moduleTraitCor_ssGSEA, nSamples)

# FDR
moduleTraitPvalue_adj_ssGSEA <- matrix(
  p.adjust(as.vector(moduleTraitPvalue_ssGSEA), method = "BH"),
  nrow = nrow(moduleTraitPvalue_ssGSEA),
  dimnames = dimnames(moduleTraitPvalue_ssGSEA)
)

# maker
signif_labels_ssGSEA <- ifelse(moduleTraitPvalue_adj_ssGSEA < 0.001, "***",
                               ifelse(moduleTraitPvalue_adj_ssGSEA < 0.01, "**",
                                      ifelse(moduleTraitPvalue_adj_ssGSEA < 0.05, "*", "")))

# Correct labels (rows = modules = columns of MEs_subset)
yLabs <- colnames(MEs_subset)
if (is.null(yLabs)) yLabs <- rownames(moduleTraitCor_ssGSEA)  # safe fallback

# Optional labels: number + color if you have LUT
# (you can # this part if you are not using lut_label2color)
ySyms <- yLabs
ySyms_num <- gsub("^ME", "", yLabs)  # if columns are "ME1","ME2",...
if (exists("lut_label2color")) {
  ySyms <- paste0(ySyms_num, " (", unname(lut_label2color[ySyms_num]), ")")
}

# Make sure the matrix of the stars have the same dimensions
stopifnot(identical(dim(moduleTraitCor_ssGSEA), dim(signif_labels_ssGSEA)))

# ------------- 18. Module-ssGSEA correlation Heatmap ------------------
pdf("./YOUR_DIRECTORY/Fig4C_Module_ssGSEA_correlation.pdf", width = 10, height = 8)

labeledHeatmap(
  Matrix = moduleTraitCor_ssGSEA,
  xLabels = colnames(ssgsea_traits_subset),
  yLabels = yLabs,
  ySymbols = ySyms,
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = signif_labels_ssGSEA,
  setStdMargins = FALSE,
  cex.text = 1.0,
  cex.lab.x = 0.5,
  cex.lab.y = 0.9,
  zlim = c(-1, 1),
  main = "Module-ssGSEA Category Correlation"
)

dev.off()

# Save files
write.csv(moduleTraitCor_ssGSEA, "./YOUR_DIRECTORY/Module_ssGSEA_correlation.csv")
write.csv(moduleTraitPvalue_adj_ssGSEA, "./YOUR_DIRECTORY/Module_ssGSEA_pvalue_adj.csv")

# ============ Correlation with WGCNA module eigengenes ============
library(WGCNA)

# ------------- 0. Load data ------------------
#moduleTraitPvalue_adj_ssGSEA <- read.csv("./YOUR_DIRECTORY/Module_ssGSEA_pvalue_adj.csv")
#moduleTraitCor_ssGSEA <- read.csv("./YOUR_DIRECTORY/Module_ssGSEA_correlation.csv")

# ------------- 19. Extract strong correlations ------------------
# strong correlation（|r| > 0.6, padj < 0.05）
strong_correlations <- which(
  abs(moduleTraitCor_ssGSEA) > 0.6 & moduleTraitPvalue_adj_ssGSEA < 0.05,
  arr.ind = TRUE
)

if (nrow(strong_correlations) > 0) {
  cat("\n=== Strong Module-ssGSEA correlations ===\n")
  strong_cor_df <- data.frame(
    Module = rownames(moduleTraitCor_ssGSEA)[strong_correlations[, 1]],
    ssGSEA_Category = colnames(moduleTraitCor_ssGSEA)[strong_correlations[, 2]],
    Correlation = moduleTraitCor_ssGSEA[strong_correlations],
    P_adj = moduleTraitPvalue_adj_ssGSEA[strong_correlations]
  ) %>%
    arrange(desc(abs(Correlation)))
  
  print(strong_cor_df)
  write.csv(strong_cor_df, "./YOUR_DIRECTORY/Strong_Module_ssGSEA_correlations.csv", row.names = FALSE)
}

# ------------- 20. Make scatterplots of key modules and ssGSEA categories ------------------
# Example: 21 module vs. Xenobiotic/Detox ssGSEA
if ("21" %in% colnames(MEs) && "ssGSEA_Xenobiotic/Detox" %in% colnames(ssgsea_traits_subset)) {
  
  plot_data <- data.frame(
    Module_eigengene = MEs_subset[, "21"],
    ssGSEA_score = ssgsea_traits_subset[, "ssGSEA_Xenobiotic/Detox"],
    Site = metadata$Site[match(rownames(MEs_subset), metadata$Sample)],
    Group = metadata$Group[match(rownames(MEs_subset), metadata$Sample)]
  )
  
  cor_test <- cor.test(plot_data$Module_eigengene, plot_data$ssGSEA_score, 
                       method = "spearman")
  
  fig4d <- ggplot(plot_data, aes(x = Module_eigengene, y = ssGSEA_score, 
                                 color = Group)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    scale_color_manual(values = c("Farm" = "red", "Wild" = "lightblue")) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
             label = paste0("Spearman's rho = ", round(cor_test$estimate, 3),
                            "\np = ", format.pval(cor_test$p.value, digits = 2))) +
    labs(title = "Module 21 vs Xenobiotic/Detox ssGSEA",
         x = "Module 21 eigengene",
         y = "Xenobiotic/Detox ssGSEA score") +
    theme_bw(base_size = 12)
  
  ggsave("./YOUR_DIRECTORY/Fig4D_Module_ssGSEA_scatter_example.pdf", fig4d, width = 7, height = 6)
}

###############################################
## EXTRA – ADD Xenobiotic/Detox ssGSEA
## (to integrate with existing OxStress / Immune traits)
## This script is a example to integrate traits even when 
## not strongly correlated but might be of interest
## we use Xenobiotic/Detox as example
###############################################

library(msigdbr)
library(GSVA)
library(WGCNA)
library(dplyr)
library(stringr)
library(ggplot2)

# ------------- 1. Build Xenobiotic / Detox gene set ----------------
# Get all human gene sets from MSigDB
msig_all <- msigdbr(species = "Homo sapiens")

# Keep Hallmark + KEGG + Reactome-type collections
msig_sub <- msig_all %>%
  dplyr::filter(
    gs_collection %in% c("H", "C2"),
    gs_subcollection %in% c(NA, "CP:KEGG", "CP:REACTOME")
  )

# Select sets related to xenobiotics / drug metabolism / CYP
xeno_raw <- msig_sub %>%
  dplyr::filter(
    str_detect(gs_name, regex("XENOBIOTIC|DRUG_METABOLISM|CYTOCHROME_P450|METABOLISM_OF_XENOBIOTICS",
               ignore_case = TRUE)) |
      str_detect(gs_description, regex("xenobiotic|detoxif",
                 ignore_case = TRUE))
  )

# (optional) inspect which sets were selected
unique(xeno_raw$gs_name)

# Keep only genes present in your expression matrix
xeno_mapped <- xeno_raw %>%
  dplyr::filter(gene_symbol %in% rownames(expression_matrix_human))

# Collapse all xenobiotic-related sets into a single union
xeno_union <- unique(xeno_mapped$gene_symbol)

length(xeno_union)   # should be >= 5

genesets_xeno <- list("Xenobiotic_Detox" = xeno_union)

# ------------- 2. Run ssGSEA ONLY for Xenobiotic_Detox ------------------
# gsva/ssgsea expects a matrix (genes x samples)
stopifnot(all(xeno_union %in% rownames(expression_matrix_human)))

ssgsea_param_xeno <- ssgseaParam(x,
  geneSets = genesets_xeno,
  minSize  = 5,
  maxSize  = 500
)

# if you use BiocParallel in the main script, keep the same BPPARAM
# example:
bp <- BiocParallel::SnowParam(workers = 4)

ssgsea_xeno <- gsva(ssgsea_param_xeno,
                    BPPARAM = bp,
                    verbose = TRUE)

# ssgsea_xeno: 1 x N_samples
# rowname = "Xenobiotic_Detox"
rownames(ssgsea_xeno)

# ------------- 3. Append Xenobiotic/Detox to your existing ssGSEA traits ------------------
# Matrix of ssGSEA categories (samples x categories)
ssgsea_traits <- t(ssgsea_by_category)
colnames(ssgsea_traits) <- paste0("ssGSEA_", colnames(ssgsea_traits))

# Make sure the names of samples match
common_samples <- intersect(rownames(MEs), rownames(ssgsea_traits))
common_samples <- intersect(common_samples, colnames(ssgsea_xeno))

# Subsets aligned
MEs_subset <- MEs[common_samples, , drop = FALSE]

ssgsea_traits_subset <- ssgsea_traits[common_samples, , drop = FALSE]

# Add the column Xenobiotic/Detox as numeric
ssgsea_xeno_vec <- as.numeric(ssgsea_xeno["Xenobiotic_Detox", common_samples])

ssgsea_traits_subset <- cbind(
  ssgsea_traits_subset,
  ssGSEA_XenobioticDetox = ssgsea_xeno_vec
)

# Convert to numerical matrix
ssgsea_traits_subset <- as.matrix(ssgsea_traits_subset)
storage.mode(ssgsea_traits_subset) <- "double"

cat("Dim MEs_subset:", dim(MEs_subset), "\n")
cat("Dim ssgsea_traits_subset:", dim(ssgsea_traits_subset), "\n")

# Save the matrix of ssGSEA by category (including Xenobiotic/Detox)
save(ssgsea_traits_subset,
     file = "./YOUR_DIRECTORY/ssgSEA_traits_by_category_with_Xeno.RData")

# ------------- 4. Recompute module-trait correlations including Xenobiotic/Detox ------------------
moduleTraitCor_ssGSEA <- cor(
  MEs_subset,
  ssgsea_traits_subset,
  use = "pairwise.complete.obs"
)

moduleTraitPvalue_ssGSEA <- WGCNA::corPvalueStudent(
  moduleTraitCor_ssGSEA,
  nSamples = nrow(MEs_subset)
)

# FDR across the full matrix
moduleTraitPvalue_adj_ssGSEA <- matrix(
  p.adjust(as.vector(moduleTraitPvalue_ssGSEA), method = "BH"),
  nrow     = nrow(moduleTraitPvalue_ssGSEA),
  dimnames = dimnames(moduleTraitPvalue_ssGSEA)
)

# Significance labels
signif_labels_ssGSEA <- ifelse(moduleTraitPvalue_adj_ssGSEA < 0.001, "***",
                               ifelse(moduleTraitPvalue_adj_ssGSEA < 0.01,  "**",
                                      ifelse(moduleTraitPvalue_adj_ssGSEA < 0.05,  "*",  "")))

## (Optional) heatmap WGCNA
pdf("./YOUR_DIRECTORY/Fig5_Module_ssGSEA_with_Xeno.pdf", width = 8, height = 8)
labeledHeatmap(
  Matrix     = moduleTraitCor_ssGSEA,
  xLabels    = colnames(ssgsea_traits_subset),
  yLabels    = colnames(MEs_subset),
  ySymbols   = colnames(MEs_subset),
  colorLabels = FALSE,
  colors     = blueWhiteRed(50),
  textMatrix = signif_labels_ssGSEA,
  setStdMargins = FALSE,
  cex.text   = 0.8,
  cex.lab.x  = 0.5,
  cex.lab.y  = 0.8,
  zlim       = c(-1, 1),
  main       = "Module–ssGSEA Category Correlation (incl. Xenobiotic/Detox)"
)
dev.off()

# Save data
write.csv(moduleTraitCor_ssGSEA, "./YOUR_DIRECTORY/Module_ssGSEA_correlation_Xeno.csv")
write.csv(moduleTraitPvalue_adj_ssGSEA, "./YOUR_DIRECTORY/Module_ssGSEA_pvalue_adj_Xeno.csv")

### Extract strong correlations ########
# strong correlation（|r| > 0.6, padj < 0.05）
strong_correlations <- which(
  abs(moduleTraitCor_ssGSEA) > 0.6 & moduleTraitPvalue_adj_ssGSEA < 0.05,
  arr.ind = TRUE
)

if (nrow(strong_correlations) > 0) {
  cat("\n=== Strong Module-ssGSEA correlations ===\n")
  strong_cor_df <- data.frame(
    Module = rownames(moduleTraitCor_ssGSEA)[strong_correlations[, 1]],
    ssGSEA_Category = colnames(moduleTraitCor_ssGSEA)[strong_correlations[, 2]],
    Correlation = moduleTraitCor_ssGSEA[strong_correlations],
    P_adj = moduleTraitPvalue_adj_ssGSEA[strong_correlations]
  ) %>%
    arrange(desc(abs(Correlation)))
  
  print(strong_cor_df)
  write.csv(strong_cor_df, "./YOUR_DIRECTORY/Strong_Module_ssGSEA_correlations_Xeno.csv", row.names = FALSE)
}

# ------------- 5. Example scatter for a module vs Xenobiotic/Detox ------------------
# choose any module that looks interesting in the heatmap
mod_of_interest <- "12"   # <-- change to "21", "13", etc.

if (mod_of_interest %in% colnames(MEs) &&
    "ssGSEA_XenobioticDetox" %in% colnames(ssgsea_traits_subset)) {
  
  plot_data_xeno <- data.frame(
    Module_eigengene = MEs_subset[, mod_of_interest],
    ssGSEA_score      = ssgsea_traits_subset[, "ssGSEA_XenobioticDetox"],
    Site  = meta_df$site[match(rownames(MEs_subset), meta_df$Sample)],
    Group = meta_df$condition[match(rownames(MEs_subset), meta_df$Sample)]
  )
  
  cor_test_xeno <- cor.test(plot_data_xeno$Module_eigengene,
                            plot_data_xeno$ssGSEA_score,
                            method = "spearman")
  
  fig_xeno <- ggplot(plot_data_xeno,
                     aes(x = Module_eigengene, y = ssGSEA_score,
                         color = Group)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE,
                color = "black", linetype = "dashed") +
    scale_color_manual(values = c("Farm" = "red", "Wild" = "lightblue")) +
    annotate("text", x = Inf, y = Inf,
             hjust = 1.1, vjust = 1.1,
             label = paste0("Spearman's rho = ",
                            round(cor_test_xeno$estimate, 3),
                            "\np = ",
                            format.pval(cor_test_xeno$p.value, digits = 2))) +
    labs(
      title = paste("Module", mod_of_interest, "vs Xenobiotic/Detox ssGSEA"),
      x     = paste("Module", mod_of_interest, "eigengene"),
      y     = "Xenobiotic/Detox ssGSEA score"
    ) +
    theme_bw(base_size = 12)
  
  ggsave(paste0("./YOUR_DIRECTORY/Fig5_Module",
                mod_of_interest, "_Xeno_ssGSEA_scatter.pdf"),
         fig_xeno, width = 7, height = 6)
}

#### -- EXTRA: Include the Xenobiotic/Detox module in the Boxplots -- ####
### Create the category-level data (incl. Xenobiotic/Detox) ###

# Load the matrix of the ssGSEA categories (samples x categories, with Xenobiotic)
load("./YOUR_DIRECTORY/ssgSEA_traits_by_category_with_Xeno.RData")
# object: ssgsea_traits_subset

ssgsea_category_long <- ssgsea_traits_subset %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  tidyr::pivot_longer(-Sample,
                      names_to  = "Category",
                      values_to = "Category_score") %>%
  dplyr::left_join(meta_df, by = "Sample") %>%
  dplyr::mutate(
    # Remove prefix and name in the facet
    Category = gsub("^ssGSEA_", "", Category),
    Category = dplyr::recode(
      Category,
      "OxidativeStress"  = "Oxidative stress",
      "Immune"           = "Immune / Inflammatory",
      "XenobioticDetox"  = "Xenobiotic / Detox"
    )
  )

####### Make boxplot (by category, 3 facets) ######
fig4a <- ggplot(ssgsea_category_long, 
                aes(x = site, y = Category_score, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  facet_wrap(~ Category, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Farm" = "red", "Wild" = "lightblue")) +
  labs(title = "Category-level ssGSEA scores by site",
       x = "Site",
       y = "ssGSEA score") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "gray90"),
        legend.position = "top")

ggsave("./YOUR_DIRECTORY/Fig4A_ssGSEA_boxplot_by_site_with_Xeno.pdf",
       fig4a, width = 12, height = 10)

########### Statistical test（Farm vs Wild, by category） #####################

stat_results <- ssgsea_category_long %>%
  group_by(Category) %>%
  wilcox_test(Category_score ~ condition) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

write.csv(stat_results,
          "./YOUR_DIRECTORY/ssGSEA_category_statistics_with_Xeno.csv",
          row.names = FALSE)

cat("\n=== Statistical tests (Farm vs Wild; incl. Xenobiotic/Detox) ===\n")
print(stat_results)
