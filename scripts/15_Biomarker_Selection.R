# ==============================================================================
# STEP6 Biomarker Candidate Selection: WGCNA Hub Genes + DESeq2 Integration
# ==============================================================================
# Purpose: Identify biomarker candidates from WGCNA modules with rigorous
#          statistical validation using DESeq2 differential expression results
#
# Target modules: darkred (21), gray60 (17), black (7), greenyellow (11) (Example)
#
# Pipeline:
# 1. Run DESeq2 for Farm vs Wild comparison
# 2. Calculate WGCNA module membership (kME)
# 3. Integrate DESeq2 + WGCNA results
# 4. Apply objective biomarker selection criteria
# 5. Evaluate classification performance (ROC-AUC)
# 6. Construct and validate multi-gene panel
# ==============================================================================

# Load required packages
library(DESeq2)        # Differential expression
library(WGCNA)         # Module membership
library(dplyr)
library(tidyr)
library(ggplot2)
library(pROC)          # ROC curve analysis
library(caret)         # Cross-validation
library(glmnet)        # Elastic net
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(ggrepel)

enableWGCNAThreads()

### ===== Step 1: Run DESeq2 for Farm vs Wild comparison =====
## ----------- 1. Prepare count matrix and metadata ------------
# Prepare count matrix and metadata
# Assuming you have:
# - raw_counts: raw count matrix (genes x samples), NOT normalized
# - metadata: data frame with columns: Sample, Group, Site
#   - Group: "Farm" or "Wild"
#   - Sample names must match column names in raw_counts

dds_full   <- readRDS("dds_full.rds")
raw_counts <- counts(dds_full)   # genes x samples, enteros ≥ 0

## (Optional) Remove suffix "countREV.out" to the columns (if any)
colnames(raw_counts) <- sub("countREV\\.out$", "", colnames(raw_counts))

metadata   <- meta_df
colnames(metadata)[colnames(metadata) == "condition"] <- "Group"
colnames(metadata)[colnames(metadata) == "site"]      <- "Site"
rownames(metadata) <- metadata$Sample

common_samples <- intersect(colnames(raw_counts), metadata$Sample)

cat("Samples counts:", length(colnames(raw_counts)), "\n")
cat("Samples in metadata:", nrow(metadata), "\n")
cat("Common samples:", length(common_samples), "\n")

raw_counts_sub <- raw_counts[, common_samples]
metadata_sub   <- metadata[common_samples, ]

# Sanity check
stopifnot(all(colnames(raw_counts_sub) == rownames(metadata_sub)))

# Verify data
cat("Checking input data...\n")
cat("Raw counts dimensions:", dim(raw_counts_sub), "\n")
cat("Metadata dimensions:", dim(metadata_sub), "\n")
cat("Sample name overlap:", 
    sum(colnames(raw_counts_sub) %in% metadata_sub$Sample), "/", 
    ncol(raw_counts_sub), "\n\n")

# Ensure metadata row names match count matrix column names
if (!all(colnames(raw_counts_sub) == metadata_sub$Sample)) {
  cat("Reordering metadata to match count matrix...\n")
  metadata_sub <- metadata_sub[match(colnames(raw_counts_sub), metadata_sub$Sample), ]
  rownames(metadata_sub) <- metadata_sub$Sample
}

# Ensure Group is a factor with correct reference level
metadata_sub$Group <- factor(metadata_sub$Group, levels = c("Farm", "Wild"))
cat("Group levels:", levels(metadata_sub$Group), "(reference: Farm)\n\n")

## ----------- 2. Create DESeq2 dataset ------------
cat("Creating DESeq2 dataset...\n")
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts_sub,
  colData = metadata_sub,
  design = ~ Group
)

cat("DESeq2 dataset created\n")
cat("  Samples:", ncol(dds), "\n")
cat("  Genes:", nrow(dds), "\n\n")

## ----------- 3. Pre-filtering: remove genes with very low counts ------------
# Keep genes with at least 10 counts in at least X samples
# X = minimum group size
min_samples <- min(table(metadata_sub$Group))
keep <- rowSums(counts(dds) >= 10) >= min_samples
dds <- dds[keep, ]

cat("After filtering low-count genes:\n")
cat("  Genes remaining:", nrow(dds), "\n\n")

## ----------- 4. Run DESeq2 ------------
cat("Running DESeq2 analysis...\n")
cat("This may take several minutes...\n\n")

dds <- DESeq(dds)

cat("DESeq2 analysis completed!\n\n")

# Extract results
deseq_results <- results(dds, 
                         contrast = c("Group", "Wild", "Farm"),
                         alpha = 0.05)  # Significance level

# Summary
cat("DESeq2 Results Summary:\n")
summary(deseq_results)
cat("\n")

## ----------- 5. Convert to data frame ------------
deseq_results_df <- as.data.frame(deseq_results)
deseq_results_df$Gene <- rownames(deseq_results_df)

# Count significant genes
sig_genes <- deseq_results_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

cat("Significant DEGs (padj < 0.05, |log2FC| > 1):", nrow(sig_genes), "\n")
cat("  Up in Wild:", sum(sig_genes$log2FoldChange > 0), "\n")
cat("  Down in Wild (Up in Farm):", sum(sig_genes$log2FoldChange < 0), "\n\n")

# Get normalized counts for downstream analysis
normalized_counts <- counts(dds, normalized = TRUE)

# Rlog transformation (for visualization and correlation)
rlog <- rlog(dds, blind = FALSE)
rlog_counts <- assay(rlog)

cat("Normalized counts extracted\n")
cat("rlog transformation completed\n\n")

## ----------- 6. Save DESeq2 results ------------
write.csv(deseq_results_df, "./YOUR_DIRECTORY/DESeq2_results_all_genes.csv", row.names = FALSE)
write.csv(sig_genes, "./YOUR_DIRECTORY/DESeq2_results_significant.csv", row.names = FALSE)

save(deseq_results_df, file = "./YOUR_DIRECTORY/deseq_results_df.dat")
save(sig_genes, file = "./YOUR_DIRECTORY/sig_genes.dat")
save(normalized_counts, file="./YOUR_DIRECTORY/normalized_counts.dat")
save(rlog_counts, file="./YOUR_DIRECTORY/rlog_counts.dat")
save(metadata_sub, file="./YOUR_DIRECTORY/metadata_sub.dat")

### ===== Step 2: WGCNA Module Membership Calculation =====
# Assuming you already have:
# - MEs: module eigengenes
# - moduleColors: module assignments (whether by colors or numbers)

## ----------- 7. Use transformed counts for WGCNA calculations ------------
# (consistent with how WGCNA network was built)
load("./YOUR_DIRECTORY/net.dat")
load("./YOUR_DIRECTORY/rlog_counts.dat")
load("./YOUR_DIRECTORY/moduleColors.dat")

datExpr <- t(rlog_counts)
MEs <- net$MEs

## ----------- 8. Calculating module membership (kME) ------------
cat("Calculating module membership (kME)...\n")
geneModuleMembership <- cor(datExpr, MEs, use = "p")
MMPvalue <- corPvalueStudent(geneModuleMembership, nrow(datExpr))

# Rename columns
colnames(geneModuleMembership) <- gsub("ME", "kME", colnames(geneModuleMembership))
colnames(MMPvalue) <- gsub("ME", "kME_pval", colnames(MMPvalue))

cat("Module membership calculated for", nrow(geneModuleMembership), "genes\n\n")

## ----------- 9. Save data ------------
save(geneModuleMembership, file="./YOUR_DIRECTORY/geneModuleMembership.dat")
save(MMPvalue, file="./YOUR_DIRECTORY/MMPvalue.dat")

### ===== Step 3: Define Selection Criteria =====
## ----------- 10. Target modules ------------
target_modules_colors <- c("darkred", "gray60", "black", "greenyellow", 
                    "red", "pink", "salmon", "turquoise", "blue", 
                    "brown", "yellow", "magenta", "tan", "cyan", "royalblue")

cat("Target modules:", paste(target_modules_colors, collapse = ", "), "\n\n")

# 2) Convert to numeric LABEL using LUT created before
target_modules <- unname(lut_color2label[target_modules_colors])
target_modules <- target_modules[!is.na(target_modules)]

cat("Target modules (labels):", paste(target_modules, collapse = ", "), "\n\n")

## ----------- 11. Biomarker selection criteria (objective, pre-defined) ------------
biomarker_criteria <- list(
  # Module membership
  kME_threshold = 0.6,              # Strong correlation with module eigengene (0.7 stringent)
  
  # Differential expression
  padj_threshold = 0.05,            # Stringent FDR threshold (0.01 is too stringent)
  log2FC_threshold = 1.0,           # Meaningful fold change (1.5 too stringent)
  log2FC_SE_ratio = 3,              # Effect size reliability (5 too stringent)
  
  # Expression level
  baseMean_threshold = 50,         # Sufficient expression for detection (100 too stringent)
  
  # Classification performance
  AUC_threshold = 0.75,             # Minimum discrimination ability 
  
  # Top N per module
  top_n_per_module = 20             # Maximum candidates per module
)

cat("Selection criteria:\n")
cat("  Module membership (kME) ≥", biomarker_criteria$kME_threshold, "\n")
cat("  Adjusted p-value <", biomarker_criteria$padj_threshold, "\n")
cat("  |log2 Fold Change| >", biomarker_criteria$log2FC_threshold, "\n")
cat("  |log2FC/lfcSE| >", biomarker_criteria$log2FC_SE_ratio, "\n")
cat("  Base mean expression >", biomarker_criteria$baseMean_threshold, "\n")
cat("  ROC-AUC ≥", biomarker_criteria$AUC_threshold, "\n\n")

## ----------- 12. Save data ------------
save.image("./YOUR_DIRECTORY/Criteria.RData")

### ===== Step 4: Integrate WGCNA + DESeq2 Results =====
## ----------- 13. Prepare the function to calculate AUC ------------
# Function to calculate AUC
calculate_auc <- function(gene_name, expr_matrix, labels) {
  expr_values <- expr_matrix[gene_name, ]
  common_samples <- intersect(names(expr_values), names(labels))
  expr_values <- expr_values[common_samples]
  labels <- labels[common_samples]
  
  roc_obj <- roc(labels, expr_values, direction = "auto", quiet = TRUE)
  auc_value <- as.numeric(auc(roc_obj))
  ci_obj <- ci.auc(roc_obj, conf.level = 0.95)
  
  return(list(
    AUC = auc_value,
    AUC_lower = ci_obj[1],
    AUC_upper = ci_obj[3]
  ))
}

## ----------- 14. Load data ------------
load(file="./YOUR_DIRECTORY/moduleColors.dat")
load(file="./YOUR_DIRECTORY/metadata_sub.dat")
load(file="./YOUR_DIRECTORY/Criteria.RData")
load(file="./YOUR_DIRECTORY/geneModuleMembership.dat")
load(file="./YOUR_DIRECTORY/deseq_results_df.dat")
load(file="./YOUR_DIRECTORY/rlog_counts.dat")

## ----------- 15. Group labels for AUC calculation ------------
group_labels <- as.numeric(metadata_sub$Group) - 1  # Farm=0, Wild=1
names(group_labels) <- metadata_sub$Sample

## ----------- 16. Process each target module ------------
biomarker_candidates_list <- list()
module_summaries <- list()

common_rows_kME_MMP <- intersect(rownames(geneModuleMembership),
                                 rownames(MMPvalue))
cat("Common genes between geneModuleMembership and MMPvalue:",
    length(common_rows_kME_MMP), "\n")

for (module in target_modules) {

  # Get genes in this module
  module_num <- as.numeric(module)
  module_genes_all <- names(moduleLabels)[moduleLabels == module_num]

  module_genes <- intersect(module_genes_all, common_rows_kME_MMP)
  n_genes <- length(module_genes)
  
  cat("Module", module, "- total genes in labels:", length(module_genes_all),
      "| present in kME+MMPvalue:", n_genes, "\n")
  
  if (n_genes == 0) {
    warning("No genes found in module ", module, " after intersecting with kME/MMPvalue rows")
    next
  }
  
  # Get kME for this module
  kME_column <- paste0("kME", module)
  
  if (!kME_column %in% colnames(geneModuleMembership)) {
    warning("kME column not found for module ", module, " (", kME_column, ")" )
    next
  }
  
  # Create integrated data frame
  module_data <- data.frame(
    Gene = module_genes,
    Module = module,
    stringsAsFactors = FALSE
  )
  
  # Add WGCNA metrics
  module_data$kME <- geneModuleMembership[module_genes, kME_column]
  # p-values use "kME_pval<module>"
  pval_column <- paste0("kME_pval", module)
  
  if (!pval_column %in% colnames(MMPvalue)) {
    warning("p-value column not found for module ", module, " (", pval_column, ")")
    module_data$kME_pval <- NA
  } else {
    module_data$kME_pval <- MMPvalue[module_genes, pval_column]
  }

  # Add DESeq2 results
  deseq_subset <- deseq_results_df[match(module_genes, deseq_results_df$Gene), ]
  
  module_data$baseMean <- deseq_subset$baseMean
  module_data$log2FC <- deseq_subset$log2FoldChange
  module_data$lfcSE <- deseq_subset$lfcSE
  module_data$stat <- deseq_subset$stat
  module_data$pvalue <- deseq_subset$pvalue
  module_data$padj <- deseq_subset$padj
  
  # Calculate effect size ratio
  module_data$log2FC_SE_ratio <- abs(module_data$log2FC / module_data$lfcSE)
  
  # Remove genes with NA values (failed in DESeq2)
  module_data <- module_data %>% filter(!is.na(padj))
  
  cat("Genes with valid DESeq2 results:", nrow(module_data), "\n\n")
  
  # Apply selection criteria (step by step for transparency)
  cat("Applying selection criteria:\n")
  
  # Criterion 1: Module membership
  pass_kME <- module_data %>% filter(kME >= biomarker_criteria$kME_threshold)
  cat("  Pass kME ≥", biomarker_criteria$kME_threshold, ":", nrow(pass_kME), "\n")
  
  # Criterion 2: Statistical significance
  pass_sig <- pass_kME %>% 
    filter(padj < biomarker_criteria$padj_threshold,
           abs(log2FC) > biomarker_criteria$log2FC_threshold)
  cat("  Pass padj <", biomarker_criteria$padj_threshold, 
      "& |log2FC| >", biomarker_criteria$log2FC_threshold, ":", nrow(pass_sig), "\n")
  
  # Criterion 3: Expression level
  pass_expr <- pass_sig %>% 
    filter(baseMean > biomarker_criteria$baseMean_threshold)
  cat("  Pass baseMean >", biomarker_criteria$baseMean_threshold, ":", 
      nrow(pass_expr), "\n")
  
  # Criterion 4: Effect size reliability
  pass_effect <- pass_expr %>% 
    filter(log2FC_SE_ratio > biomarker_criteria$log2FC_SE_ratio)
  cat("  Pass |log2FC/lfcSE| >", biomarker_criteria$log2FC_SE_ratio, ":", 
      nrow(pass_effect), "\n\n")
  
  if (nrow(pass_effect) == 0) {
    cat("  No candidates pass all criteria for module", module, "\n")
    next
  }
  
  # Calculate AUC for remaining candidates
  cat("Calculating ROC-AUC for", nrow(pass_effect), "candidates...\n")
  
  auc_results <- lapply(pass_effect$Gene, function(gene) {
    result <- calculate_auc(gene, rlog_counts, group_labels)
    data.frame(
      Gene = gene,
      AUC = result$AUC,
      AUC_lower = result$AUC_lower,
      AUC_upper = result$AUC_upper
    )
  })
  
  auc_df <- bind_rows(auc_results)
  
  # Merge with candidate data
  candidates <- pass_effect %>%
    left_join(auc_df, by = "Gene")
  
  # Criterion 5: Classification performance
  high_performance <- candidates %>%
    filter(AUC >= biomarker_criteria$AUC_threshold) %>%
    arrange(desc(AUC))
  
  cat("  Pass AUC ≥", biomarker_criteria$AUC_threshold, ":", 
      nrow(high_performance), "\n\n")
  
  # Limit to top N per module
  if (nrow(high_performance) > biomarker_criteria$top_n_per_module) {
    high_performance <- high_performance %>%
      slice_head(n = biomarker_criteria$top_n_per_module)
    cat("  Limited to top", biomarker_criteria$top_n_per_module, "by AUC\n\n")
  }
  
  # Store results
  biomarker_candidates_list[[module]] <- high_performance
  
  # Module summary
  module_summary <- data.frame(
    Module = module,
    Total_genes = n_genes,
    Pass_kME = nrow(pass_kME),
    Pass_DE = nrow(pass_sig),
    Pass_expression = nrow(pass_expr),
    Pass_effect_size = nrow(pass_effect),
    Final_candidates = nrow(high_performance),
    Mean_kME = mean(high_performance$kME),
    Mean_log2FC = mean(abs(high_performance$log2FC)),
    Mean_AUC = mean(high_performance$AUC)
  )
  
  module_summaries[[module]] <- module_summary
  
  # Print top candidates
  cat("Top 5 candidates for module", module, ":\n")
  print(high_performance %>% 
          select(Gene, kME, log2FC, padj, baseMean, AUC) %>%
          head(5))
  cat("\n")
}

## ----------- 17. Combine results ------------
all_candidates <- bind_rows(biomarker_candidates_list)
module_summary_table <- bind_rows(module_summaries)

## ----------- 18. Save data ------------
save(all_candidates, file="./YOUR_DIRECTORY/all_candidates.dat")
save(module_summary_table, file="./YOUR_DIRECTORY/module_summary_table.dat")
save(group_labels, file="./YOUR_DIRECTORY/group_labels.dat")

write.csv(module_summary_table,
          "./YOUR_DIRECTORY/Module_summary.csv",
          row.names = FALSE)

### ===== Step 5: Prioritize Candidates =====
## ----------- 19. Load data ------------
load(file="./YOUR_DIRECTORY/all_candidates.dat")

## ----------- 20. Create composite score for ranking ------------
all_candidates <- all_candidates %>%
  mutate(
    # Normalize metrics to 0-1 scale
    kME_norm = (kME - min(kME)) / (max(kME) - min(kME)),
    AUC_norm = (AUC - min(AUC)) / (max(AUC) - min(AUC)),
    log2FC_norm = (abs(log2FC) - min(abs(log2FC))) / 
      (max(abs(log2FC)) - min(abs(log2FC))),
    padj_norm = 1 - (log10(padj) / log10(max(padj))),  # Lower padj = higher score
    
    # Composite score (weighted average)
    # AUC is most important for biomarkers
    Composite_score = 0.40 * AUC_norm + 
      0.25 * kME_norm + 
      0.20 * log2FC_norm +
      0.15 * padj_norm,
    
    # Priority classification
    Priority = case_when(
      AUC >= 0.90 & kME >= 0.85 & abs(log2FC) >= 2.0 ~ "High",
      AUC >= 0.80 & kME >= 0.75 & abs(log2FC) >= 1.5 ~ "Medium",
      TRUE ~ "Low"
    )
  ) %>%
  arrange(desc(Composite_score))

cat("Priority Distribution:\n")
print(table(all_candidates$Priority))
cat("\n")

cat("Top 10 biomarker candidates overall:\n")
print(all_candidates %>%
        select(Gene, Module, kME, log2FC, padj, AUC, Priority, Composite_score) %>%
        head(10))
cat("\n")

high_priority <- all_candidates %>% filter(Priority == "High")

## ----------- 20A. Panel for multigene biomarker (less strict, but still strong) ------------
# Panel filter: still good markers, but not as extreme as "High"
panel_candidates <- all_candidates %>%
  dplyr::filter(
    AUC   >= 0.80,          # keeps good classifiers
    kME   >= 0.60,          # reasonably central in the module
    abs(log2FC) >= 1.2      # meaningful effect size
  )

# From those, choose top genes per module by composite score
panel_genes2 <- panel_candidates %>%
  dplyr::group_by(Module) %>%
  dplyr::slice_max(order_by = Composite_score, n = 3, with_ties = FALSE) %>%
  dplyr::ungroup()

cat("Number of genes in panel (all modules):", nrow(panel_genes2), "\n")
print(panel_genes2 %>%
        dplyr::select(Gene, Module, kME, log2FC, padj, AUC, Composite_score))

## ----------- 21. Save data ------------
save(all_candidates, file="./YOUR_DIRECTORY/all_candidates_ranked.dat")
save(panel_genes2,
     file = "./YOUR_DIRECTORY/Biomarker_panel_multigene.dat")

write.csv(all_candidates,
          "./YOUR_DIRECTORY/Biomarker_candidates_all.csv",
          row.names = FALSE)
write.csv(high_priority,
          "./YOUR_DIRECTORY/Biomarker_candidates_high_priority.csv",
          row.names = FALSE)
write.csv(panel_genes2,
          "./YOUR_DIRECTORY/Biomarker_panel_multigene.csv",
          row.names = FALSE)

### ===== Step 6: Construct Multi-gene Panel =====
set.seed(123)

## ----------- 22. Load data ------------
load("./YOUR_DIRECTORY/all_candidates_ranked.dat")
load("./YOUR_DIRECTORY/rlog_counts.dat")

## ----------- 23. Select top candidates for panel ------------
# Aim for 10-15 genes from high and medium priority
panel_genes <- all_candidates %>%
  filter(Priority %in% c("High", "Medium")) %>%
  arrange(desc(Composite_score)) %>%
  slice_head(n = 10) %>%                  #change the number 
  pull(Gene)

cat("Panel genes selected:", length(panel_genes), "\n")
cat("Genes:\n")
for (i in seq_along(panel_genes)) {
  gene_info <- all_candidates %>% filter(Gene == panel_genes[i])
  cat(sprintf("  %d. %s (Module: %s, AUC: %.3f, log2FC: %.2f)\n",
              i, gene_info$Gene, gene_info$Module, gene_info$AUC, gene_info$log2FC))
}
cat("\n")

## ----------- 24. Panel strategy 1: Simple mean expression ------------
cat("Strategy 1: Mean expression\n")
panel_expr <- rlog_counts[panel_genes, ]
panel_score_mean <- colMeans(panel_expr)

panel_roc_mean <- roc(group_labels, panel_score_mean, quiet = TRUE)
panel_auc_mean <- auc(panel_roc_mean)
panel_ci_mean <- ci.auc(panel_roc_mean)

cat("  AUC:", round(panel_auc_mean, 3), "\n")
cat("  95% CI:", round(panel_ci_mean[1], 3), "-", round(panel_ci_mean[3], 3), "\n\n")

## ----------- 25. Panel Strategy 2: Weighted by individual AUC ------------
cat("Strategy 2: AUC-weighted\n")
gene_weights <- all_candidates %>%
  filter(Gene %in% panel_genes) %>%
  arrange(match(Gene, panel_genes)) %>%
  pull(AUC)

gene_weights <- gene_weights / sum(gene_weights)

panel_score_weighted <- apply(panel_expr, 2, function(x) sum(x * gene_weights))

panel_roc_weighted <- roc(group_labels, panel_score_weighted, quiet = TRUE)
panel_auc_weighted <- auc(panel_roc_weighted)
panel_ci_weighted <- ci.auc(panel_roc_weighted)

cat("  AUC:", round(panel_auc_weighted, 3), "\n")
cat("  95% CI:", round(panel_ci_weighted[1], 3), "-", round(panel_ci_weighted[3], 3), "\n\n")

## ----------- 26. Panel Strategy 3: Elastic Net ------------
cat("Strategy 3: Elastic Net (penalized regression)\n")

X_panel <- t(panel_expr)
y_panel <- group_labels

# Cross-validation to find optimal lambda
cv_glmnet <- cv.glmnet(X_panel, y_panel,
                       family = "binomial",
                       alpha = 0.5,  # Elastic net
                       nfolds = 10,
                       type.measure = "auc")

cat("  Optimal lambda:", cv_glmnet$lambda.min, "\n")

# Predict with optimal model
panel_score_glmnet <- predict(cv_glmnet, 
                              newx = X_panel,
                              s = "lambda.min",
                              type = "response")

panel_roc_glmnet <- roc(y_panel, as.vector(panel_score_glmnet), quiet = TRUE)
panel_auc_glmnet <- auc(panel_roc_glmnet)
panel_ci_glmnet <- ci.auc(panel_roc_glmnet)

cat("  AUC:", round(panel_auc_glmnet, 3), "\n")
cat("  95% CI:", round(panel_ci_glmnet[1], 3), "-", round(panel_ci_glmnet[3], 3), "\n\n")

# Get coefficients
panel_coef <- coef(cv_glmnet, s = "lambda.min")
panel_coef_df <- data.frame(
  Gene = rownames(panel_coef)[-1],  # Exclude intercept
  Coefficient = as.vector(panel_coef)[-1]
) %>%
  filter(Coefficient != 0) %>%
  arrange(desc(abs(Coefficient)))

cat("  Non-zero coefficients:", nrow(panel_coef_df), "/", length(panel_genes), "\n")
cat("  Active genes in panel:\n")
print(panel_coef_df)
cat("\n")

## ----------- 27. Panel Summary ------------
panel_summary <- all_candidates %>%
  filter(Gene %in% panel_genes) %>%
  left_join(panel_coef_df, by = "Gene") %>%
  mutate(Weight_AUC = gene_weights) %>%
  select(Gene, Module, kME, log2FC, padj, baseMean, AUC, AUC_lower, AUC_upper,
         Composite_score, Priority, Weight_AUC, Coefficient) %>%
  arrange(desc(AUC))

# Export 4: Panel scores for all samples
panel_scores_df <- data.frame(
  Sample = colnames(rlog_counts),
  Score_mean = panel_score_mean,
  Score_weighted = panel_score_weighted,
  Score_glmnet = as.vector(panel_score_glmnet),
  Group = metadata$Group[match(colnames(rlog_counts), metadata$Sample)],
  Site = metadata$Site[match(colnames(rlog_counts), metadata$Sample)]
)

## ----------- 28. Save data ------------
save(y_panel, file = "./YOUR_DIRECTORY/y_panel.dat")
save(X_panel, file = "./YOUR_DIRECTORY/X_panel.dat")

write.csv(panel_summary,
          "./YOUR_DIRECTORY/Biomarker_panel_summary.csv",
          row.names = FALSE)
write.csv(panel_scores_df,
          "./YOUR_DIRECTORY/Panel_scores_all_samples.csv",
          row.names = FALSE)

### ===== Step 7: Cross-validation =====
## ----------- 29. Load data ------------
load(file="./YOUR_DIRECTORY/y_panel.dat")
load(file="./YOUR_DIRECTORY/X_panel.dat")

## ----------- 30. Performing Leave-One-Out Cross-Validation (LOOCV) ------------
loocv_predictions <- numeric(length(y_panel))

for (i in seq_along(y_panel)) {
  # Leave one out
  X_train <- X_panel[-i, ]
  y_train <- y_panel[-i]
  X_test <- X_panel[i, , drop = FALSE]
  
  # Train model
  cv_temp <- cv.glmnet(X_train, y_train,
                       family = "binomial",
                       alpha = 0.5,
                       nfolds = min(10, length(y_train)))
  
  # Predict
  loocv_predictions[i] <- predict(cv_temp,
                                  newx = X_test,
                                  s = "lambda.min",
                                  type = "response")
}

# Calculate CV AUC
loocv_roc <- roc(y_panel, loocv_predictions, quiet = TRUE)
loocv_auc <- auc(loocv_roc)

cat("Cross-validated AUC (LOOCV):", round(loocv_auc, 3), "\n\n")

## ----------- 31. Save data ------------
save(loocv_roc, file = "./YOUR_DIRECTORY/loocv_roc.dat")
save(loocv_auc, file = "./YOUR_DIRECTORY/loocv_auc.dat")

### ===== Step 8: Visualizations (optional) =====
## ----------- 32. Load data ------------
load(file="./YOUR_DIRECTORY/module_summary_table.dat")
load(file="./YOUR_DIRECTORY/module_summary_table.dat")
load(file="./YOUR_DIRECTORY/deseq_results_df.dat") 
load(file="./YOUR_DIRECTORY/all_candidates.dat") 
load(file="./YOUR_DIRECTORY/rlog_counts.dat")
load(file="./YOUR_DIRECTORY/group_labels.dat")
load(file="./YOUR_DIRECTORY/panel_scores_df.dat")
load(file="./YOUR_DIRECTORY/panel_expr.dat")
load(file="./YOUR_DIRECTORY/metadata.dat")

## ----------- 33. Figure A: Selection funnel ------------
funnel_data <- module_summary_table %>%
  select(Module, Total_genes, Pass_kME, Pass_DE, Pass_expression, 
         Pass_effect_size, Final_candidates) %>%
  pivot_longer(-Module, names_to = "Stage", values_to = "N_genes") %>%
  mutate(Stage = factor(Stage, levels = c("Total_genes", "Pass_kME", "Pass_DE",
                                          "Pass_expression", "Pass_effect_size",
                                          "Final_candidates")))

figA <- ggplot(funnel_data, aes(x = Stage, y = N_genes, fill = Module, group = Module)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = setNames(target_modules, target_modules)) +
  labs(title = "Biomarker Candidate Selection Funnel",
       subtitle = "Progressive filtering by objective criteria",
       x = "Selection Stage",
       y = "Number of Genes") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## ----------- 34. Figure B: Volcano plot with candidates highlighted ------------
volcano_data <- deseq_results_df %>%
  mutate(
    Significance = case_when(
      Gene %in% all_candidates$Gene ~ "Biomarker candidate",
      padj < 0.05 & abs(log2FoldChange) > 1 ~ "Significant DEG",
      TRUE ~ "Not significant"
    )
  )

figB <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(data = filter(volcano_data, Gene %in% panel_genes),
                  aes(label = Gene),
                  size = 3,
                  max.overlaps = 15,
                  color = "black") +
  geom_vline(xintercept = c(-biomarker_criteria$log2FC_threshold, 
                            biomarker_criteria$log2FC_threshold),
             linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(biomarker_criteria$padj_threshold),
             linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Biomarker candidate" = "red",
                                "Significant DEG" = "blue",
                                "Not significant" = "gray80")) +
  labs(title = "Differential Expression: Farm vs Wild",
       subtitle = "Biomarker candidates highlighted",
       x = "log2 Fold Change (Wild / Farm)",
       y = "-log10(adjusted p-value)") +
  theme_bw(base_size = 12)

ggsave("./YOUR_DIRECTORY/Fig6B_Volcano_plot.pdf", figB, width = 10, height = 8)

## ----------- 35. Figure C: AUC vs log2FC scatter ------------
figC <- ggplot(all_candidates, aes(x = abs(log2FC), y = AUC, 
                                   color = Priority, size = kME)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(data = filter(all_candidates, Priority == "High"),
                  aes(label = Gene),
                  size = 3,
                  max.overlaps = 15) +
  scale_color_manual(values = c("High" = "red", "Medium" = "orange", "Low" = "gray70")) +
  scale_size_continuous(range = c(2, 8)) +
  labs(title = "Biomarker Performance vs Effect Size",
       x = "|log2 Fold Change|",
       y = "ROC-AUC",
       color = "Priority",
       size = "Module\nMembership") +
  theme_bw(base_size = 12)

ggsave("./YOUR_DIRECTORY/Fig6C_AUC_vs_log2FC.pdf", figC, width = 10, height = 7)

## ----------- 36. Figure D: ROC curves ------------
figD <- ggplot() +
  # Top 5 individual genes
  lapply(1:min(5, nrow(all_candidates)), function(i) {
    gene <- all_candidates$Gene[i]
    gene_expr <- rlog_counts[gene, ]
    roc_obj <- roc(group_labels, gene_expr, quiet = TRUE)
    
    roc_data <- data.frame(
      Specificity = roc_obj$specificities,
      Sensitivity = roc_obj$sensitivities,
      Gene = paste0(gene, " (", round(auc(roc_obj), 2), ")")
    )
    
    geom_line(data = roc_data,
              aes(x = 1 - Specificity, y = Sensitivity, color = Gene),
              linewidth = 0.8, alpha = 0.7)
  }) +
  
  # Panel ROC
  geom_line(data = data.frame(
    Specificity = panel_roc_glmnet$specificities,
    Sensitivity = panel_roc_glmnet$sensitivities
  ),
  aes(x = 1 - Specificity, y = Sensitivity),
  color = "black", linewidth = 1.5) +
  
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  
  annotate("text", x = 0.7, y = 0.3,
           label = paste0("Panel (n=", length(panel_genes), ")\nAUC: ",
                          round(panel_auc_glmnet, 3)),
           size = 4, fontface = "bold") +
  
  scale_color_brewer(palette = "Set1") +
  labs(title = "ROC Curves: Individual Genes and Multi-gene Panel",
       x = "1 - Specificity",
       y = "Sensitivity",
       color = "Top Genes") +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.7, 0.3))

ggsave("./YOUR_DIRECTORY/Fig6D_ROC_curves.pdf", figD, width = 8, height = 8)

## ----------- 37. Figure E: Panel score distribution ------------
figE <- ggplot(panel_scores_df, aes(x = Group, y = Score_glmnet, fill = Group)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 2) +
  scale_fill_manual(values = c("Farm" = "red", "Wild" = "lightblue")) +
  stat_compare_means(method = "wilcox.test", size = 5) +
  labs(title = "Multi-gene Panel Score Distribution",
       subtitle = paste0("Panel: ", length(panel_genes), " genes | ",
                         "AUC: ", round(panel_auc_glmnet, 3), " | ",
                         "CV-AUC: ", round(loocv_auc, 3)),
       x = "Group",
       y = "Panel Score (Elastic Net)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave("./YOUR_DIRECTORY/Fig6E_Panel_distribution.pdf", figE, width = 7, height = 6)

## ----------- 38. Figure F: Heatmap of panel genes ------------
panel_expr_zscore <- t(scale(t(panel_expr)))
sample_order <- order(panel_scores_df$Score_glmnet)
panel_expr_zscore <- panel_expr_zscore[, sample_order]

col_anno <- HeatmapAnnotation(
  Group = metadata$Group[match(colnames(panel_expr_zscore), metadata$Sample)],
  Site = metadata$Site[match(colnames(panel_expr_zscore), metadata$Sample)],
  Panel_Score = anno_barplot(panel_scores_df$Score_glmnet[sample_order],
                             gp = gpar(fill = "darkgreen"),
                             height = unit(2, "cm")),
  col = list(
    Group = c("Farm" = "red", "Wild" = "lightblue"),
    Site = setNames(
      colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(unique(metadata$Site))),
      unique(metadata$Site)
    )
  )
)

row_anno <- rowAnnotation(
  Module = panel_summary$Module[match(rownames(panel_expr_zscore), panel_summary$Gene)],
  AUC = anno_barplot(panel_summary$AUC[match(rownames(panel_expr_zscore), panel_summary$Gene)],
                     gp = gpar(fill = "steelblue"),
                     width = unit(2, "cm")),
  col = list(Module = setNames(target_modules, target_modules))
)

figF <- Heatmap(
  panel_expr_zscore,
  name = "Z-score",
  top_annotation = col_anno,
  right_annotation = row_anno,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 10),
  column_title = "Samples (ordered by Panel Score)",
  row_title = "Panel Genes"
)

pdf("./YOUR_DIRECTORY/Fig6F_Panel_heatmap.pdf", width = 14, height = 8)
draw(figF)
dev.off()

### ===== Step 9: Create Manuscript Tables (optional) =====
## ----------- 39. Table 1: Top biomarker candidates ------------
table1 <- all_candidates %>%
  filter(Priority %in% c("High", "Medium")) %>%
  select(Gene, Module, log2FC, padj, baseMean, kME, AUC, AUC_lower, AUC_upper, Priority) %>%
  mutate(
    log2FC = round(log2FC, 2),
    padj = format(padj, scientific = TRUE, digits = 2),
    baseMean = round(baseMean, 0),
    kME = round(kME, 3),
    AUC = round(AUC, 3),
    CI_95 = paste0("[", round(AUC_lower, 3), "-", round(AUC_upper, 3), "]")
  ) %>%
  select(-AUC_lower, -AUC_upper) %>%
  arrange(desc(AUC))

write.csv(table1, "./YOUR_DIRECTORY/Table1_Top_candidates.csv", row.names = FALSE)

## ----------- 40. Table 2: Panel performance ------------
table2 <- data.frame(
  Strategy = c("Mean expression", "AUC-weighted", "Elastic Net", "Elastic Net (LOOCV)"),
  N_genes = rep(length(panel_genes), 4),
  AUC = c(panel_auc_mean, panel_auc_weighted, panel_auc_glmnet, loocv_auc),
  CI_lower = c(panel_ci_mean[1], panel_ci_weighted[1], panel_ci_glmnet[1], NA),
  CI_upper = c(panel_ci_mean[3], panel_ci_weighted[3], panel_ci_glmnet[3], NA)
) %>%
  mutate(
    AUC = round(AUC, 3),
    CI_95 = ifelse(is.na(CI_lower), "N/A",
                   paste0("[", round(CI_lower, 3), "-", round(CI_upper, 3), "]"))
  ) %>%
  select(Strategy, N_genes, AUC, CI_95)

write.csv(table2, "./YOUR_DIRECTORY/Table2_Panel_performance.csv", row.names = FALSE)

