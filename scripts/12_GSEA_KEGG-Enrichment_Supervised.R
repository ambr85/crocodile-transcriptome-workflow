#############################################################
#Programas:    GSEA & Cluster Profiler
#Description:  Biological Pathways affected by GSEA and KEGG
#Author:       Marisol Buenfil Rojas
#Date:         2025
#Notes:        Example with 1 Site. Supervised Analysis
#############################################################

library("ggplot2")
library("VennDiagram")
library("clusterProfiler")

install.packages("BiocManager")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")

library("AnnotationDbi")
library("org.Hs.eg.db")
library("AnnotationHub")

library(ggpubr)
library(dplyr)
library(zoo)
library(tidyr)

install.packages("e1071")
library(e1071)

##################
########## GSEA

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fgsea")

library(fgsea)
############################
directory <- "Your-path/GSEA"
setwd(directory)
getwd()

####Import Data

res_Site1 <- read.csv("DESeq2_Results_CC_vs_SL.csv", header = TRUE,
                         as.is=TRUE,row.names=1,  quote= "")
head(res_Site1)

res_Site2 <- read.csv("DESeq2_Results_CC_vs_CH.csv", header = TRUE,
                    as.is=TRUE,row.names=1,  quote= "")
head(res_Site2)

############################
# Use GSEA to assess if the top genes are enriched in biological pathways

## SITE 1 & 2
# Remove rows with NA in log2FoldChange for res_Site1
res_Site1_clean <- res_Site1[!is.na(res_Site1$X.log2FoldChange), ]
res_Site2_clean <- res_Site2[!is.na(res_Site2$X.log2FoldChange), ]

# Create a ranked gene list
ranked_gene_list_Site1 <- sort(res_Site1_clean$X.log2FoldChange, decreasing = TRUE)
ranked_gene_list_Site2 <- sort(res_Site2_clean$X.log2FoldChange, decreasing = TRUE)

# Clean up the gene symbols by removing the quotes
cleaned_gene_symbols_Site1 <- gsub("\"", "", rownames(res_Site1_clean))
cleaned_gene_symbols_Site2 <- gsub("\"", "", rownames(res_Site2_clean))

# Map cleaned gene symbols to Entrez IDs
entrez_ids_Site1 <- mapIds(org.Hs.eg.db, keys = cleaned_gene_symbols_Site1, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
entrez_ids_Site2 <- mapIds(org.Hs.eg.db, keys = cleaned_gene_symbols_Site2, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

# Remove NA values from the ranked gene list
ranked_gene_list_Site1_clean <- ranked_gene_list_Site1[!is.na(entrez_ids_Site1)]
names(ranked_gene_list_Site1_clean) <- entrez_ids_Site1[!is.na(entrez_ids_Site1)]

ranked_gene_list_Site2_clean <- ranked_gene_list_Site2[!is.na(entrez_ids_Site2)]
names(ranked_gene_list_Site2_clean) <- entrez_ids_Site2[!is.na(entrez_ids_Site2)]

# Check the length and the first few entries
length(ranked_gene_list_Site1_clean)
head(ranked_gene_list_Site1_clean)

length(ranked_gene_list_Site2_clean)
head(ranked_gene_list_Site2_clean)

###############
#####GSEA (2 SITES FACETED AND COMPARED, BUT INCLUDING ALL PATHWAYS) 
# Run GSEA for Site1
gsea_result_Site1 <- gseGO(geneList = ranked_gene_list_Site1_clean, 
                              OrgDb = org.Hs.eg.db, 
                              ont = "BP", 
                              pvalueCutoff = 0.05, 
                              minGSSize = 10, 
                              maxGSSize = 500)

# Run GSEA for Site2
gsea_result_Site2 <- gseGO(geneList = ranked_gene_list_Site2_clean, 
                           OrgDb = org.Hs.eg.db, 
                           ont = "BP", 
                           pvalueCutoff = 0.05, 
                           minGSSize = 10, 
                           maxGSSize = 500)

# Extracting results for res_Site1
gsea_data_Site1 <- gsea_result_Site1@result
gsea_data_Site1$dataset <- "SL"  # Add a column to label the dataset

# Extracting results for res_Site2
gsea_data_Site2 <- gsea_result_Site2@result
gsea_data_Site2$dataset <- "CH"  # Add a column to label the dataset

# Combine the result data frames into one
gsea_combined <- rbind(gsea_data_Site1, gsea_data_Site2)

# Check the combined data
head(gsea_combined)

# Sort by p-value and select the top 10 terms
gsea_combined_sorted <- gsea_combined[order(gsea_combined$pvalue), ]
top_gsea_terms <- head(gsea_combined_sorted, 10)

# Visualize top pathways for Site1
ggplot(gsea_data_Site1, aes(x = reorder(Description, NES), y = NES)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  xlab("Pathway Description") +
  ylab("Normalized Enrichment Score (NES)") +
  ggtitle("Top Pathways for Sabana Lake Dataset") +
  theme_minimal()

# Visualize top pathways for Site2
ggplot(gsea_data_Site2, aes(x = reorder(Description, NES), y = NES)) +
  geom_bar(stat = "identity", fill = "tomato") +
  coord_flip() +
  xlab("Pathway Description") +
  ylab("Normalized Enrichment Score (NES)") +
  ggtitle("Top Pathways for Chichankanaab Lake Dataset") +
  theme_minimal()

# Summarize top pathways for all Sites
head(gsea_data_Site1)
head(gsea_data_Site2)

# Create a combined plot using facets, but keep the datasets separate
gsea_combined_separate <- rbind(
  transform(gsea_data_Site1, dataset = "SL"),
  transform(gsea_data_Site2, dataset = "CH"),
)

# Plot with facets
ggplot(gsea_combined_separate, aes(x = reorder(Description, NES), y = NES)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~dataset) +
  xlab("Pathway Description") +
  ylab("Normalized Enrichment Score (NES)") +
  ggtitle("Top Pathways for Each Dataset") +
  theme_minimal()

##################
##FINAL CODE GSEA (SEPARATED TOP 10 PATHWAYS)
# Run GSEA for Site1
gsea_result_Site1 <- gseGO(geneList = ranked_gene_list_Site1_clean, 
                              OrgDb = org.Hs.eg.db, 
                              ont = "BP", 
                              pvalueCutoff = 0.05,     # Relaxed p-value 0.1 cutoff in case no pathways express
                              minGSSize = 5,          # Lower minimum gene set size
                              maxGSSize = 1000)       # Larger gene set size

# Run GSEA for Site2
gsea_result_Site2 <- gseGO(geneList = ranked_gene_list_Site2_clean, 
                           OrgDb = org.Hs.eg.db, 
                           ont = "BP", 
                           pvalueCutoff = 0.05,     # Relaxed p-value 0.1 cutoff in case no pathways express
                           minGSSize = 5,          # Lower minimum gene set size
                           maxGSSize = 1000)       # Larger gene set size

# Function to get top 10 pathways based on FDR-adjusted p-value
get_top_10_pathways <- function(gsea_result) {
  # Extract results from the gseaResult object (GSEA results table)
  gsea_table <- gsea_result@result
  
  # Make sure the table has a 'padj' column
  if ("padj" %in% colnames(gsea_table)) {
    # Sort by adjusted p-value (padj) in ascending order (best pathways at the top)
    sorted_gsea_table <- gsea_table[order(gsea_table$padj), ]
    
    # Take the top 10 pathways
    top_10_pathways <- head(sorted_gsea_table, 10)
  } else {
    # If there is no adjusted p-value (padj), fall back to p-value (pvalue)
    if ("pvalue" %in% colnames(gsea_table)) {
      sorted_gsea_table <- gsea_table[order(gsea_table$pvalue), ]
      top_10_pathways <- head(sorted_gsea_table, 10)
    } else {
      stop("No p-value or adjusted p-value column found in the GSEA results.")
    }
  }
  
  return(top_10_pathways)
}

# Get top 10 pathways for each dataset
top_10_Site1 <- get_top_10_pathways(gsea_result_Site1)
top_10_Site2 <- get_top_10_pathways(gsea_result_Site2)

# Check the top 10 pathways
print(top_10_Site1)
print(top_10_Site2)

# Load required libraries
library(reshape2)
library(pheatmap)

# Plot for Site1 (Top 10 significant pathways)
ggplot(top_10_Site1, aes(x = reorder(Description, NES), y = NES)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  xlab("Pathway Description") +
  ylab("Normalized Enrichment Score (NES)") +
  ggtitle("Top 10 Pathways for Sabana Lake") +
  theme_minimal()

# Plot for Site2 (Top 10 significant pathways)
ggplot(top_10_Site2, aes(x = reorder(Description, NES), y = NES)) +
  geom_bar(stat = "identity", fill = "tomato") +
  coord_flip() +
  xlab("Pathway Description") +
  ylab("Normalized Enrichment Score (NES)") +
  ggtitle("Top 10 Pathways for Chichankanaab Lake Dataset") +
  theme_minimal()

# Dotplot for Site1 (showing top 10 significant pathways)
dotplot(gsea_result_Site1, showCategory = 10) + 
  ggtitle("Top 10 Pathways for Sabana Lake Dataset")

# Dotplot for Site2 (showing top 10 significant pathways)
dotplot(gsea_result_Site2, showCategory = 10) + 
  ggtitle("Top 10 Pathways for Chichankanaab Lake Dataset")

#############################################
### KEGG enrichment pathways ###
## SUPERVISED ANALYSIS ##
#############################################

## It is necessary to use significant genes from REFvsSiteX. From that file, we need
## to extract from the annotation the Entrez ID values and log2foldchange values in a separate txt.
## Check the script of annotation after DESeq2. 

#in the case of needing to use KEGG from your own computer memory
library(devtools)

if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager::install("clusterProfiler")

#after everything is installed
library(enrichplot)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

directory <- "YourPath/KEGG_Enrichment"
setwd(directory)
getwd()

#Load Universe genes generated in the PC1 loadings step
load("./universe_genes.dat")

#################################
## Format for heatmaps 
# USE IT FOR ALL THE HEATMAPS
theme_heatmap_labels <- theme(
  axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4),
  axis.ticks.length.x = unit(1, "pt"),
  plot.margin  = margin(5, 14, 5, 5) 
)

############################
##Helpers for Sankeys and Macro-category plots
## ===== PATCH FOR SITES (KEGG classes + plots) =====
library(dplyr); library(tidyr); library(forcats)
library(ggplot2); library(ggalluvial); library(patchwork)
library(enrichplot); library(KEGGREST); library(stringr)

.kegg_cols_site <- c(
  "Metabolism"                           = "#F05A4A",
  "Genetic Information Processing"       = "#67C9D0",
  "Environmental Information Processing" = "#13A39A",
  "Cellular Processes"                   = "#2F4D8C",
  "Organismal Systems"                   = "#F4A48E",
  "Human Diseases"                       = "#7E89A4",
  "Other"                                = "grey70"
)

.safe_kegg_class_lookup <- function(ids_char) {
  # ids_char: vector type "hsa04110" o "path:hsa04110"
  ids <- ids_char[!is.na(ids_char)]
  if (!length(ids)) return(setNames(character(0), character(0)))
  ids_pref <- ifelse(grepl("^path:", ids), ids, paste0("path:", ids))
  out <- setNames(rep("Other", length(ids_pref)), sub("^path:", "", ids_pref))
  for (i in seq_along(ids_pref)) {
    res <- tryCatch(KEGGREST::keggGet(ids_pref[i]),
                    error = function(e) NULL,
                    warning = function(w) NULL)
    if (!is.null(res) && length(res) > 0) {
      cls <- res[[1]]$CLASS
      if (!is.null(cls) && length(cls) > 0) {
        top <- str_trim(strsplit(cls[1], ";")[[1]][1])
        if (!is.na(top) && nzchar(top)) {
          out[sub("^path:", "", ids_pref[i])] <- top
        }
      }
    }
    Sys.sleep(0.1)
  }
  out
}

add_kegg_class_site <- function(kk_readable_topN_df) {
  df <- kk_readable_topN_df
  if (!nrow(df)) return(df %>% mutate(kegg_class = "Other"))
  class_map <- .safe_kegg_class_lookup(df$ID)
  df$kegg_class <- unname(class_map[df$ID])
  df$kegg_class[is.na(df$kegg_class) | !(df$kegg_class %in% names(.kegg_cols_site))] <- "Other"
  df
}

kegg_bar_by_class_site <- function(kk_readable, show = 18, title = NULL){
  base_df <- as.data.frame(kk_readable) |>
    dplyr::arrange(p.adjust, dplyr::desc(Count))
  show_n <- min(show, nrow(base_df))
  if (show_n == 0) return(ggplot2::ggplot() + ggplot2::labs(title = paste0(title %||% "KEGG macro-categories", " (sin términos)")))
  df <- dplyr::slice_head(base_df, n = show_n) |>
    add_kegg_class_site() |>
    dplyr::mutate(
      Description = forcats::fct_reorder(Description, Count),
      kegg_class  = factor(kegg_class, levels = names(.kegg_cols_site))
    )
  
  ggplot2::ggplot(df, ggplot2::aes(x = Count, y = Description, fill = kegg_class)) +
    ggplot2::geom_col(width = .8) +
    ggplot2::geom_text(ggplot2::aes(label = Count), hjust = -0.25, size = 3) +
    ggplot2::scale_fill_manual(values = .kegg_cols_site, guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, .10))) +
    ggplot2::labs(x = "Gene Count", y = NULL, title = title) +
    ggplot2::facet_grid(kegg_class ~ ., scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold", size = 11, margin = ggplot2::margin(r=6)),
      axis.text.y = ggplot2::element_text(size = 9, lineheight = 0.9),
      panel.grid.major.y = ggplot2::element_blank()
    ) +
    ggplot2::scale_y_discrete(labels = \(s) stringr::str_wrap(s, width = 35))
}

kegg_sankey_with_dot_site <- function(kk_readable,
                                      fc_sym = NULL,
                                      show = 10,
                                      max_genes = 40,
                                      gene_label_size = 2.6,       
                                      widths = c(7, 1.2),          
                                      show_category_legend = TRUE, 
                                      title = NULL) {
  
  base_df <- as.data.frame(kk_readable) |>
    dplyr::arrange(p.adjust, dplyr::desc(Count))
  
  show_n <- min(show, nrow(base_df))
  if (show_n == 0) {
    if (is.null(title)) title <- "Genes ↔ KEGG pathways"
    return(ggplot2::ggplot() + ggplot2::labs(title = paste0(title, " (no terms)")))
  }
  
  base_df <- dplyr::slice_head(base_df, n = show_n)
  df <- add_kegg_class_site(base_df)
  
  # pairs gene ↔ pathway
  edges <- df |>
    dplyr::select(Description, kegg_class, geneID) |>
    tidyr::separate_rows(geneID, sep = "/") |>
    dplyr::rename(gene = geneID, term = Description)
  
  # add FC (if available) and limit # genes for legibility
  if (!is.null(fc_sym)) {
    fc_df <- tibble::tibble(gene = names(fc_sym), fc = as.numeric(fc_sym))
    edges <- dplyr::left_join(edges, fc_df, by = "gene")
  } else {
    edges$fc <- 0
  }
  
  if (is.finite(max_genes)) {
    keep_genes <- edges |>
      dplyr::group_by(gene) |>
      dplyr::summarise(maxabs = max(abs(fc), na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(dplyr::desc(maxabs)) |>
      dplyr::slice_head(n = max_genes) |>
      dplyr::pull(gene)
    edges <- dplyr::filter(edges, gene %in% keep_genes)
  }
  
  # order by genes |FC|
  gene_levels <- edges |>
    dplyr::group_by(gene) |>
    dplyr::summarise(maxabs = max(abs(fc), na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(maxabs)) |>
    dplyr::pull(gene)
  edges$gene <- factor(edges$gene, levels = rev(gene_levels))
  
  legend_guide <- if (show_category_legend) ggplot2::waiver() else "none"
  
  # Sankey
  p_sankey <- ggplot2::ggplot(edges, ggplot2::aes(axis1 = gene, axis2 = term, y = 1, fill = kegg_class)) +
    ggalluvial::geom_alluvium(width = 0.12, alpha = 0.9, size = 0.15, color = NA, show.legend = TRUE) +
    ggalluvial::geom_stratum(width = 0.12, fill = "grey95", color = "grey60", show.legend = FALSE) +
    ggplot2::geom_text(stat = "stratum", ggplot2::aes(label = ggplot2::after_stat(stratum)),
                       size = gene_label_size, check_overlap = TRUE) +
    ggplot2::scale_fill_manual(values = .kegg_cols_site, name = "KEGG category") +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 1))) +
    ggplot2::scale_x_discrete(limits = c("Genes", "Pathways"), expand = c(.05, .05)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(5,5,5,5))
  
  # Dotplot (small)
  show_dot <- min(show_n, nrow(as.data.frame(kk_readable)))
  p_dot <- enrichplot::dotplot(kk_readable, showCategory = show_dot) +
    ggplot2::labs(x = "Gene.Ratio") +
    ggplot2::theme_minimal(base_size = 11)
  
  (p_sankey | p_dot) +
    patchwork::plot_layout(widths = widths, guides = "collect") &  # merge legends
    ggplot2::theme(legend.position = "right")
  
  patch <- (p_sankey | p_dot) +
    patchwork::plot_layout(widths = widths, guides = "collect")
  
  patch +
    patchwork::plot_annotation(
      title = title,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 16, face = "bold"),
        plot.title.position = "plot"
      )
    ) &
    ggplot2::theme(legend.position = "right")
}

#################################
###Example: Site REF1 (repeat for all the sites)
mink<-read.table("SigGenes_REF1vsCC.txt", header= F)
datamink=mink[,2]
names(datamink)=as.character(mink[,1])
datamink = sort(datamink, decreasing =T)
genemink <-names(datamink)

KEGGmink <- enrichKEGG(gene         = genemink,
                       organism     = 'hsa',
                       universe     = universe_genes,
                       pvalueCutoff = 1,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 1)
head(KEGGmink)

# Convert to SYMBOL in the object
KEGGmink_sym <- setReadable(KEGGmink, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Rename the fold-change to SYMBOL (so it matches with the object)
entrez2sym <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                    keys = genemink,
                                    keytype = "ENTREZID",
                                    column = "SYMBOL",
                                    multiVals = "first")
keep <- !is.na(entrez2sym)
fc_sym <- datamink[!is.na(entrez2sym)]
names(fc_sym) <- entrez2sym[!is.na(entrez2sym)]

write.table (KEGGmink_sym, "KEGGSigGenes_REF1vsCC.txt", sep="\t")

heatplot(KEGGmink_sym, foldChange= fc_sym, showCategory = 20)+
  scale_fill_gradientn(name = "log2FC", 
                       colours = c("blue","white","red"), 
                       values=c(0, 0.5, 1), 
                       limits = c(-6, 7))+
  labs(title = "KEGG enrichment heatmap - REF1 vs CC", 
       subtitle = "Top terms; fold-change mapped to color")+
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_heatmap_labels

barplot(KEGGmink_sym, showCategory=20)+
  labs(title = "Top 20 KEGG pathways - REF1 vs CC")+
  labs(fill = "Adjusted p", color = "Adjusted p")

# GO REF1
GOmink <- enrichGO(genemink,
                   universe     = universe_genes,
                   OrgDb = org.Hs.eg.db, 
                   keyType = "ENTREZID",
                   ont = "BP", 
                   pvalueCutoff = 1,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 1,
                   readable = TRUE) # Returns gene in SYMBOL
head(GOmink)

write.table (GOmink, "GOSigGenes_REF1vsCC.txt", sep="\t")

heatplot(GOmink, foldChange= fc_sym, showCategory=20)+
  scale_fill_gradientn(colours = c("blue","white","red"), 
                       values=c(0, 0.5, 1), 
                       limits = c(-6, 7))+
  labs(title = "GO heatmap - REF1 vs CC", 
       subtitle = "Top terms; fold-change mapped to color")+
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_heatmap_labels +
  theme(axis.text.y = element_text(size = 6)) +                       
  scale_y_discrete(labels = \(s) stringr::str_wrap(s, width = 35))    

barplot(GOmink, showCategory=20)+
  labs(title = "Top 20 GO - REF1 vs CC")+
  labs(fill = "Adjusted p", color = "Adjusted p")

# EXTRAS REF1 
kegg_bar_by_class_site(KEGGmink_sym, show = 18,
                       title = "KEGG macro-categories – REF1 vs CC")

kegg_sankey_with_dot_site(
  KEGGmink_sym, fc_sym,
  show = 10,
  max_genes = 30,           # less genes = more readable labels
  gene_label_size = 3.2,    
  widths = c(8, 1),         # Sankey bigger than the dotplot
  show_category_legend = TRUE,
  title = "Genes ↔ KEGG pathways – REF1 vs CC"
)

##########################################
### KEGG and GO in Union and Intersect Genes

# Union (Repeat for Intersect by replacing the txt file (FinalIntersect_DEGs.txt). Make sure you change labels.
union <-read.table("FinalUnion_DEGs.txt", header= F)
dataunion=union[,1]
names(dataunion)=as.character(union[,1])
dataunion = sort(dataunion, decreasing =T)
geneunion <-names(dataunion)

# Run enrichKEGG
ekegg_union <- enrichKEGG(gene = geneunion,
                          universe     = universe_genes,
                          organism     = 'hsa',
                          pvalueCutoff = 1,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 1)

# Visualize
head(ekegg_union)

# Convert to SYMBOL in the object
KEGGmink_sym <- setReadable(ekegg_union, 
                            OrgDb = org.Hs.eg.db, 
                            keyType = "ENTREZID")

# Rename the FC to SYMBOL (to match with the object)
entrez2sym <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                    keys = geneunion,
                                    keytype = "ENTREZID",
                                    column = "SYMBOL",
                                    multiVals = "first")
keep   <- !is.na(entrez2sym)
fc_sym <- dataunion[!is.na(entrez2sym)]
names(fc_sym) <- entrez2sym[!is.na(entrez2sym)]

write.table (KEGGmink_sym, "KEGGSigGenes_Union.txt", sep="\t")

heatplot(KEGGmink_sym, foldChange= fc_sym, showCategory=20)+
  scale_fill_gradientn(colours = c("blue","white","red"), values=c(0, 0.5, 1), limits = c(-6, 7))+
  labs(title = "KEGG enrichment heatmap - Union Genes", subtitle = "Top terms; fold-change mapped to color")+ #No tiene caso por que no tenemos el dato log2fc
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_heatmap_labels

barplot(KEGGmink_sym, showCategory=20)+
  labs(title = "Top 20 KEGG pathways - Union Genes")+
  labs(fill = "Adjusted p", color = "Adjusted p")

# GO Union
GOunion <- enrichGO(geneunion,
                    universe     = universe_genes,
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID", 
                    ont = "BP", 
                    pvalueCutoff = 1,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 1,
                    readable = TRUE) 
head(GOunion)

write.table (GOunion, "GOSigGenes_UnionDEGs.txt", sep="\t")

heatplot(GOunion, foldChange= fc_sym, showCategory=20)+
  scale_fill_gradientn(colours = c("blue","white","red"), values=c(0, 0.5, 1), limits = c(-6, 7))+
  labs(title = "GO heatmap - Union Genes", subtitle = "Top terms; fold-change mapped to color")+ #For Union and Intersect there is no log2FC available
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_heatmap_labels +
  theme(axis.text.y = element_text(size = 6)) +                       
  scale_y_discrete(labels = \(s) stringr::str_wrap(s, width = 35))    

barplot(GOunion, showCategory=20)+
  labs(title = "Top 20 GO - Union Genes")+
  labs(fill = "Adjusted p", color = "Adjusted p")

#####################################################################
### EXTRA PLOTS: Macro-categories KEGG and Sankey/Dot (Union / Intersect)
#####################################################################
# Unión/Intersect
# Make sure you have this object (replace by ekegg_intersect for that analysis, and union word for intersect):
kegg_union_sym <- setReadable(ekegg_union, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Barplot by macro-category (top 18)
kegg_bar_by_class_site(kegg_union_sym, show = 18,
                       title = "KEGG macro-categories – Union DEGs")

# Sankey big + dotplot small (no FC)
kegg_sankey_with_dot_site(
  kegg_union_sym, fc_sym = NULL,
  show = 12,               # nº pathways in sankey and dot
  max_genes = 40,          # limit genes for readability
  gene_label_size = 3.2,   
  widths = c(8, 1),        # sankey:dot
  show_category_legend = TRUE,
  title = "Genes ↔ KEGG pathways – Union DEGs"
)
