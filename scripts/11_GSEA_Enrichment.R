#############################################################
#Programas:    GSEA & Cluster Profiler
#Description:  Biological Pathways affected by GSEA and KEGG
#Author:       Marisol Buenfil Rojas
#Date:         2025
#Notes:        Example with 2 Sites
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

#################################

###Site 1

# KEGG Site 1

mink<-read.table("SigGenes_REFvsSite1.txt", header= F)
datamink=mink[,2]
names(datamink)=as.character(mink[,1])
datamink = sort(datamink, decreasing =T)
genemink <-names(datamink)

KEGGmink <- enrichKEGG(gene         = genemink,
                       organism     = 'hsa', 
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "fdr",
                       qvalueCutoff = 0.05)
head(KEGGmink)

write.table (KEGGmink, "KEGGSigGenes_REFvsSite1.txt", sep="\t")

heatplot(KEGGmink, foldChange=datamink)+
  scale_fill_gradientn(colours = c("blue","white","red"), values=c(0, 0.5, 1), limits = c(-6, 7))

heatplot(KEGGmink, foldChange=datamink)+
  scale_fill_gradientn(colours = c("#440154FF","#1F968BFF","#FDE725FF"), values=c(0, 0.5, 1), limits = c(-6, 7))

barplot(KEGGmink, showCategory=20)

# GO Site 1

GOmink <- enrichGO(genemink,'org.Hs.eg.db', 
                   ont = "BP", 
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "fdr",
                   qvalueCutoff = 0.05)
head(GOmink)

write.table (GOmink, "GOSigGenes_REFvsSite1.txt", sep="\t")

heatplot(GOmink, foldChange=datamink)+
  scale_fill_gradientn(colours = c("blue","white","red"), values=c(0, 0.5, 1), limits = c(-6, 7))

barplot(GOmink, showCategory=20)

###Site 2

# KEGG Site 2

mink<-read.table("SigGenes_REFvsSite2.txt", header= F)
datamink=mink[,2]
names(datamink)=as.character(mink[,1])
datamink = sort(datamink, decreasing =T)
genemink <-names(datamink)

KEGGmink <- enrichKEGG(gene         = genemink,
                       organism     = 'hsa', 
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "fdr",
                       qvalueCutoff = 0.05)
head(KEGGmink)

write.table (KEGGmink, "KEGGSigGenes_REFvsSite2.txt", sep="\t")

heatplot(KEGGmink, foldChange=datamink)+
  scale_fill_gradientn(colours = c("blue","white","red"), values=c(0, 0.5, 1), limits = c(-6, 7))

heatplot(KEGGmink, foldChange=datamink)+
  scale_fill_gradientn(colours = c("#440154FF","#1F968BFF","#FDE725FF"), values=c(0, 0.5, 1), limits = c(-6, 7))

barplot(KEGGmink, showCategory=20)

# GO Site 2

GOmink <- enrichGO(genemink,'org.Hs.eg.db', 
                   ont = "BP", 
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "fdr",
                   qvalueCutoff = 0.05)
head(GOmink)

write.table (GOmink, "GOSigGenes_REFvsSite2.txt", sep="\t")

heatplot(GOmink, foldChange=datamink)+
  scale_fill_gradientn(colours = c("blue","white","red"), values=c(0, 0.5, 1), limits = c(-6, 7))

barplot(GOmink, showCategory=20)

