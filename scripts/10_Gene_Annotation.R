#########################################################
#Program: Gene Annotacion after DESeq2
#Author: Marisol Buenfil Rojas
#Date: 2025
#Notes: Example Site 1
#########################################################

## Install libraries
install.packages("BiocManager")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")

library("AnnotationDbi")
library("org.Hs.eg.db")
library("AnnotationHub")

# Load data
directory <- "YourPath/Gene_Annotation"
setwd(directory)
getwd()

######################################################################

deseq_Results <- read.csv("Significant_Genes_REFvsSite1.csv", header = TRUE,
                          as.is=TRUE,row.names=1)
head(deseq_Results)

# Names of the annotation
columns(org.Hs.eg.db) #human annotation
keytypes(org.Hs.eg.db)

# Annotation in different formats: Gene name, ENSEMBL, Entrez ID, Gene Onthology and KEGG pathways
deseq_Results$genename <- mapIds(org.Hs.eg.db, keys=row.names(deseq_Results),column="GENENAME", keytype="SYMBOL", multiVals="first") #Name of the gene column
deseq_Results$entrezid <- mapIds(org.Hs.eg.db, keys=row.names(deseq_Results), column="ENTREZID",keytype="SYMBOL",multiVals="first")#Inserts Entrez ID column
deseq_Results$ensembl <- mapIds(org.Hs.eg.db, keys=row.names(deseq_Results), column="ENSEMBL",keytype="SYMBOL",multiVals="first")#Ensembl column
deseq_Results$go <- mapIds(org.Hs.eg.db, keys=row.names(deseq_Results), column="GO",keytype="SYMBOL",multiVals="first")#Gene onthology column
deseq_Results$path <- mapIds(org.Hs.eg.db, keys=row.names(deseq_Results), column="PATH",keytype="SYMBOL",multiVals="first")#KEGG pathway column

# Check the new annotation columns added
head(deseq_Results)

# Save the new table with the annotations
write.csv(deseq_Results, file="FinalSigGenes_REFvsSite1.csv")

