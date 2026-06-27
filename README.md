# 🐊 RNA-seq analysis workflow for non-lethal scute transcriptomics in *Crocodylus moreletii*

This repository contains the bioinformatic pipeline used to process and analyze RNA-seq data from *Crocodylus moreletii*, integrating transcriptomics tools for biomarker discovery.

## 🔬 Project Overview
The aim of this project is to identify molecular biomarkers associated with environmental stress in wild populations of *C. moreletii*. The analysis includes quality control, mapping, differential expression (DESeq2), functional enrichment (GSEA, GO/KEGG), Weighted Gene Co-expression Network Analysis (WGCNA), single-sample GSEA (ssGSEA) and biomarker selection by machine learning (ML).

### Repository scope
This repository contains the RNA-seq pipeline and scripts used for the transcriptomic analyses of *Crocodylus moreletii* (Publication 1).

*Note: Additional scripts for multiblock integration (mixOmics, coinertia, etc.) will be added in future publications once contaminant datasets become available.*

## 🖥️ System Requirements

- Unix-based system (tested on macOS and Linux)
- R (≥ 4.0.0) with Bioconductor packages
- Python 3.x (for some wrapper scripts)
- Conda (for managing tool environments)
- Software:
  - FastQC 0.12.1
  - Trimmomatic 0.40
  - Bowtie2 2.5.5
  - STAR 2.6.1a_08-27
  - RSeQC 5.0.4
  - Qualimap 2.3
  - HTSeq 2.1.2
  - MultiQC 1.14
  - DESeq2 1.46.0 (R packages)
  - GSVA 2.0.7 (R packages)
  - ClusterProfiler 4.15.0 (R packages)
  - WGCNA 1.73 (R packages)
  - msigdbr 25.1.1 (R packages)
  - pROC 1.19.0.1 (R packages)
  - glmnet 4.1.10 (R packages)
  - org.Hs.eg.db 3.20.0 (R packages)
  - AnnotationDbi 1.68.0 (R packages)


## 📁 Folder Structure

Each folder corresponds to a pipeline step:

- `01_QC_Raw/`: Raw data quality check
- `02_Trimming/`: Adapter and quality trimming
- `03_rRNA_Removal/`: rRNA filtering via Bowtie2
- `04_Mapping_STAR/`: STAR genome alignment
- `05_Mapping_QC/`: Quality control of alignments
- `06_Counts_HTSeq/`: Gene-level quantification
- `07_MultiQC_Summary_QC/`: Summarize analysis results for multiple tools and samples in a single report
- `09_DESeq2_Analysis/`: Differential expression and PCA
- `10_Gene_Annotation/`: Gene annotation including gene names, entrez id, ensembl id, go id and kegg id
- `11_PC1_loading-based_enrich_analysis/`: Visualize the overall transcriptional structure and interpret the biological meaning of the first principal component (PC1)
- `12_GSEA_KEGG-Enrichment_Supervised/`: Confirm the pathways of PC1 by doing supervised analysis between six sites (wild) and captive group
- `13_WGCNA/`: Identify co-expressed gene modules using top variable genes
- `14_ssGSEA/`: Assess individual-level activation of specific stress pathways and explore site-specific stress signatures
- `15_Biomarker_Selection/`: Identify potential biomarker genes indicative of environmental stress by machine learning.


## ⚙️ Reproducibility

All scripts are annotated and modular. You can reproduce each step independently. System-specific scripts (Mac/Linux) are included where applicable.

## 📌 Notes

- **Data**: FASTQ, count tables, and DE results are not included in this repository. They may be available upon request after publication.
- **License**: MIT License – you may reuse and adapt the code with attribution.
- **Status**: This repository is under development and will be made public upon manuscript acceptance.

### 🧪 Additional Modules

#### Midori Tool – Transcript Coverage QC

This optional module performs transcript coverage quality checks using mapped BAM files and GTF annotations. It uses `bedtools`, `awk`, and `R` for normalization and visualization.

📁 Directory: `08_Midori_QC/`  
📄 Documentation: [`README_Midori_Tool_transcript_coverage.md`](./08_Midori_QC/README_Midori_Tool_transcript_coverage.md)

**To run:**

```bash```

```Rscript coverage_analysis.R -g input.gtf -b input.bam -o output_plot.png```

## 👩‍🔬 Authors

**Asela Marisol Buenfil-Rojas**  
Center for Marine Environmental Studies (CMES) – Ehime University  
ORCID: [0000-0003-0999-7428](https://orcid.org/0000-0003-0999-7428)  
Contact: buenfil-rojas.am@outlook.com

**Midori Iida**  
Kyushu Institute of Technology  
ORCID: [0000-0002-5832-6148](https://orcid.org/0000-0002-5832-6148)

**Mirella Kanerva**  
Center for Marine Environmental Studies (CMES) - Ehime University  
ORCID: [0000-0002-2463-2497](https://orcid.org/0000-0002-2463-2497)


## 📖 Citation

Once published, please cite the associated article and this repository DOI: 10.5281/zenodo.18724008 (via Zenodo).
