# 🐊 Crocodile Transcriptome Workflow
Transcriptomic analysis of *Crocodylus moreletii*. RNAseq, DESeq2 and MixOmics Workflow

This repository contains the bioinformatic pipeline used to process and analyze RNA-seq data from *Crocodylus moreletii*, integrating transcriptomics with contaminant exposure data for biomarker discovery.

## 🔬 Project Overview

The aim of this project is to identify molecular biomarkers associated with environmental contaminants in wild populations of *C. moreletii*. The analysis includes quality control, mapping, differential expression (DESeq2), multiblock integration (mixOmics DIABLO & CIA), and functional enrichment (GSEA, GO/KEGG).

## 🖥️ System Requirements

- Unix-based system (tested on macOS and Linux)
- R (≥ 4.0.0) with Bioconductor packages
- Python 3.x (for some wrapper scripts)
- Conda (for managing tool environments)
- Software:
  - FastQC
  - Trimmomatic
  - Bowtie2
  - STAR
  - RSeQC
  - Qualimap
  - HTSeq
  - MultiQC
  - DESeq2 (R packages)
  - GSEA (R packages)
  - ClusterProfiler (R packages)
  - CIA (Co-inertia Analysis) (R packages)
  - MixOmics (R packages)


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
- `10_GSEA_Enrichment/`: Pathway and KEGG enrichment analysis
- `11_mixOmics_CIA/`: Integration analyses (DIABLO, CIA)


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

**Mirella Kanerva**  
Affiliation: [por confirmar]  
ORCID: [0000-0002-2463-2497](https://orcid.org/0000-0002-2463-2497)

**Midori Iida**  
Kyushu Institute of Technology  
ORCID: [0000-0002-5832-6148](https://orcid.org/0000-0002-5832-6148)


## 📖 Citation

Once published, please cite the associated article and this repository DOI (via Zenodo).
