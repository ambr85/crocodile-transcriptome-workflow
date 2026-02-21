# 📚 Software References

This repository uses several bioinformatic tools and R packages. Below are the recommended references for each.

## 🔬 RNA-seq Preprocessing and Alignment

- **FastQC**  
  Andrews, S. (2010). *FastQC: a quality control tool for high throughput sequence data*.
  [https://www.bioinformatics.babraham.ac.uk/projects/fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc)

- **Trimmomatic**  
  Bolger, A. M., Lohse, M., & Usadel, B. (2014). *Trimmomatic: a flexible trimmer for Illumina sequence data*. Bioinformatics, 30(15), 2114–2120. doi:10.1093/bioinformatics/btu170

- **SortMeRNA**  
  Kopylova, E., Noé, L., & Touzet, H. (2012). *SortMeRNA: fast and accurate filtering of ribosomal RNAs in metatranscriptomic data*. Bioinformatics, 28(24), 3211–3217. doi:10.1093/bioinformatics/bts611

- **Bowtie2**  
  Langmead, B., & Salzberg, S. L. (2012). *Fast gapped-read alignment with Bowtie 2*. Nature Methods, 9, 357–359. doi:10.1038/nmeth.1923

- **STAR**  
  Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., & Gingeras, T. R. (2013). *STAR: ultrafast universal RNA-seq aligner*. Bioinformatics, 29(1), 15–21. doi:10.1093/bioinformatics/bts635

- **HTSeq**  
  Anders, S., Pyl, P. T., & Huber, W. (2015). *HTSeq—a Python framework to work with high-throughput sequencing data*. Bioinformatics, 31(2), 166–169. doi:10.1093/bioinformatics/btu638

- **MultiQC**  
  Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). *MultiQC: Summarize analysis results for multiple tools and samples in a single report*. Bioinformatics, 32(19), 3047–3048. doi:10.1093/bioinformatics/btw354
  
- **Qualimap**  
  Okonechnikov, K., Conesa, A., & García-Alcalde, F. (2016). *Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data*. Bioinformatics, 32(2), 292–294. doi:10.1093/bioinformatics/btv566
 
- **RSeQC**  
  Wang, L., Wang, S., & Li, W. (2012). *RSeQC: quality control of RNA-seq experiments*. Bioinformatics, 28(16), 2184–2185. doi:10.1093/bioinformatics/bts356
 
## 📈 Differential Expression, Network Analysis and Biomarker Selection

- **DESeq2**  
  Love, M. I., Huber, W., & Anders, S. (2014). *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2*. Genome Biology, 15, 550. doi:10.1186/s13059-014-0550-8

- **WGCNA**  
  Langfelder, P., Horvath, S. (2008). *WGCNA: an R package for weighted correlation network analysis*. BMC Bioinformatics, 559.
doi:10.1186/1471-2105-9-559

  Langfelder, P., Horvath, S. (2012). *Fast R Functions for Robust Correlations and Hierarchical Clustering*. Journal of Statistical Software, 46(11), 1-17.
  [https://www.jstatsoft.org/v46/i11](https://www.jstatsoft.org/v46/i11)
  
- **pROC**  
  Robin, X., Turck, N., Hainard, A., Tiberti, N., Lisacek, F., Sanchez J.-C. Müller, M. (2011). *pROC: an open-source package for R and S+ to analyze and compare ROC curves*. BMC Bioinformatics, 12, 77.
doi:10.1186/1471-2105-12-77
 
- **caret**  
  Kuhn, M. (2008). *Building Predictive Models in R Using the caret Package*. Journal of Statistical Software, 28(5), 1–26.
doi:10.18637/jss.v028.i05
  
- **glmnet**  
  Friedman, J., Hastie, T., Tibshirani, R. (2010). *Regularization Paths for Generalized Linear Models via Coordinate Descent*. Journal of Statistical Software, 33(1), 1-22.
doi:10.18637/jss.v033.i01

  Tay, J.K., Narasimhan, B., Hastie, T. (2023). *Elastic Net Regularization Paths for All Generalized Linear Models*. Journal of Statistical Software, 106(1), 1-31. doi:10.18637/jss.v106.i01

## 🧬 Functional Enrichment

- **clusterProfiler**  
  Yu, G., Wang, L. G., Han, Y., & He, Q. Y. (2012). *clusterProfiler: an R package for comparing biological themes among gene clusters*. OMICS, 16(5), 284–287. doi:10.1089/omi.2011.0118

- **AnnotationDbi**  
  Pagès, H., Carlson, M., Falcon, S., Li, N. (2024). *AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor*.
doi:10.18129/B9.bioc.AnnotationDbi
  
- **GSVA**  
  Hänzelmann, S., Castelo, R. and Guinney, A. (2013). *GSVA: gene set variation analysis for microarray and RNA-seq data*. BMC Bioinformatics, 14:7.
doi:10.1186/1471-2105-14-7

- **msigdbr**  
  Dolgalev, I. (2025). *msigdbr: MSigDB Gene Sets for Multiple Organisms in a Tidy Data Format*. R package version 25.1.1
  [http://biorxiv.org/content/early/2016/06/20/060012](https://CRAN.R-project.org/package=msigdbr)
  
- **GSEA**  
  Korotkevich, G., Sukhov, V., & Sergushichev, A. (2019). *Fast gene set enrichment analysis*. bioRxiv. doi:10.1101/060012  
  [http://biorxiv.org/content/early/2016/06/20/060012](http://biorxiv.org/content/early/2016/06/20/060012)

## 🧪 Additional Tools

- **Midori Tool (Transcript Coverage QC)**  
  Iida, M. (2024). *Midori Tool: Coverage-based transcript QC*. [Documentation included in repo]
