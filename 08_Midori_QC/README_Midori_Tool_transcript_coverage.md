# Midori Tool - Transcript Coverage Analysis

This script processes BAM and GTF files to analyze transcript coverage across gene length groups. It uses bedtools and awk for intermediate data processing, and R packages for coverage calculation, filtering, and visualization.

## Pre-requisites

Ensure you have the following installed:

- Software:

  - bedtools
  - awk
  - R (with the necessary libraries)

-	R Libraries:

```install.packages(c("dplyr", "tidyr", "ggplot2", "reshape2", "argparse"))```

```BiocManager::install(c("GenomicAlignments", "GenomicRanges", "BiocParallel"))```

### Input Files

1.	GTF File: A file containing gene annotation data.
2.	BAM File: A file containing mapped read data.

### Output
The script generates a plot showing the average coverage per segment across different gene length groups.

### Usage
Run the script as follows:

```bash```

```Rscript coverage_analysis.R -g input.gtf -b input.bam -o output_plot.png```

Arguments

- -g or --gtf: Path to the input GTF file.
- -b or --bam: Path to the input BAM file.
- -o or --output: Path to save the output plot (PNG format).

### Workflow

1.	Data Processing:
    - The script uses bedtools and awk to process GTF and BAM files into a filtered BED file.
    - Only exonic regions are extracted, and transcript IDs are formatted.

2.	Filtering:
    - Transcripts with zero coverage are removed.
    - Transcripts with a total length ≤ 100bp are excluded.

3.	Gene Length Grouping:
    - Transcripts are grouped into length categories:
      - 	0-500
      - 	501-1000
      - 	1001-2000
      - 	2001-5000
      - 	5001-10000
      - 	10001+

4.	Coverage Quantile Calculation:
    - Coverage values are summarized into 100 quantiles for each transcript.

5.	Normalization:
    - Coverage values are normalized such that the maximum value per row is scaled to 100.

6.	Visualization:
    - The plot shows the average coverage per segment for each gene length group.

### Example Workflow

- Example input files<br>
  - bash<br> input.gtf = "example.gtf"<br> input.bam = "example.bam"<br> output_plot = "coverage_plot.png"<br>

- Run the script<br>
```Rscript coverage_analysis.R -g $input.gtf -b $input.bam -o $output_plot```<br>
```# The resulting plot will be saved as coverage_plot.png.```

- Example output files
  - The plot shows:
   - X-axis: Segment (1-100).
   - Y-axis: Normalized coverage (%) across all transcripts in each group.
   - Color: Gene length groups.
 
#### Notes
 - Ensure bedtools and awk are properly installed and available in your PATH.
 - Large BAM files may take significant time to process. Consider optimizing resources or working with subsets.
