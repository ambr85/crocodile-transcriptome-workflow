###########################
#  Transcripts coverage   #
###########################

#Usage
#Rscript coverage_analysis.R -g input.gtf -b input.bam -o output_plot.png

# 必要なパッケージをロード
library(GenomicAlignments)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(BiocParallel)
library(argparse)
library(purrr)

# Process command line arguments
parser <- ArgumentParser(description = "Process and visualize transcript coverage")
parser$add_argument("-g", "--gtf", required = TRUE, help = "Input GTF file")
parser$add_argument("-b", "--bam", required = TRUE, help = "Input BAM file")
parser$add_argument("-o", "--output", required = TRUE, help = "Output plot file")
args <- parser$parse_args()

# Get input and output file names
gtf_file <- args$gtf
bam_file <- args$bam
output_file <- args$output

# Set temporary file name
coverage_bed_file <- "coverage.bed"
input_bed_file    <- "input_bed_file.bed"

### 1. rrun bedtools & awk ##############################################
system(paste("bedtools coverage -a", gtf_file, "-b", bam_file, ">", coverage_bed_file))
system(paste(
  "less", coverage_bed_file, "|",
  "awk -F'[\\t]' '$3 == \"exon\" {split($9, arr, \";\"); $9 = arr[2]; gsub(/transcript_id \"/, \"\", $9); gsub(/^ /, \"\", $9); gsub(/\"/, \"\", $9); print $1, $4, $5, $9, $7, $11}' >",
  input_bed_file
))

### 2. load data ##########################################################
# read bed
bed_file <- input_bed_file
bed      <- read.delim(bed_file, header = FALSE, comment.char = "#", sep = " ")
colnames(bed) <- c("seqname", "start", "end", "id", "strand", "coverage")

### 2. filter transcripts ################################################
# Calculate total coverage for each id and exclude those with a total of 0
filtered_bed <- bed %>%
  group_by(id) %>%
  filter(sum(coverage) > 0) %>%
  ungroup()

# Exclude total lengths less than 100 bp
filtered_len_bed <- filtered_bed %>% 
  group_by(id) %>% mutate(range=abs(start-end)+1, 
                          total_length = sum(range)) %>%
  filter(total_length > 100) 

# grouping the transripts by gene
grouped_filtered_len_bed <- filtered_len_bed %>%
  mutate(
    transcript_length_group = factor(
      dplyr::case_when(
        total_length <= 500 ~ '0-500',
        total_length <= 1000 ~ '501-1000',
        total_length <= 2000 ~ '1001-2000',
        total_length <= 5000 ~ '2001-5000',
        total_length <= 10000 ~ '5001-10000',
        TRUE ~ '10001+'
      ),
      levels = c('0-500', '501-1000', '1001-2000', '2001-5000', '5001-10000', '10001+')
    )
  )

### 3. get coverage in transcripts #######################################
result <- grouped_filtered_len_bed %>%
  group_by(id) %>%
  summarise(
    quantile_coverage = list({
      # Processed within each group
      x <- seq(min(start), max(end))
      v <- start
      if (first(strand) == "+") {
        quantile(coverage[findInterval(x, sort(v))], probs = seq(0, 1, length.out = 101))
      } else {
        quantile(rev(coverage)[findInterval(x, sort(v))], probs = seq(0, 1, length.out = 101))
      }
    })
  )

# Expand results
result_expanded <- result %>%
  mutate(quantile_coverage = map(quantile_coverage, as.numeric)) %>%
  unnest_wider(quantile_coverage, names_sep = "_")

### 5. normalize coverage in transcripts #######################################
# combined results by id
combined_result <- left_join(result_expanded, grouped_filtered_len_bed %>%
                               select(id, transcript_length_group) %>%
                               distinct(), by = "id")

gropuedsum_combined_result <-
combined_result %>% 
  group_by(transcript_length_group)  %>%
  summarise(
    across(starts_with("quantile_"), sum, .names = "sum_{col}")
  )

# Normalize quantile_ columns in each row (so that the maximum value is 100)
normalized_result <-apply(gropuedsum_combined_result[ 2:101], 1, function(x){x/max(x)*100})
colnames(normalized_result) <- gropuedsum_combined_result[,1]$transcript_length_group
rownames(normalized_result) <- gsub("sum_quantile_coverage_", "", rownames(normalized_result))

### 6. Average coverage per group is illustrated ##################################
# Convert data to long format (tidy format)
plot_data <- melt(normalized_result)

# Set Var2 as factor and make the order explicit
plot_data$Var2 <- factor(
  plot_data$Var2,
  levels = c('0-500', '501-1000', '1001-2000', '2001-5000', '5001-10000', '10001+')
)

#  Average coverage per group is illustrated
plot <- 
ggplot(plot_data, aes(x = Var1, y = value, color = Var2, group = Var2)) +
  geom_line() +
  labs(title = "Average Coverage per Segment by Gene Length Group",
       x = "Segment (1-100)",
       y = "Mean Coverage (%)",
       color = "Gene Length Group") +
  theme_minimal()

# save plot
ggsave(output_file, plot)

##Author: Iida Midori
