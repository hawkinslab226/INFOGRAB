
# Load required libraries
library(dplyr)
library(rtracklayer)
library(VariantAnnotation)
library(knitr)
library(kableExtra)
library(tispec)

# Set working directory
setwd("~/Documents/Internship/InfoGrab")

# Load data
comparison_results <- readRDS("data/comparison_results.rds")
RPKM_data <- read.delim("data/logRPKM0625_filter.txt")
phenodata <- read.csv("data/Pheno_data071523.csv")
cnt_data <- read.delim("data/all_cnt_0625.txt")

# Preliminary Data
num_genes <- nrow(RPKM_data)
num_samples <- nrow(phenodata)
num_tissues <- length(unique(phenodata$Tissue)) # Assuming 'tissue' column exists in phenodata
num_comparisons <- length(comparison_results)

# Differential Expression (DE) Analysis
significance_threshold <- 0.05
unique_de_genes <- unique(unlist(lapply(comparison_results, function(comparison_data) {
  rownames(comparison_data[!is.na(comparison_data$padj) & comparison_data$padj < significance_threshold, ])
})))
total_unique_de_genes <- length(unique_de_genes)
average_de_genes <- mean(sapply(comparison_results, function(comparison_data) {
  sum(!is.na(comparison_data$padj) & comparison_data$padj < significance_threshold)
}))
median_de_genes <- median(sapply(comparison_results, function(comparison_data) {
  sum(!is.na(comparison_data$padj) & comparison_data$padj < significance_threshold)
}))
sd_de_genes <- sd(sapply(comparison_results, function(comparison_data) {
  sum(!is.na(comparison_data$padj) & comparison_data$padj < significance_threshold)
}))
min_de_genes <- min(sapply(comparison_results, function(comparison_data) {
  sum(!is.na(comparison_data$padj) & comparison_data$padj < significance_threshold)
}))
max_de_genes <- max(sapply(comparison_results, function(comparison_data) {
  sum(!is.na(comparison_data$padj) & comparison_data$padj < significance_threshold)
}))
total_genes_analyzed <- length(unique(unlist(lapply(comparison_results, rownames))))
pairs_with_significant_de_genes <- sum(sapply(comparison_results, function(comparison_data) {
  any(!is.na(comparison_data$padj) & comparison_data$padj < significance_threshold)
}))

# Tissue-Specific (TS) Analysis
tau_scores <- calcTau(RPKM_data)
n_tau_greater_zero <- length(tau_scores$tau[tau_scores$tau > 0])
significant_tau_count <- length(tau_scores$tau[tau_scores$tau > 0.85])
proportion_significant_tau <- significant_tau_count / length(tau_scores$tau)
mean_tau_score <- mean(tau_scores$tau, na.rm = TRUE)
median_tau_score <- median(tau_scores$tau, na.rm = TRUE)
sd_tau_score <- sd(tau_scores$tau, na.rm = TRUE)

# Allele-Specific Expression (ASE) Analysis
gene_ase_path <- "browser_data_for_app/gene_ase_converted_cleaned_with_strand.bed"
gene_ase_data <- rtracklayer::import(gene_ase_path, format = "bed")
gene_ase_df <- as.data.frame(gene_ase_data)
n_gene_ase <- length(unique(gene_ase_df$name))
proportion_ase_gene <- n_gene_ase / length(cnt_data$gene)
vcf_file <- "browser_data_for_app/ASE.SNP.gallus.updated02.vcf"
vcf_data <- readVcf(vcf_file, genome = "Galgal6")
snp_data <- rowRanges(vcf_data)
fixed_fields_s4 <- VariantAnnotation::fixed(vcf_data)
fixed_fields <- as.data.frame(fixed_fields_s4)
snp_df <- data.frame(
  chr = as.character(seqnames(snp_data)),
  start = start(snp_data),
  end = end(snp_data),
  REF = fixed_fields$REF,
  ALT = sapply(fixed_fields$ALT, function(x) paste(x, collapse = ","))
)
num_ase_snps <- length(snp_df$REF)
avg_ase_snps_per_gene <- num_ase_snps / n_gene_ase
median_ase_snps_per_gene <- median(table(gene_ase_df$name))

# Create the main summary table
# Create the main summary table without the initial NA row
# Create the main summary table without headers
# Create the main summary table without headers
table_of_values <- data.frame(
  Measure = c(
    # Preliminary Data
    "Number of Genes",
    "Number of Samples",
    "Number of Tissues",
    "Number of Comparisons",
    
    # DE Analysis
    "Total Unique DE Genes",
    "Average DE Genes per Pair",
    "Median DE Genes per Pair",
    "Standard Deviation of DE Genes per Pair",
    "Minimum DE Genes in a Pair",
    "Maximum DE Genes in a Pair",
    "Total Number of Genes Analyzed",
    "Number of Pairs with Significant DE Genes",
    
    # TS Analysis
    "Number of Genes with Tau Score Greater than 0",
    "Number of Significant TS Genes (Tau > 0.85)",
    "Proportion of Significant TS Genes (>0.85)",
    "Mean Tau Score",
    "Median Tau Score",
    "Standard Deviation of Tau Scores",
    
    # ASE Analysis
    "Number of ASE Genes",
    "Proportion of ASE Genes",
    "Number of ASE SNPs",
    "Average Number of ASE SNPs per Gene",
    "Median Number of ASE SNPs per Gene"
  ),
  Value = c(
    # Preliminary Data values
    format(num_genes, big.mark = ",", scientific = FALSE),
    format(num_samples, big.mark = ",", scientific = FALSE),
    format(num_tissues, big.mark = ",", scientific = FALSE),
    format(num_comparisons, big.mark = ",", scientific = FALSE),
    
    # DE Analysis values
    format(total_unique_de_genes, big.mark = ",", scientific = FALSE),
    format(round(average_de_genes, 4), nsmall = 4), # Average values with 4 decimal places
    format(median_de_genes, big.mark = ",", scientific = FALSE),
    format(round(sd_de_genes, 4), nsmall = 4),      # Standard deviation with 4 decimal places
    format(min_de_genes, big.mark = ",", scientific = FALSE),
    format(max_de_genes, big.mark = ",", scientific = FALSE),
    format(total_genes_analyzed, big.mark = ",", scientific = FALSE),
    format(pairs_with_significant_de_genes, big.mark = ",", scientific = FALSE),
    
    # TS Analysis values
    format(n_tau_greater_zero, big.mark = ",", scientific = FALSE),
    format(significant_tau_count, big.mark = ",", scientific = FALSE),
    format(round(proportion_significant_tau, 4), nsmall = 4),
    format(round(mean_tau_score, 4), nsmall = 4),
    format(round(median_tau_score, 4), nsmall = 4),
    format(round(sd_tau_score, 4), nsmall = 4),
    
    # ASE Analysis values
    format(n_gene_ase, big.mark = ",", scientific = FALSE),
    format(round(proportion_ase_gene, 4), nsmall = 4),
    format(num_ase_snps, big.mark = ",", scientific = FALSE),
    format(round(avg_ase_snps_per_gene, 4), nsmall = 4),
    format(median_ase_snps_per_gene, big.mark = ",", scientific = FALSE)
  )
)

# Display table without headers and with integer formatting
kable(table_of_values, col.names = NULL, align = "l", format = "html") %>%
  kable_styling(full_width = F, bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  pack_rows("Datasets Overview", 1, 4) %>%
  pack_rows("Differential Expression (DE) Analysis", 5, 12) %>%
  pack_rows("Tissue-Specific (TS) Analysis", 13, 18) %>%
  pack_rows("Allele-Specific Expression (ASE) Analysis", 19, 23) %>%
  column_spec(1, bold = TRUE) %>%
  column_spec(2, width = "6em")


