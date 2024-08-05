##loop for DE analysis 

##Oliver Brown

# Load packages
library(limma)
library(DESeq2)
library(edgeR)
library(dplyr)
library(combinat)

# Change
setwd("~/Documents/Internship")

# Reads data
cntdata <- read.table("data/all_cnt_0625.txt", 
                      header = TRUE, 
                      na.strings = "NA")

phenodata <- read.table("data/pheno_2_edit.txt", 
                        header = TRUE)

# Filter samples
samplematches <- colnames(cntdata) %in% phenodata$Sample_name
mRNAdata <- cntdata[, samplematches]

# Get gene IDs
gene <- as.data.frame(cntdata$gene)
colnames(gene) <- "gene_id"
mRNAdata <- cbind(gene, mRNAdata)

# Clean gene IDs
mRNAdata$gene_id <- sub("^gene-", "", mRNAdata$gene_id)

# Ensure unique genes in data
mRNAdata <- mRNAdata[!duplicated(mRNAdata$gene_id), ]
row.names(mRNAdata) <- mRNAdata$gene_id

# Remove gene ID column
mRNAdata <- mRNAdata[,-1]

# Factors tissues
phenodata$Tissue <- as.factor(phenodata$Tissue)

# All tissue comparison combinations
tissue_pairs <- combn(levels(phenodata$Tissue), 2, simplify = FALSE)

dds_list <- list()
comparison_results <- list()

for (pair in tissue_pairs) {
  tissue1 <- pair[1]
  tissue2 <- pair[2]
  
  subset_pheno <- phenodata[phenodata$Tissue %in% pair,]
  
  sample_names <- subset_pheno$Sample_name
  sample_names1 <- subset_pheno$Sample_name[subset_pheno$Tissue == tissue1]
  sample_names2 <- subset_pheno$Sample_name[subset_pheno$Tissue == tissue2]
  
  subset_counts <- mRNAdata[, colnames(mRNAdata) %in% sample_names]
  
  # Ensure genes are expressed in both tissue types
  subset_counts <- subset_counts[rowSums(subset_counts[, colnames(subset_counts) %in% sample_names1]) > 0 & 
                                   rowSums(subset_counts[, colnames(subset_counts) %in% sample_names2]) > 0, ]
  
  # Ensure subset_pheno matches subset_counts
  subset_pheno <- subset_pheno[match(colnames(subset_counts), subset_pheno$Sample_name), ]
  
  dds <- DESeqDataSetFromMatrix(countData = subset_counts, 
                                colData = subset_pheno, 
                                design = ~Tissue)
  
  dds <- DESeq(dds, betaPrior = FALSE)
  
  # Adds to list of differential expression data
  dds_list[[paste("Tissue", tissue1, "vs", tissue2, sep = "_")]] <- dds
  
  res <- results(dds, 
                 contrast = c("Tissue", tissue1, tissue2), 
                 pAdjustMethod = "fdr")
  
  # Adds to list of comparison data
  comparison_results[[paste("Tissue", tissue1, "vs", tissue2, sep = "_")]] <- res
}

# Write to CSV files
for (name in names(comparison_results)) {
  write.csv(as.data.frame(comparison_results[[name]]), file = paste0(name, "_results.csv"))
}
