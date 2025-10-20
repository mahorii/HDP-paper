#R version 4.4.2 (2024-10-31) -- "Pile of Leaves"
#Copyright (C) 2024 The R Foundation for Statistical Computing
#Platform: aarch64-apple-darwin20

#GSVA
#Utilities
library(DESeq2)
library(GSVA)
library(limma)
library(msigdbr)
library(tidyverse)

# Load raw count data (rows: genes, columns: samples)
countData <- read.csv("FeatureCounts_rawcounts_filt.csv", row.names = 1)
metadata <- read.csv("design.csv", row.names = 1)

# Ensure sample order matches in count matrix and metadata
all(colnames(countData) == rownames(metadata))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = countData, colData = metadata, design = ~  batch + GA + sex + BMI_cat + disease)

# Normalize data using variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)
mm <- model.matrix(~disease, colData(vsd))
mat2 <- limma::removeBatchEffect(assay(vsd), batch = interaction(vsd$batch, vsd$sex), design = mm)
assay(vsd) <- mat2
norm_counts <- as.matrix(assay(vsd))

#Convert ensembl gene id to hgnc
library(biomaRt)

# Load Ensembl database
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Convert ENSEMBL to Gene Symbols
gene_conversion <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                         filters = "ensembl_gene_id",
                         values = rownames(norm_counts),
                         mart = mart)

# Remove empty gene symbols
gene_conversion <- gene_conversion[gene_conversion$hgnc_symbol != "", ]

# Map ENSEMBL to Gene Symbols
rownames(norm_counts) <- gene_conversion$hgnc_symbol[match(rownames(norm_counts), gene_conversion$ensembl_gene_id)]

# Remove genes that could not be mapped
norm_counts <- norm_counts[!is.na(rownames(norm_counts)), ]

# Check for duplicate row names
duplicated_rows <- duplicated(rownames(norm_counts))
if (any(duplicated_rows)) {
  print("Duplicate rows found!")
  print(norm_counts[duplicated_rows, ])  # This shows the duplicate rows
  }

# Combine any symbols that are the same by averaging them
norm_counts2 <- as.data.frame(norm_counts) %>%
  rownames_to_column("SYMBOL") %>%
  na.omit() %>%
  distinct() %>%
  group_by(SYMBOL) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  # Averaging duplicate SYMBOLs

# Retain the SYMBOL as rownames and convert to matrix (removing SYMBOL column)
norm_counts2_matrix <- as.matrix(norm_counts2[, -which(names(norm_counts2) == "SYMBOL")])

# Set the rownames of the matrix to the SYMBOL column
rownames(norm_counts2_matrix) <- norm_counts2$SYMBOL

# Load gene sets from MSigDB (Hallmark gene sets as an example)
gene_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Define GSVA parameters (now using the matrix form)
params <- gsvaParam(expr = norm_counts2_matrix, geneSets = gene_sets)

# Perform GSVA
gsva_scores <- gsva(params)

# save the GSVA results
write.csv(gsva_scores, "gsva_scores.csv")

#plot box plot
# Specify the pathway you want to plot
top_pathway <- "HALLMARK_APICAL_JUNCTION"  # Replace with your specific pathway name

# Get the GSVA scores for the specific pathway
top_pathway_scores <- gsva_scores[top_pathway, ]

# Create a data frame for ggplot
plot_data <- data.frame(
  condition = metadata$disease,
  score = top_pathway_scores
)