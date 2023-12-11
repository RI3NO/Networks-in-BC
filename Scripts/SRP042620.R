setwd("/Users/matviimykhailichenko/Documents/GitHub/Networks-in-BC/Scripts")

library(recount3)
library(DESeq2)
library(tidyverse)
library(STRINGdb)
library(igraph)
library(biomaRt)
library(org.Hs.eg.db)


# Loading dataset from recount3
human_projects <- available_projects()
proj_info <- subset(
  human_projects,
  project == "SRP042620" & project_type == "data_sources"
)
SRP042620 <- create_rse(proj_info)

SRP042620_expanded_sample_attributes <- expand_sra_attributes(SRP042620)

# Checking if dataset has miRNA, lncRNA and mRNA
rowData(SRP042620_expanded_sample_attributes) %>%
  .$gene_type %>%
  unique()

# Computing read counts from base-pair counts
read_counts <- compute_read_counts(SRP042620_expanded_sample_attributes)

colData(SRP042620_expanded_sample_attributes) %>%
  .$sra_attribute.source_name %>%
  unique()

# Create a function to determine the condition
determine_condition <- function(source_name) {
  if (source_name %in% c("ER+ Breast Cancer Primary Tumor", "Breast Cancer Cell Line", 
                         "Triple Negative Breast Cancer Primary Tumor")) {
    return("cancer")
  } else if (source_name %in% c("Uninvolved Breast Tissue Adjacent to ER+ Primary Tumor", 
                                "Reduction Mammoplasty - No known cancer", 
                                "Uninvolved Breast Tissue Adjacent to TNBC Primary Tumor")) {
    return("normal")
  } else {
    return(NA)  # or some default value
  }
}

# Apply the function to create the new column
colData(SRP042620_expanded_sample_attributes)$condition <- sapply(colData(SRP042620_expanded_sample_attributes)$`sra_attribute.source_name`, determine_condition)

# Making DeSeqDataSet object from matrix
dds <- DESeqDataSetFromMatrix(countData = read_counts,
                              colData = colData(SRP042620_expanded_sample_attributes),
                              design = ~ condition)

# Doing DeSeq function (finding DEGs) on DDS object
deseq_result <- DESeq(subset_dds)


