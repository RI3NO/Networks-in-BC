# The dataset has lncRNAs, miRNAs and mRNAs

setwd("/Users/matviimykhailichenko/Library/CloudStorage/OneDrive-UniwersytetWrocławski/Carrier/Articles/Networks_of_Biomarkers_in_BC/Code/Networks-in-BC")

library(recount3)
library(DESeq2)
library(tidyverse)
library(STRINGdb)
library(igraph)
library(biomaRt)
library(org.Hs.eg.db)
library(DEGreport)


# Loading dataset from recount3
human_projects <- available_projects()
proj_info <- subset(
  human_projects,
  project == "SRP042620" & project_type == "data_sources"
)
SRP042620 <- create_rse(proj_info)
rownames(SRP042620) <- rowData(SRP042620)$gene_name
SRP042620_expanded_sample_attributes <- expand_sra_attributes(SRP042620)

# Computing read counts from base-pair counts
assay(SRP042620_expanded_sample_attributes) <- compute_read_counts(SRP042620_expanded_sample_attributes)

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




