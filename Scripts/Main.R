# Change to directory where you cloned the repository
setwd("/Users/matviimykhailichenko/Library/CloudStorage/OneDrive-UniwersytetWrocławski/Carrier/Articles/Networks_of_Biomarkers_in_BC/Code/Networks-in-BC")

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(enrichR)
library(recount3)
library(DESeq2)
library(tidyverse)
library(STRINGdb)
library(igraph)
library(biomaRt)
library(org.Hs.eg.db)
library(DEGreport)
library(enrichplot)
library(clusterProfiler)
library(magrittr)
source("Scripts/Functions.R")

dds_SRP042620 <- retrieve_dds("SRP042620")
rowData(dds_SRP042620)$gene_type %>%
  unique()
# Create a vector of tissue types you want to keep
# tissue_types_to_keep <- c("ER+ Breast Cancer Primary Tumor", "Uninvolved Breast Tissue Adjacent to ER+ Primary Tumor")

# Subset the cells based on the tissue types
# rse_ER_pos <- dataset[, colData(dataset)$sra_attribute.tissue %in% tissue_types_to_keep]

# rse_ER_pos %>%
#  colData() %>%
#  .$sra_attribute.tissue

# colData_ER_pos <- colData(rse_ER_pos) %>% 
#  as.data.frame() %>%
#  .$condition %>%
#  data.frame(condition=.)

# dds_ER_pos <- DESeqDataSetFromMatrix(countData = assay(rse_ER_pos),
#                                      colData = colData_ER_pos,
#                                     design = ~ condition)

#Performing variance stabilizing transformation (vst)
# vsd_ER_pos <- vst(dds_ER_pos)

# Perform PCA
# PCA_plot_ER_pos <- plotPCA(vsd_ER_pos, intgroup = "condition")
# PCA_plot_ER_pos

# dataset %>%
#   colData() %>%
#   .$sra_attribute.tissue

# Filtering TNBC cells
tissue_types_to_keep <- c("Uninvolved Breast Tissue Adjacent to TNBC Primary Tumor", "Triple Negative Breast Cancer Primary Tumor")

dds_TNBC <- dds_SRP042620[, colData(dds_SRP042620)$sra_attribute.tissue %in% tissue_types_to_keep]

"""

# Creating a plot to ranking genes to make a treshold  

# Calculate the average expression of genes across samples
avg_expression <- rowMeans(counts(dds_TNBC))



# Create a data frame with gene IDs and average expression values
avg_expression_df <- data.frame(Gene_ID = rownames(counts(dds_TNBC)), Avg_Expression = avg_expression)

# Create a scatter plot
ggplot(avg_expression_df, aes(x = Gene_ID, y = Avg_Expression)) +
  geom_point() +
  labs(x = "Gene ID", y = "Average Expression")

"""



#Performing variance stabilizing transformation (vst)
vsd_TNBC <- vst(dds_TNBC)

# Perform PCA
PCA_plot_TNBC <- plotPCA(vsd_TNBC, intgroup = "condition")
PCA_plot_TNBC

# dataset %>%
#  colData() %>%
#  .$sra_attribute.tissue

# Create a vector of tissue types you want to keep
# tissue_types_to_keep <- c("Breast Cancer Cell Line", "Reduction Mammoplasty - No known cancer")

# Subset the cells based on the tissue types
# rse_BC_cell_line <- dataset[, colData(dataset)$sra_attribute.tissue %in% tissue_types_to_keep]

# rse_BC_cell_line %>%
#  colData() %>%
#  .$sra_attribute.tissue

#colData_BC_cell_line <- colData(rse_BC_cell_line) %>% 
#  as.data.frame() %>%
#  .$condition %>%
#  data.frame(condition=.)

#dds_BC_cell_line <- DESeqDataSetFromMatrix(countData = assay(rse_BC_cell_line),
#                                   colData = colData_BC_cell_line,
#                                   design = ~ condition)

#Performing variance stabilizing transformation (vst)
#vsd_BC_cell_line <- vst(dds_BC_cell_line)

# Perform PCA
#PCA_plot_BC_cell_line <- plotPCA(vsd_BC_cell_line, intgroup = "condition")
#PCA_plot_BC_cell_line

# ER_pos DEGs Analysis

#dds_ER_pos <- DESeq(dds_ER_pos)

#results_ER_pos <- results(dds_ER_pos, contrast = c('condition','cancer','normal'))

#DEGs_ER_pos <- subset(results_ER_pos, padj < 0.05 & abs(log2FoldChange) > 2 & baseMean > 1000 )

#DEGs_ER_pos <- DEGs_ER_pos[order(DEGs_ER_pos$log2FoldChange, decreasing = TRUE),]

#DEGs_dds_ER_pos <- dds_ER_pos[rownames(dds_ER_pos) %in% rownames(DEGs_ER_pos),]

#DEGs_counts_ER_pos<-assay(DEGs_dds_ER_pos)[rownames(DEGs_dds_ER_pos), rownames(colData(DEGs_dds_ER_pos))] #sig genes x samples


# Function to generate heatmap, taking a DEGs_counts matrix and a dds object
#Heatmap_ER_pos <- generate_heatmap(DEGs_counts_ER_pos, dds_ER_pos)
#Heatmap_ER_pos

# TNBC DEGs Analysis

dds_TNBC <- DESeq(dds_TNBC)
results_TNBC <- results(dds_TNBC, contrast = c('condition','cancer','normal'))
DEGs_TNBC <- subset(results_TNBC, padj < 0.05)


# Print the row containing "MALAT1"
# print(result_row)
# DEGs_TNBC <- DEGs_TNBC[order(DEGs_TNBC$log2FoldChange, decreasing = TRUE),]
#DEGs_dds_TNBC <- dds_TNBC[rownames(dds_TNBC) %in% rownames(DEGs_TNBC),]
#DEGs_counts_TNBC <-assay(DEGs_dds_TNBC)[rownames(DEGs_dds_TNBC), rownames(colData(DEGs_dds_TNBC))] #sig genes x samples


# Using tresholds of padj < 0.05 & abs(log2FoldChange) > 3 & baseMean > 2000 to generate a Heatmap
Heatmap_DEGs_TNBC <- subset(results_TNBC, padj < 0.05 & abs(log2FoldChange) > 3 & baseMean > 2000) %>%
  .[order(Heatmap_DEGs_TNBC$log2FoldChange, decreasing = TRUE),]
Heatmap_DEGs_dds_TNBC <- dds_TNBC[rownames(dds_TNBC) %in% rownames(Heatmap_DEGs_TNBC),]
Heatmap_TNBC <- generate_heatmap(assay(Heatmap_DEGs_dds_TNBC), Heatmap_DEGs_dds_TNBC)
Heatmap_TNBC



# BC_cell_line DEGs Analysis

#dds_BC_cell_line<- DESeq(dds_BC_cell_line)
#results_BC_cell_line <- results(dds_BC_cell_line, contrast = c('condition','cancer','normal'))
#DEGs_BC_cell_line <- subset(results_BC_cell_line, padj < 0.05 & abs(log2FoldChange) > 2 & baseMean > 1000 )
#DEGs_BC_cell_line <- DEGs_BC_cell_line[order(DEGs_BC_cell_line$log2FoldChange, decreasing = TRUE),]
#DEGs_dds_BC_cell_line <- dds_BC_cell_line[rownames(dds_BC_cell_line) %in% rownames(DEGs_BC_cell_line),]
#DEGs_counts_BC_cell_line <-assay(DEGs_dds_BC_cell_line)[rownames(DEGs_dds_BC_cell_line), rownames(colData(DEGs_dds_BC_cell_line))] #sig genes x samples

#Heatmap_BC_cell_line  <- generate_heatmap(DEGs_counts_BC_cell_line , dds_BC_cell_line )
#Heatmap_BC_cell_line 


# Enrichment analysis
# ER+
#enrichment_ER_pos <- enrichGO(gene = rownames(DEGs_ER_pos), 
#                OrgDb = org.Hs.eg.db, # replace with the appropriate organism database
#                keyType = "SYMBOL", # change based on your gene ID type
#                ont = "BP", # Biological Process (BP), Molecular Function (MF), or Cellular Component (CC)
#                pAdjustMethod = "BH", # method for adjusting p values
#                qvalueCutoff = 0.05, # cutoff for q value
#                readable = TRUE) # convert entrez IDs to gene symbols

#enrichment_pairwise_ER_pos <- pairwise_termsim(enrichment_ER_pos)

# Barplot
#barplot_ER_pos <- barplot(enrichment_ER_pos, showCategory=20)
#barplot_ER_pos

# Dotplot
#dotplot_ER_pos <- dotplot(enrichment_ER_pos, showCategory = 15) 
#dotplot_ER_pos

# Treeplot
#treeplot_ER_pos <- treeplot(enrichment_pairwise_ER_pos)
#treeplot_ER_pos

# Enrichment map
#emapplot_ER_pos <- emapplot(enrichment_pairwise_ER_pos)
#emapplot_ER_pos

# Cnetplot
#cnetplot_ER_pos <- cnetplot(enrichment_ER_pos, cex_label_gene=0.5)
#cnetplot_ER_pos

# Upsetplot
#upsetplot_ER_pos <- upsetplot(enrichment_ER_pos)
#upsetplot_ER_pos


# TNBC
enrichment_TNBC <- enrichGO(gene = rownames(DEGs_TNBC), 
                              OrgDb = org.Hs.eg.db, # replace with the appropriate organism database
                              keyType = "SYMBOL", # change based on your gene ID type
                              ont = "BP", # Biological Process (BP), Molecular Function (MF), or Cellular Component (CC)
                              pAdjustMethod = "BH", # method for adjusting p values
                              qvalueCutoff = 0.05, # cutoff for q value
                              readable = TRUE) # convert entrez IDs to gene symbols

enrichment_pairwise_TNBC <- pairwise_termsim(enrichment_TNBC)

# Not sure now if this code is needed as I discovered there is GSEA in clusterProfiler
"""
colnames(enrichment_TNBC@result)

enrichment_TNBC@result$Description[1:10]

extracellular_matrix_organization_genes <- strsplit(enrichment_TNBC@result$geneID[1], "/")[[1]]
"""

# Barplot
barplot_TNBC <- barplot(enrichment_TNBC, showCategory=10, font.size=8)
barplot_TNBC

# Dotplot
#dotplot_TNBC <- dotplot(enrichment_TNBC, showCategory = 15) 
#dotplot_TNBC

# Treeplot
#treeplot_TNBC <- treeplot(enrichment_pairwise_TNBC)
#treeplot_TNBC

# Enrichment map
#emapplot_TNBC <- emapplot(enrichment_pairwise_TNBC)
#emapplot_TNBC

# Cnetplot
cnetplot_TNBC <- cnetplot(enrichment_TNBC,showCategory = 5, cex_label_gene=0.5, node_label='category')
cnetplot_TNBC

# Upsetplot
#upsetplot_TNBC <- upsetplot(enrichment_TNBC)
#upsetplot_TNBC

ontology_types <- c("BP", "MF", "CC")

perform_enrichment_analysis(ontology_types, DEGs_TNBC)

# BC_cell_line
#enrichment_BC_cell_line <- enrichGO(gene = rownames(DEGs_BC_cell_line), 
#                            OrgDb = org.Hs.eg.db, # replace with the appropriate organism database
#                            keyType = "SYMBOL", # change based on your gene ID type
#                            ont = "BP", # Biological Process (BP), Molecular Function (MF), or Cellular Component (CC)
#                            pAdjustMethod = "BH", # method for adjusting p values
#                            qvalueCutoff = 0.05, # cutoff for q value
#                            readable = TRUE) # convert entrez IDs to gene symbols

#enrichment_pairwise_BC_cell_line <- pairwise_termsim(enrichment_BC_cell_line)

# Barplot
#barplot_BC_cell_line <- barplot(enrichment_BC_cell_line, showCategory=20)
#barplot_BC_cell_line

# Dotplot
#dotplot_BC_cell_line <- dotplot(enrichment_BC_cell_line, showCategory = 15) 
#dotplot_BC_cell_line

# Treeplot
#treeplot_BC_cell_line <- treeplot(enrichment_pairwise_BC_cell_line)
#treeplot_BC_cell_line

# Enrichment map
#emapplot_BC_cell_line <- emapplot(enrichment_pairwise_BC_cell_line,cex_label_category = 0.5)
#emapplot_BC_cell_line

# Cnetplot
#cnetplot_BC_cell_line <- cnetplot(enrichment_BC_cell_line, cex_label_gene=0.5)
#cnetplot_BC_cell_line

# Upsetplot
#upsetplot_BC_cell_line <- upsetplot(enrichment_BC_cell_line)
#upsetplot_BC_cell_line

#rownames(DEGs_ER_pos)
#subset(dataset, rowData(dataset)$gene_type == 'protein_coding') %>%
#  .[rownames(.) %in% rownames(DEGs_ER_pos),]

#subset_rse <- rse[rownames(rse) %in% selected_genes, ]

#rownames(DEGs_ER_pos)
#mapping <- data.frame(gene_id = rowData(dataset)$gene_id, gene_name = rowData(dataset)$gene_name)
#result <- merge(data.frame(gene_name = rownames(DEGs_ER_pos)), mapping, by = "gene_name", all.x = TRUE)
#result$gene_id

#data_modified <- sapply(strsplit(result$gene_id,"\\."), function(x) x[1])

# Translating genes to proteins using BiomaRt
#ensembl <- useEnsembl(biomart = "genes")
#ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
#protein_DEGs <- getBM(attributes = c("ensembl_peptide_id", "ensembl_gene_id"),
#                      filters = "ensembl_gene_id",
#                      values = data_modified,
#                      mart = ensembl)

# Making Protein-Protein Interation (PPI) network using STRING
#string_db <- STRINGdb$new(version = "11.5", species = 96a06, score_threshold=200, input_directory="")
#example1_mapped <- string_db$map(protein_DEGs, "ensembl_peptide_id", removeUnmappedRows = TRUE )
#ppi_data <- string_db$plot_network(example1_mapped$STRING_id[1:2000])


# Write the DEGs to a CSV file
write.csv(DEGs_TNBC, file = "DEGs_TNBC_log2FoldChange_1_.csv", row.names = TRUE)
setwd('/Users/matviimykhailichenko/Library/CloudStorage/OneDrive-UniwersytetWrocławski/Carrier/Articles/Networks_of_Biomarkers_in_BC/Code/Networks-in-BC/Materials')
DEGs_TNBC %>%
  as.data.frame() %>% 
  .[, "log2FoldChange", drop = FALSE] %>%
  write.csv(., file = "DEGs_TNBC.csv", row.names = TRUE)

# Set working directory
setwd("/Users/matviimykhailichenko/Library/CloudStorage/OneDrive-UniwersytetWrocławski/Carrier/Articles/Networks_of_Biomarkers_in_BC/Code/Networks-in-BC/Files")

# Load necessary packages
library(readxl)
library(dplyr)

# Read the Excel file
data <- read_excel("PPI_network_MCODE_clusters.xls", sheet = "Sheet1") %>%
  na.omit()

# Ensure log2FoldChange is numeric
data$log2FoldChange <- as.numeric(data$log2FoldChange)

# Filter data for each cluster
cluster1 <- data %>% filter(MCODE_clusters == "Cluster 1")
cluster2 <- data %>% filter(MCODE_clusters == "Cluster 2")
cluster3 <- data %>% filter(MCODE_clusters == "Cluster 3")
cluster4 <- data %>% filter(MCODE_clusters == "Cluster 4")

# Determine the maximum length
max_length <- max(nrow(cluster1), nrow(cluster2), nrow(cluster3), nrow(cluster4))

# Function to extend the vector length with NA
extend_vector <- function(cluster) {
  extended_vector <- rep(NA, max_length)
  extended_vector[1:length(cluster$log2FoldChange)] <- cluster$log2FoldChange
  return(extended_vector)
}

# Extend the length of log2FoldChange vectors
cluster1_lf <- extend_vector(cluster1)
cluster2_lf <- extend_vector(cluster2)
cluster3_lf <- extend_vector(cluster3)
cluster4_lf <- extend_vector(cluster4)

# Calculate Pearson correlation
correlation_1_2 <- cor(cluster1_lf, cluster2_lf, use = "complete.obs")
correlation_1_3 <- cor(cluster1_lf, cluster3_lf, use = "complete.obs")
correlation_1_4 <- cor(cluster1_lf, cluster4_lf, use = "complete.obs")
correlation_2_3 <- cor(cluster2_lf, cluster3_lf, use = "complete.obs")
correlation_2_4 <- cor(cluster2_lf, cluster4_lf, use = "complete.obs")
correlation_3_4 <- cor(cluster3_lf, cluster4_lf, use = "complete.obs")

# Display the correlations
print(list("Correlation between Cluster 1 and 2" = correlation_1_2,
           "Correlation between Cluster 1 and 3" = correlation_1_3,
           "Correlation between Cluster 1 and 4" = correlation_1_4,
           "Correlation between Cluster 2 and 3" = correlation_2_3,
           "Correlation between Cluster 2 and 4" = correlation_2_4,
           "Correlation between Cluster 3 and 4" = correlation_3_4))

# Example of testing the significance of the correlation between Cluster 1 and 2
test_result <- cor.test(cluster3_lf, cluster4_lf, method = "pearson")

pmcplot(query, period, proportion = TRUE)