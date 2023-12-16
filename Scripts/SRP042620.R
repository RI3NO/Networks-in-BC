setwd("/Users/matviimykhailichenko/Documents/GitHub/Networks-in-BC/Scripts")

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

dds_coldata <- colData(dds)$sra_attribute.tissue %>%
  data.frame(condition = colData(dds)$condition, tissue = ., run = colData(dds)$sra.run_alias)

dds_coldata %>% 
  class()

# Making DeSeqDataSet object from matrix
dds <- DESeqDataSetFromMatrix(countData = read_counts,
                              colData = dds_coldata,
                              design = ~ condition + tissue)

# Doing DeSeq function (finding DEGs) on DDS object
dds <- DESeq(dds)

vsd <- vst(dds)

# Perform PCA
pca_result <- plotPCA(vsd, intgroup = "condition")

pca_result

# Hypothesise that sra.sample_name, sra_attribute.tissue and sra.run_alias could also introduce variance, checking sra_attribute.tissue


results <- results(dds)


degQC(counts(dds), design[["condition"]], pvalue = res[["pvalue"]])


dds_coldata <- colData(dds)$condition

dds_coldata <- colData(dds)$sra_attribute.tissue %>%
  data.frame(condition = dds_coldata, tissue = ., run = colData(dds)$sra.run_alias)

dds_coldata

colData(dds)$sra.sample_title %>%
  unique()

results <- results(dds, contrast = c('condition','cancer','normal'))

alpha <- 0.05 # Adjust according to your significance threshold
log2FC_threshold <- 1 # Adjust according to your desired log2 fold change threshold

# Filter DEGs
DEGs <- subset(results, padj < alpha & abs(log2FoldChange) > log2FC_threshold)

#Volcanoplot

results <- as.data.frame(results)

ggplot(results, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05, "Significant", "Not Significant")), size = 2) +
  labs(x = "log2 Fold Change", y = "-log10(Adjusted p-value)") +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  ggtitle("Volcano Plot") +
  geom_point(aes(color = ifelse(abs(log2FoldChange) > 1, "Significant", "Not Significant")), size = 2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = -1, linetype = "dashed", color = "red")


# Translating genes to proteins using BiomaRt
ensembl <- useEnsembl(biomart = "genes")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
protein_DEGs <- getBM(attributes = c("ensembl_peptide_id"),
                      filters = "ensembl_gene_id_version",
                      values = rownames(DEGs),
                      mart = ensembl)


# Making Protein-Protein Interation (PPI) network using STRING
string_db <- STRINGdb$new(version = "11.0", species = 9606, score_threshold=200, input_directory="")
example1_mapped <- string_db$map(protein_DEGs, "ensembl_peptide_id", removeUnmappedRows = TRUE )
ppi_data <- string_db$plot_network(example1_mapped)

# Subset upregulated and downregulated DEGs
upregulated_DEGs <- subset(DEGs, log2FoldChange > 0)
downregulated_DEGs<- subset(DEGs, log2FoldChange < 0)

DEGs_names <- rownames(DEGs)
rowData(SRP042620) %>%
  .$gene_type

DEGs_names %>%
  class()

mapped_DEGs <- rowData(SRP042620)[DEGs_names, ]

mapped_DEGs$gene_type %>%
  unique()

DEmRNAs <- subset(mapped_DEGs, gene_type == "protein_coding")

DEmiRNAs <- subset(mapped_DEGs, gene_type == "miRNA")

DElncRNAs <- subset(mapped_DEGs, gene_type == c("lincRNA","lincRNA"))

# Install and load the "gplots" package (if not already installed)
install.packages("gplots")
library(gplots)

DEmRNAs_counts <- counts(dds)[rownames(DEmRNAs), ]

# Create the heatmap
heatmap.2(as.matrix(DEmRNAs_counts), 
          scale = "row",  # Scale rows (genes)
          dendrogram = "both",  # Add dendrograms for rows and columns
          Rowv = TRUE, Colv = TRUE,  # Reorder rows and columns based on clustering
          col = colorRampPalette(c("blue", "white", "red"))(50),  # Choose a color palette
          main = "DEmRNAs Heatmap")

