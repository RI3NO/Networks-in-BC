# The dataset has lncRNAs, miRNAs and mRNAs

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

# Computing read counts from base-pair counts
read_counts <- compute_read_counts(SRP042620_expanded_sample_attributes)

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

colData <- colData(SRP042620_expanded_sample_attributes) %>% 
  as.data.frame() %>%
  .$condition %>%
  data.frame(condition=.)

rownames(colData) <- rownames(colData(SRP042620_expanded_sample_attributes))



 

# Making DeSeqDataSet object from matrix
dds <- DESeqDataSetFromMatrix(countData = read_counts,
                              colData = colData,
                              design = ~ condition)

# Doing DeSeq function (finding DEGs) on DDS object
dds <- DESeq(dds)

vsd <- vst(dds)

# Perform PCA
pca_result <- plotPCA(vsd, intgroup = "condition")

pca_result

# PCAplot is okay 

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

mapped_DEGs <- rowData(SRP042620)[DEGs_names, ]

# Merge data frames by row names
merged_df <- merge(mapped_DEGs, DEGs, by = 0)

rownames(merged_df) <- rownames(mapped_DEGs)

upregulated_DEmRNAs <- subset(merged_df, gene_type == "protein_coding" & log2FoldChange > 0)
downregulated_DEmRNAs <- subset(merged_df, gene_type == "protein_coding" & log2FoldChange < 0)

upregulated_DEmiRNAs <- subset(merged_df, gene_type == "miRNA" & log2FoldChange > 0)
downregulated_DEmiRNAs <- subset(merged_df, gene_type == "miRNA" & log2FoldChange < 0)

upregulated_DElncRNAs <- subset(merged_df, gene_type == c("macro_lncRNA","lincRNA") & log2FoldChange > 0)
downregulated_DElncRNAs <- subset(merged_df, gene_type == c("macro_lncRNA","lincRNA") & log2FoldChange < 0)



DEmiRNAs <- subset(mapped_DEGs, gene_type == "miRNA")

DElncRNAs <- subset(mapped_DEGs, gene_type == c("lincRNA","lincRNA"))


upregulated_DEmiRNAs_counts <- counts(dds)[rownames(upregulated_DEmiRNAs), ]

denoised_upregulated_DEmRNAs <- subset(upregulated_DEmRNAs, baseMean>50 & log2FoldChange > 2)

denoised_upregulated_DEmRNAs <- denoised_upregulated_DEmRNAs[order(denoised_upregulated_DEmRNAs$log2FoldChange, decreasing = TRUE),]

mat<-assay(vsd)[rownames(denoised_upregulated_DEmRNAs), rownames(colData(dds))] #sig genes x samples

base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat.scaled)<-colnames(mat)

# Keeping 25 most upregulated genes

l2_val <- as.matrix(denoised_upregulated_DEmRNAs[1:25,]$log2FoldChange) #getting log2 value for each gene we are keeping
colnames(l2_val)<-"logFC"

mean <- as.matrix(denoised_upregulated_DEmRNAs[1:25,]$baseMean) #getting mean value for each gene we are keeping
colnames(mean)<-"AveExpr"

