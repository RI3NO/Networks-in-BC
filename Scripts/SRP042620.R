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

DEmRNAs <- subset(merged_df, gene_type == "protein_coding")

DEmiRNAs <- subset(mapped_DEGs, gene_type == "miRNA")

DElncRNAs <- subset(mapped_DEGs, gene_type == c("lincRNA","lincRNA"))


upregulated_DEmiRNAs_counts <- counts(dds)[rownames(upregulated_DEmiRNAs), ]

denoised_DEmRNAs <- subset(DEmRNAs, baseMean>1000 & log2FoldChange > 2)

denoised_DEmRNAs <- denoised_DEmRNAs[order(denoised_DEmRNAs$log2FoldChange, decreasing = TRUE),]

mat<-assay(vsd)[rownames(denoised_DEmRNAs), rownames(colData(dds))] #sig genes x samples

base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat.scaled)<-colnames(mat)

num_keep <- 25
#1 to num_keep len-num_keep to len
rows_keep <- c(seq(1:num_keep), seq((nrow(mat.scaled)-num_keep), nrow(mat.scaled)) )

l2_val <- as.matrix(denoised_DEmRNAs[rows_keep,]$log2FoldChange) #getting log2 value for each gene we are keeping
colnames(l2_val)<-"logFC"

mean <- as.matrix(denoised_DEmRNAs[rows_keep,]$baseMean) #getting mean value for each gene we are keeping
colnames(mean)<-"AveExpr"

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#maps values between b/w/r for min and max l2 values
col_logFC <- colorRamp2(c(min(l2_val),0, max(l2_val)), c("blue", "white", "red")) 

#maps between 0% quantile, and 75% quantile of mean values --- 0, 25, 50, 75, 100
col_AveExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red"))

ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2), 
                                               height = unit(2, "cm")))

h1 <- Heatmap(mat.scaled[rows_keep,], cluster_rows = F, 
              column_labels = colnames(mat.scaled), name="Z-score",
              cluster_columns = T)
h2 <- Heatmap(l2_val, row_labels = denoised_DEmRNAs$gene_name[rows_keep], 
              cluster_rows = F, name="logFC", top_annotation = ha, col = col_logFC,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(l2_val[i, j],2), x, y)
              })
h3 <- Heatmap(mean, row_labels = denoised_upregulated_DEmRNAs$gene_name[rows_keep], 
              cluster_rows = F, name = "AveExpr", col=col_AveExpr,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(mean[i, j],2), x, y)
              })

h<-h1+h2+h3
h
