library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

# Create a vector of tissue types you want to keep
tissue_types_to_keep <- c("ER+ Breast Cancer Primary Tumor", "Uninvolved Breast Tissue Adjacent to ER+ Primary Tumor")

# Subset the cells based on the tissue types
rse_ER_pos <- SRP042620_expanded_sample_attributes[, colData(SRP042620_expanded_sample_attributes)$sra_attribute.tissue %in% tissue_types_to_keep]

rse_ER_pos %>%
  colData() %>%
  .$sra_attribute.tissue

colData_ER_pos <- colData(subset_rse) %>% 
  as.data.frame() %>%
  .$condition %>%
  data.frame(condition=.)

dds_ER_pos <- DESeqDataSetFromMatrix(countData = assay(rse_ER_pos),
                                     colData = colData_ER_pos,
                                     design = ~ condition)

#Performing variance stabilizing transformation (vst)
vsd_ER_pos <- vst(dds_ER_pos)

# Perform PCA
PCA_plot_ER_pos <- plotPCA(vsd_ER_pos, intgroup = "condition")
PCA_plot_ER_pos

SRP042620_expanded_sample_attributes %>%
  colData() %>%
  .$sra_attribute.tissue

# Create a vector of tissue types you want to keep
tissue_types_to_keep <- c("Uninvolved Breast Tissue Adjacent to TNBC Primary Tumor", "Triple Negative Breast Cancer Primary Tumor")

# Subset the cells based on the tissue types
rse_TBNC <- SRP042620_expanded_sample_attributes[, colData(SRP042620_expanded_sample_attributes)$sra_attribute.tissue %in% tissue_types_to_keep]

rse_TBNC %>%
  colData() %>%
  .$sra_attribute.tissue

colData_TBNC <- colData(rse_TBNC) %>% 
  as.data.frame() %>%
  .$condition %>%
  data.frame(condition=.)

dds_TBNC <- DESeqDataSetFromMatrix(countData = assay(rse_TBNC),
                                     colData = colData_TBNC,
                                     design = ~ condition)

#Performing variance stabilizing transformation (vst)
vsd_TBNC <- vst(dds_TBNC)

# Perform PCA
PCA_plot_TBNC <- plotPCA(vsd_TBNC, intgroup = "condition")
PCA_plot_TBNC

SRP042620_expanded_sample_attributes %>%
  colData() %>%
  .$sra_attribute.tissue

# Create a vector of tissue types you want to keep
tissue_types_to_keep <- c("Breast Cancer Cell Line", "Reduction Mammoplasty - No known cancer")

# Subset the cells based on the tissue types
rse_BC_cell_line <- SRP042620_expanded_sample_attributes[, colData(SRP042620_expanded_sample_attributes)$sra_attribute.tissue %in% tissue_types_to_keep]

rse_BC_cell_line %>%
  colData() %>%
  .$sra_attribute.tissue

colData_BC_cell_line <- colData(rse_BC_cell_line) %>% 
  as.data.frame() %>%
  .$condition %>%
  data.frame(condition=.)

dds_BC_cell_line <- DESeqDataSetFromMatrix(countData = assay(rse_BC_cell_line),
                                   colData = colData_BC_cell_line,
                                   design = ~ condition)

#Performing variance stabilizing transformation (vst)
vsd_BC_cell_line <- vst(dds_BC_cell_line)

# Perform PCA
PCA_plot_BC_cell_line <- plotPCA(vsd_BC_cell_line, intgroup = "condition")
PCA_plot_BC_cell_line

# ER_pos DEGs Analysis

dds_ER_pos <- DESeq(dds_ER_pos)

results_ER_pos <- results(dds_ER_pos, contrast = c('condition','cancer','normal'))

DEGs_ER_pos <- subset(results_ER_pos, padj < 0.05 & abs(log2FoldChange) > 2 & baseMean > 1000 )

DEGs_ER_pos <- DEGs_ER_pos[order(DEGs_ER_pos$log2FoldChange, decreasing = TRUE),]

DEGs_dds_ER_pos <- dds_ER_pos[rownames(dds_ER_pos) %in% rownames(DEGs_ER_pos),]

DEGs_counts_ER_pos<-assay(DEGs_dds_ER_pos)[rownames(DEGs_dds_ER_pos), rownames(colData(DEGs_dds_ER_pos))] #sig genes x samples


# Function to generate heatmap, taking a DEGs_counts matrix and a dds object
generate_heatmap <- function(DEGs_counts, dds) {
  DEGs_scaled_counts<- t(apply(DEGs_counts, 1, scale)) #center and scale each column (Z-score) then transpose
  colnames(DEGs_scaled_counts) <- colnames(DEGs_counts)
  num_keep <- 25
  rows_keep <- c(seq(1:num_keep), seq((nrow(DEGs_scaled_counts) - num_keep + 1), nrow(DEGs_scaled_counts)))
  
  condition_vector <- dds %>%
    colData() %>%
    .$condition
  
  h1 <- HeatmapAnnotation(
    Condition = condition_vector,  # Use the 'condition' column as annotation data
    col = list(Condition = c("cancer" = "red", "normal" = "blue"))  # Define colors for annotation levels
  ) %>%
    Heatmap(
      DEGs_scaled_counts[rows_keep, ],
      cluster_rows = FALSE,
      column_labels = colnames(DEGs_scaled_counts),
      name = "Z-score",
      cluster_columns = TRUE,
      top_annotation = .,
      show_column_names = FALSE
    )
  
  return(h1)
}

Heatmap_ER_pos <- generate_heatmap(DEGs_counts_ER_pos, ER_pos_dds)
Heatmap_ER_pos

# TBNC DEGs Analysis

dds_TBNC <- DESeq(dds_TBNC)
results_TBNC <- results(dds_TBNC, contrast = c('condition','cancer','normal'))
DEGs_TBNC <- subset(results_TBNC, padj < 0.05 & abs(log2FoldChange) > 2 & baseMean > 1000 )
DEGs_TBNC <- DEGs_TBNC[order(DEGs_TBNC$log2FoldChange, decreasing = TRUE),]
DEGs_dds_TBNC <- dds_TBNC[rownames(dds_TBNC) %in% rownames(DEGs_TBNC),]
DEGs_counts_TBNC <-assay(DEGs_dds_TBNC)[rownames(DEGs_dds_TBNC), rownames(colData(DEGs_dds_TBNC))] #sig genes x samples

Heatmap_TBNC <- generate_heatmap(DEGs_counts_TBNC, dds_TBNC)
Heatmap_TBNC



# Code to map ENSEMBL ids
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

gene_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
                   filters = "ensembl_gene_id_version", 
                   values = rownames(mat.scaled), 
                   mart = ensembl)
gene_info

