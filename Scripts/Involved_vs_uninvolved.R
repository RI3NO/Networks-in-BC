# Create a vector of tissue types you want to keep
tissue_types_to_keep <- c("ER+ Breast Cancer Primary Tumor", "Uninvolved Breast Tissue Adjacent to ER+ Primary Tumor")

# Subset the cells based on the tissue types
ER_pos_rse <- SRP042620_expanded_sample_attributes[, colData(SRP042620_expanded_sample_attributes)$sra_attribute.tissue %in% tissue_types_to_keep]

ER_pos_rse %>%
  colData() %>%
  .$sra_attribute.tissue

ER_pos_colData <- colData(subset_rse) %>% 
  as.data.frame() %>%
  .$condition %>%
  data.frame(condition=.)

ER_pos_dds <- DESeqDataSetFromMatrix(countData = assay(ER_pos_rse),
                                     colData = ER_pos_colData,
                                     design = ~ condition)

#Performing variance stabilizing transformation (vst)
ER_pos_vsd <- vst(ER_pos_dds)

# Perform PCA
PCA_plot <- plotPCA(ER_pos_vsd, intgroup = "condition")
PCA_plot

SRP042620_expanded_sample_attributes %>%
  colData() %>%
  .$sra_attribute.tissue

# Create a vector of tissue types you want to keep
tissue_types_to_keep <- c("Uninvolved Breast Tissue Adjacent to TNBC Primary Tumor", "Triple Negative Breast Cancer Primary Tumor")

# Subset the cells based on the tissue types
TBNC_rse <- SRP042620_expanded_sample_attributes[, colData(SRP042620_expanded_sample_attributes)$sra_attribute.tissue %in% tissue_types_to_keep]

TBNC_rse %>%
  colData() %>%
  .$sra_attribute.tissue

TBNC_colData <- colData(TBNC_rse) %>% 
  as.data.frame() %>%
  .$condition %>%
  data.frame(condition=.)

TBNC_dds <- DESeqDataSetFromMatrix(countData = assay(TBNC_rse),
                                     colData = TBNC_colData,
                                     design = ~ condition)

#Performing variance stabilizing transformation (vst)
TBNC_vsd <- vst(TBNC_dds)

# Perform PCA
PCA_plot <- plotPCA(TBNC_vsd, intgroup = "condition")
PCA_plot

SRP042620_expanded_sample_attributes %>%
  colData() %>%
  .$sra_attribute.tissue

# Create a vector of tissue types you want to keep
tissue_types_to_keep <- c("Breast Cancer Cell Line", "Reduction Mammoplasty - No known cancer")

# Subset the cells based on the tissue types
BC_cell_line_rse <- SRP042620_expanded_sample_attributes[, colData(SRP042620_expanded_sample_attributes)$sra_attribute.tissue %in% tissue_types_to_keep]

BC_cell_line_rse %>%
  colData() %>%
  .$sra_attribute.tissue

BC_cell_line_colData <- colData(BC_cell_line_rse) %>% 
  as.data.frame() %>%
  .$condition %>%
  data.frame(condition=.)

BC_cell_line_dds <- DESeqDataSetFromMatrix(countData = assay(BC_cell_line_rse),
                                   colData = BC_cell_line_colData,
                                   design = ~ condition)

#Performing variance stabilizing transformation (vst)
BC_cell_line_vsd <- vst(BC_cell_line_dds)

# Perform PCA
PCA_plot <- plotPCA(BC_cell_line_vsd, intgroup = "condition")
PCA_plot

ER_pos_dds <- DESeq(ER_pos_dds)

ER_pos_results <- results(ER_pos_dds, contrast = c('condition','cancer','normal'))

ER_pos_DEGs <- subset(ER_pos_results, padj < 0.05 & abs(log2FoldChange) > 2 & baseMean > 1000 )

ER_pos_DEGs <- ER_pos_DEGs[order(ER_pos_DEGs$log2FoldChange, decreasing = TRUE),]

DEGs_subset_ER_pos <- ER_pos_dds[rownames(ER_pos_dds) %in% rownames(ER_pos_DEGs),]

mat<-assay(DEGs_subset_ER_pos)[rownames(DEGs_subset_ER_pos), rownames(colData(DEGs_subset_ER_pos))] #sig genes x samples

base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat.scaled) <- colnames(mat)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

gene_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
                   filters = "ensembl_gene_id_version", 
                   values = rownames(mat.scaled), 
                   mart = ensembl)
gene_info

listAttributes(ensembl)



rownames(mat.scaled) <- gene_info$external_gene_name
library(biomaRt)

num_keep <- 25

rows_keep <- c(seq(1:num_keep), seq((nrow(mat.scaled)-num_keep), nrow(mat.scaled)) )

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)


heatmap_annotation <- HeatmapAnnotation(
  Condition = column_colors$condition,  # Use the 'condition' column as annotation data
  col = list(Condition = c("cancer" = "red", "normal" = "blue"))  # Define colors for annotation levels
)

h1 <- Heatmap(mat.scaled[rows_keep,], cluster_rows = F, 
              column_labels = colnames(mat.scaled), name="Z-score",
              cluster_columns = T,
              top_annotation = heatmap_annotation,
              show_column_names = FALSE)

h1

column_colors <- ER_pos_dds %>%
  colData() %>%
  .$condition %>%
  data.frame(condition=.)

rownames(column_colors) <- ER_pos_dds %>%
  colnames()

# Create a color mapping for the conditions
color_mapping <- c("cancer" = "1", "normal" = "-1")

# Assuming your dataframe is called 'df'
# Add a new column 'color' based on the 'condition' column
column_colors <- column_colors %>%
  mutate(color = color_mapping[condition])



col_fun = colorRamp2(c(-1, 1), c("blue", "red"))


column_colors <- column_colors[, 2, drop = TRUE] %>%
  as.numeric()

column_colors <- col_fun(column_colors)
