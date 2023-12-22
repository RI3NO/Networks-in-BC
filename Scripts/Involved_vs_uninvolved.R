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
