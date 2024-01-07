retrieve_dataset <- function(project_id) {
  # Set working directory (consider passing this as an argument or setting it outside the function)
  setwd("/Users/matviimykhailichenko/Library/CloudStorage/OneDrive-UniwersytetWrocławski/Carrier/Articles/Networks_of_Biomarkers_in_BC/Code/Networks-in-BC")
  
  # Loading dataset from recount3
  human_projects <- available_projects()
  proj_info <- subset(
    human_projects,
    project == project_id & project_type == "data_sources"
  )
  rse_data <- create_rse(proj_info)
  rownames(rse_data) <- rowData(rse_data)$gene_name
  rse_data_expanded_sample_attributes <- expand_sra_attributes(rse_data)
  
  # Computing read counts from base-pair counts
  assay(rse_data_expanded_sample_attributes) <- compute_read_counts(rse_data_expanded_sample_attributes)
  
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
  colData(rse_data_expanded_sample_attributes)$condition <- sapply(colData(rse_data_expanded_sample_attributes)$`sra_attribute.source_name`, determine_condition)
  
  return(rse_data_expanded_sample_attributes)
}

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
