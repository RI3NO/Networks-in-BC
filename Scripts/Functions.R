retrieve_dds <- function(project_id) {
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
  
  
  colData <- colData(rse_data_expanded_sample_attributes)
  
  dds <- DESeqDataSetFromMatrix(countData = assay(rse_data_expanded_sample_attributes),
                                     colData = colData,
                                     rowData = rowData(rse_data_expanded_sample_attributes),
                                     design = ~ condition)
  
  return(dds)
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
      cluster_rows = TRUE,
      column_labels = colnames(DEGs_scaled_counts),
      name = "Z-score",
      cluster_columns = TRUE,
      top_annotation = .,
      show_column_names = FALSE
    )
  
  return(h1)
}

perform_enrichment_analysis <- function(ontology_types, DEGs_TNBC) {
  setwd("/Users/matviimykhailichenko/Library/CloudStorage/OneDrive-UniwersytetWrocławski/Carrier/Articles/Networks_of_Biomarkers_in_BC/Code/Networks-in-BC/Plots/Enrichment_Analysis")
  version <- 1
  while (file.exists(paste("v", version, sep = ""))) {
    version <- version + 1  # Increment version if file already exists
  }
  dir_name <- paste("v",version, sep='')
  dir.create(dir_name)  # Create the versioned directory
  setwd(paste(getwd(), "/", "v",version, sep=''))
  for (ont_type in ontology_types) {
    enrichment_TNBC <- enrichGO(
      gene = rownames(DEGs_TNBC),
      OrgDb = org.Hs.eg.db, # Replace with the appropriate organism database
      keyType = "SYMBOL",   # Change based on your gene ID type
      ont = ont_type,       # Ontology type (e.g., "BP", "MF", "CC")
      pAdjustMethod = "BH", # Method for adjusting p values
      qvalueCutoff = 0.05,  # Cutoff for q value
      readable = TRUE       # Convert entrez IDs to gene symbols
    )
    
    enrichment_pairwise_TNBC <- pairwise_termsim(enrichment_TNBC)
    
    # Barplot
    barplot_TNBC <- barplot(enrichment_TNBC, showCategory = 20, font.size = 8)
    
    # Dotplot
    dotplot_TNBC <- dotplot(enrichment_TNBC, showCategory = 15)
    
    # Treeplot
    treeplot_TNBC <- treeplot(enrichment_pairwise_TNBC)
    
    # Enrichment map
    emapplot_TNBC <- emapplot(enrichment_pairwise_TNBC)
    
    # Cnetplot
    cnetplot_TNBC <- cnetplot(enrichment_TNBC,showCategory = 5, cex_label_gene=0.5, node_label='category')
    
    # Upsetplot
    upsetplot_TNBC <- upsetplot(enrichment_TNBC)
    # Create and save plots with unique file names including version number
    for (plot_name in c("barplot", "dotplot", "treeplot", "emapplot", "cnetplot", "upsetplot")) {
      file_name <- paste(plot_name,"_" , ont_type, "_TNBC_v", version, ".svg", sep='')
      svg(file_name)
      print(get(paste(plot_name, "TNBC", sep = "_")))
      dev.off()
    }
  }
  setwd("/Users/matviimykhailichenko/Library/CloudStorage/OneDrive-UniwersytetWrocławski/Carrier/Articles/Networks_of_Biomarkers_in_BC/Code/Networks-in-BC")
}

