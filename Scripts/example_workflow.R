# Example workflow for PPI network making 
# Dataset is SRP144041. Finding DEG in high/low proliferation BC groups

setwd("/Users/matviimykhailichenko/Documents/GitHub/Networks-in-BC/Scripts")

library(recount3)
library(DESeq2)
library(tidyverse)
library(STRINGdb)
library(igraph)
library(biomaRt)
library(org.Hs.eg.db)


# Loading dataset from recount3
human_projects <- available_projects()
proj_info <- subset(
  human_projects,
  project == "SRP144041" & project_type == "data_sources"
)
SRP144041 <- create_rse(proj_info)

# Subset the SummarizedExperiment to keep only genes with non-zero counts
SRP144041_filtered <- SRP144041[rowSums(assay(SRP144041)) > 0, ]

# Function to make inputed dataframe look good  
StupidCharToDataFrame <- function(char_vector) {
  # Initialize an empty list to store the rows
  rows_list <- list()

  # Loop through each line and split based on "|"
  for (line in char_vector) {
    fields <- unlist(strsplit(line, "\\|"))
    row_data <- list()

    # Loop through each field and split based on ";;"
    for (field in fields) {
      parts <- strsplit(field, ";;")[[1]]
      if (length(parts) == 2) {
        row_data[[parts[1]]] <- parts[2]
      }
    }

    # Add the row data to the list
    rows_list <- c(rows_list, list(row_data))
  }

  # Convert the list of rows to a data frame
  
  df <- as.data.frame(do.call(rbind, rows_list))
  return(df)
}


df <- StupidCharToDataFrame(colData(SRP144041_filtered)$sra.sample_attributes) %>%
  DataFrame()

rownames(df) <- colnames(SRP144041_filtered)


colData(SRP144041_filtered) <- cbind(colData(SRP144041_filtered), df)

colData(SRP144041_filtered) <- DataFrame(lapply(colData(SRP144041_filtered), function(x) gsub(" ", "_", x)))


# Making DeSeqDataSet object from matrix
dds <- DESeqDataSetFromMatrix(countData = assays(SRP144041_filtered)$raw_counts,
                              colData = colData(SRP144041_filtered),
                              design = ~ proliferation.group)


# Subsetting it to 10 sample because it's computationally heavy
subset_dds <- dds[,1:10]

# Doing DeSeq function (finding DEGs) on DDS object
deseq_result <- DESeq(subset_dds)

# Transforming to do PCA
DESeqTransform_obj <- DESeqTransform(SRP144041_filtered)

colnames(DESeqTransform_obj) <- colnames(SRP144041) 

# Perform PCA
plotPCA(DESeqTransform_obj, intgroup = "proliferation.group")

# Outputing results and filtering them 
results <- results(deseq_result, contrast = c('proliferation.group','high','low'))

alpha <- 0.5 # Adjust according to your significance threshold
log2FC_threshold <- 1 # Adjust according to your desired log2 fold change threshold

DEGs <- subset(results, padj < alpha & abs(log2FoldChange) > log2FC_threshold) %>%
  na.omit() %>%
  rownames()


results <- as.data.frame(results)

# Making Volcanoplot

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
                   values = DEGs,
                   mart = ensembl)


# Making Protein-Protein Interation (PPI) network using STRING
string_db <- STRINGdb$new(version = "11.0", species = 9606, score_threshold=200, input_directory="")
example1_mapped <- string_db$map(protein_DEGs, "ensembl_peptide_id", removeUnmappedRows = TRUE )
ppi_data <- string_db$plot_network(example1_mapped)








