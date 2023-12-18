
# Transpose the assay data
assay_data_transposed <- t(assay(dds))

# Convert it into a data frame
assay_data_df <- as.data.frame(assay_data_transposed)

# Assuming your data matrix is called 'expression_data'
tsne_result <- Rtsne(assay_data_df, dims = 2, perplexity = 30, verbose = TRUE)

# Assuming 'dds' is your DESeqDataSet object

# Extract the assay data as a dataframe
assay_data <- as.data.frame(assay(dds))

# Identify and remove duplicate rows
non_duplicates <- !duplicated(assay_data)
filtered_assay_data <- assay_data[non_duplicates, ]

# Convert it back to a DESeqDataSet object if needed
filtered_dds <- DESeqDataSetFromMatrix(countData = filtered_assay_data, colData = colData(dds), design = ~condition)

results <- as.data.frame(tsne_result$Y)

# Assuming 'results_df' is your data frame with hyperparameter results
ggplot(results, aes(x = V1, y = V2)) +
  geom_point()
