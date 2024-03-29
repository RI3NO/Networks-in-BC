library(GDCRNATools)
library(egdeR)

# Trying treshold of l2fc > 2 to see if ceRNA looks good 
DEGs_TNBC <- subset(results_TNBC, padj < 0.05 & abs(log2FoldChange) > 2)

DEGs_dds_TNBC <- dds_TNBC[rownames(dds_TNBC) %in% rownames(DEGs_TNBC),]

desired_gene_types <- c("protein_coding", "lincRNA")
filtered_DEGs_dds_TNBC <- DEGs_dds_TNBC[rowData(DEGs_dds_TNBC)$gene_type %in% desired_gene_types, ]
filtered_DEGs_dds_TNBC_voom_normalized <- gdcVoomNormalization(assay(filtered_DEGs_dds_TNBC)) %>%
  as.data.frame()

rownames(DEGs_dds_TNBC) <- rowData(DEGs_dds_TNBC)$gene_id

proteing_coding <- rowData(DEGs_dds_TNBC)[rowData(DEGs_dds_TNBC)$gene_type == "protein_coding",] %>%
  rownames()

lncrna <- rowData(DEGs_dds_TNBC)[rowData(DEGs_dds_TNBC)$gene_type == "lincRNA",] %>%
  rownames() 

dds_miRNA <- DEGs_dds_TNBC[rowData(DEGs_dds_TNBC)$gene_type == "miRNA",]
rowData(dds_miRNA)

rownames(dds_miRNA) <- rowData(dds_miRNA)$gene_name

rowData(dds_miRNA)$gene_id

# Assuming rowData(dds_miRNA)$gene_id is a character vector
rowData(dds_miRNA)$gene_id <- sub("\\..*", "", rowData(dds_miRNA)$gene_id)

miRNA <- rowData(dds_miRNA)$gene_id 

miRNA_dds_TNBC_voom_normalized <- gdcVoomNormalization(assay(dds_miRNA)) %>%
  as.data.frame()

DEGs_dds_TNBC_voom_normalized <- gdcVoomNormalization(assay(DEGs_dds_TNBC)) %>%
  as.data.frame()

library(biomaRt)

# Connect to the Ensembl database using biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

rowData(dds_miRNA)$gene_name

# Get miRNA information for the specified Ensembl ID
miRNA_BM <- getBM(
  attributes = c("mirbase_accession", "ensembl_gene_id"),
  filters = "ensembl_gene_id",
  values = miRNA,
  mart = ensembl
)
miRNA_info <- miRNA_BM$mirbase_accession %>%
  as.character()
listAttributes(ensembl)
listFilters(ensembl)

library(miRBaseConverter)


#### Step 3. miRBase Accessions to miRNA Names of the target version
result2 = miRNA_AccessionToName(miRNA_info,targetVersion = "v22")
result2
result3 <- miRNA_PrecursorToMature(result2$TargetName)
miRNA_BM

merged <- merge(result3, result2, by.x = "OriginalName", by.y = "Accession", all.x = TRUE)

# Assuming you have two DataFrames: result2 and miRNA_BM

# Use the merge() function to perform the mapping
merged_data <- merge(result3, miRNA_BM, by.x = "Accession", by.y = "mirbase_accession", all.x = TRUE)

# The resulting merged_data DataFrame will have both Accession and ensembl_gene_id columns

rowData(dds_miRNA)

# Assuming you have the gene_id and merged_data DataFrames

# Merge by "gene_id" column from gene_id and "ensembl_gene_id" column from merged_data
merged_result <- merge(rowData(dds_miRNA), merged_data, by.x = "gene_id", by.y = "ensembl_gene_id", all.x = TRUE)

# The resulting merged_result DataFrame will contain the merged data
merged_result$gene_id



duplicates <- duplicated(merged_result$gene_id)

# Filter out rows with duplicates in "gene_id" and update merged_result
merged_result <- merged_result[!duplicates, ]

rowData(dds_miRNA) <- merged_result

rownames(dds_miRNA) <- rowData(dds_miRNA)$TargetName

miRNA <- rownames(dds_miRNA)

miRNA_dds_TNBC_voom_normalized <- gdcVoomNormalization(assay(dds_miRNA)) %>%
  as.data.frame()


gdcCEAnalysis(lnc = lncrna, pc = proteing_coding, deMIR = miRNA, lnc.targets = "starBase",
              pc.targets = "starBase", rna.expr = filtered_DEGs_dds_TNBC_voom_normalized, mir.expr = miRNA_dds_TNBC_voom_normalized)

gdcCEAnalysis(lnc = lncrna, pc = proteing_coding, deMIR = NULL, lnc.targets = "starBase",
              pc.targets = "starBase", rna.expr = filtered_DEGs_dds_TNBC_voom_normalized, mir.expr = miRNA_dds_TNBC_voom_normalized)
