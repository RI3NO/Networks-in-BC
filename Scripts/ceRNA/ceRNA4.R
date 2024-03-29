
DEGs_TNBC <- subset(results_TNBC, padj < 0.05 & abs(log2FoldChange) > 1)
DEGs_dds_TNBC <- dds_TNBC[rownames(dds_TNBC) %in% rownames(DEGs_TNBC),]
dds_demiRNA <- DEGs_dds_TNBC[rowData(DEGs_dds_TNBC)$gene_type == "miRNA",]
rownames(dds_demiRNA)

getwd()
# Write the character vector to a CSV file
write.csv(rownames(dds_demiRNA), file = "demiRNA_names.csv", row.names = FALSE)

# Trying to use Targetscan to predict miRNA targets

library(targetscan.Hs.eg.db)
ls("package:targetscan.Hs.eg.db")

merged_result$Accession

fams <- sample(ls(targetscan.Hs.egMIRNA), 3)
mget("hsa-miR-6803-5p", targetscan.Hs.egMIRNA)



