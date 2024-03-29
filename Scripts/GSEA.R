# Require dds_TNBC, results_TNBC and packages that are imported in Main.R

library(pathview)

organism <- "org.Hs.eg.db"

DEGs_TNBC <- subset(results_TNBC, padj < 0.05) %>%
  as.data.frame()

DEGs_TNBC$gene_name <- rownames(DEGs_TNBC)

gene_mapping <- data.frame(gene_id=rowData(dds_TNBC)$gene_id, gene_name=rowData(dds_TNBC)$gene_name)

DEGs_TNBC <- merge(DEGs_TNBC, gene_mapping[,c("gene_name", "gene_id")], by="gene_name")

DEGs_TNBC$gene_id <- sub("\\..*", "", DEGs_TNBC$gene_id)

original_gene_list <- DEGs_TNBC$log2FoldChange

names(original_gene_list) <- DEGs_TNBC$gene_id

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list <- sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")
gse_pairwise <- pairwise_termsim(enrichment_TNBC)
require(DOSE)
#Showing 10 categories as more is unreadable
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
#Showing 30 categories as it looks good
emapplot(gse_pairwise, showCategory = 30)
# categorySize can be either 'pvalue' or 'geneNum'. Nice thing, but should be polished
# as it looks ugly now
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory=30, 
         cex_label_gene=0, cex_label_category=0.5, cex_category=0.5,
         cex_gene=0.5, layout = "kk")

# No clue what it shows or could it be useful. The distribution is quite simialr
ridgeplot(gse, core_enrichment = TRUE, label_format=40, orderBy = "NES",
          decreasing = TRUE, showCategory = 10) + labs(x = "enrichment distribution")

# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)

terms <- gse$Description[1:20]
terms_with_TNBC <- paste0(terms, " TNBC")
# Nice graph, maybe will show what terms are interesting 
pmcplot(terms_with_TNBC, 2010:2023, proportion=TRUE)


