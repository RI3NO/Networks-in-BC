library(recount3)
library(DESeq2)

 
setwd(file.path(getwd(),"Scripts"))

source("Functions.R")
dds_SRP042620 <- retrieve_dds("SRP042620")
dds_SRP042620
View(assay(dds_SRP042620))

deseq_result <- DESeq(dds_SRP042620)
deseq_result
res <- results(deseq_result)
View(res)
#!any(rownames(deseq_result) == rownames(rowData(dds_SRP042620)))
rowData(dds_SRP042620)

Symbol_ENSGid <- data.frame(SYMBOL = rownames(dds_SRP042620), 
                            gene_id = rowData(dds_SRP042620)$gene_id, 
                            gene_type = rowData(dds_SRP042620)$gene_type )

Symbol_ENSGid <- Symbol_ENSGid[Symbol_ENSGid$gene_type == "protein_coding",]

res_protein_coding <- res[Symbol_ENSGid$SYMBOL,]
res_protein_coding
res_protein_coding <- na.omit(res_protein_coding)
res_protein_coding <- res_protein_coding[res_protein_coding$padj <= 0.05,]
res_protein_coding



DEGs_Symbol_ENSGid_protein_coding <- Symbol_ENSGid[Symbol_ENSGid$SYMBOL %in% rownames(res_protein_coding),]
DEGs_Symbol_ENSGid_protein_coding


# Regular expression pattern to match up to the first dot
pattern <- "\\..*$"

# Use sub to extract the substring
DEGs_Symbol_ENSGid_protein_coding$sub_gene_id <- sub(pattern, "", DEGs_Symbol_ENSGid_protein_coding$gene_id)
sum(duplicated(DEGs_Symbol_ENSGid_protein_coding$sub_gene_id))

sum(duplicated(DEGs_Symbol_ENSGid_protein_coding$SYMBOL))

getwd()

list.files(file_path)

file_path <-  file.path("C:/Users/RTIntelektFBT/Desktop/Roman_Project/Networks-in-BC/Networks-in-BC/","Files")


# clusters_path <- file.path("C:/Users/RTIntelektFBT/Desktop/Roman_Project/Networks-in-BC/Networks-in-BC/","Files","clusters")

all_clusters <- list.files(clusters_path)
all_clusters

for (cluster in all_clusters) {
    cluster.vec <- as.vector(read.table(file.path(clusters_path,cluster)))
    print(cluster.vec)
    
    # create new directory for cluster to save the results there
    new_dir <- file.path(file_path,paste0(cluster))
    dir.create(file.path(new_dir))
    
    run_GSEA(cluster.vec, cluster, new_dir)
    #print(res[cluster.vec,])
    ENSG00000223972.5
    print(res[cluster.vec,])
    stop()

    
}

DEGs_Symbol_ENSGid_protein_coding
res_protein_coding


class(res_protein_coding)
library(clusterProfiler)
require(DOSE)
library(enrichplot)

run_GSEA <- function(cluster.vec, cluster, new_dir) {

    
    # ____ - it is something I will be using (ENSEBL or SYMBOL)
    #DEGs <- res_protein_coding[DEGs_Symbol_ENSGid_protein_coding$____ %in% cluster.vec,]
    res_protein_coding_cluster <-  res_protein_coding[cluster.vec,]

    # Prepare your ranked list of genes
    rankings <- res_protein_coding_cluster$log2FoldChange*(-log10(res_protein_coding_cluster$padj)) # we will use the signed p values from spatial DGE as ranking
    names(rankings) <- rownames(res_protein_coding_cluster) # genes as names
    #View(rankings)
    rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
    head(rankings)
    plot(rankings)
    rankings
    organism <- "org.Hs.eg.db"
    gse <- gseGO(geneList = rankings,  # rankings_gseGO
                 ont ="BP",  # ALL
                 keyType = "SYMBOL", 
                 nPerm = 10000, 
                 minGSSize = 20, 
                 maxGSSize = 1200, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = organism, 
                 pAdjustMethod = "BH" # none or BH
    ) 
    # In the context of Gene Set Enrichment Analysis (GSEA), the gene ratio refers to the ratio of the number of genes within a gene set that are differentially expressed or otherwise of interest in your dataset to the total number of genes in that gene set.
    
   # View(gse@result)
    
    write.xlsx(gse@result,
               file.path(new_dir,paste0("GSEA_",cluster)))
    
    # Get the similarity matrix. pairwise_termsim - this function add similarity matrix to the termsim slot of enrichment result
    enriched_terms <- pairwise_termsim(gse)
    View(enriched_terms@result)
    enriched_terms@termsim
    
    #Showing 10 categories as more is unreadable
    GSEA_Dotplot_TOP5 <- dotplot(gse, showCategory = 5, split=".sign",font.size = 10) + facet_grid(.~.sign)
    GSEA_Dotplot_TOP5
    ggsave(file = file.path(new_dir,paste0("GSEA_dotplot_",cluster,".svg"), plot = GSEA_Dotplot_TOP5, 
           width = 10, height = 8))
    
    
    #Showing 30 categories as it looks good
    # Enrichment Map for enrichment result of over-representation test or gene set enrichment analysis
    GSEA_Emapplot_Top5 <- emapplot(enriched_terms, showCategory = 5, cex_label_category = 0.8)
    GSEA_Emapplot_Top5
    ggsave(file = file.path(new_dir,paste0("GSEA_Emapplot_",cluster,".svg")), plot = GSEA_Emapplot_Top5, 
           width = 10, height = 8)
    
    
    
    
    # creare gene list to show LogFC for genes in cnetplot
    original_gene_list <- res_protein_coding_cluster$log2FoldChange
    names(original_gene_list) <- rownames(res_protein_coding_cluster)
    # original_gene_list
    # omit any NA values 
    gene_list <- na.omit(original_gene_list)
    gene_list = sort(gene_list, decreasing = TRUE)
    #length(original_gene_list)
    #length(gene_list)
    GSEA_cnetplot_TOP5 <- cnetplot(gse, categorySize="pvalue", foldChange = gene_list, 
                                         showCategory = 5, # 20
                                         cex_label_category = 0.8,# size for pathways
                                         cex_label_gene = 0.7, # size for genes
                                         shadowtext = "none"                               # cex_gene=0.5
                                         # cex_label_gene=0, cex_category=0.5,
                                         # cex_gene=0.5, layout = "kk"
    )
    GSEA_cnetplot_TOP5
    ggsave(file = file.path(new_dir,paste0("GSEA_cnetplot_",cluster,".svg")), plot = GSEA_cnetplot_TOP5, 
           width = 10, height = 8)
    
    ridgeplot <- ridgeplot(gse, core_enrichment = TRUE, label_format = 40, orderBy = "NES",
                           decreasing = TRUE, showCategory = 10) + labs(x = "enrichment distribution")
    ridgeplot
    ggsave(file = file.path(new_dir,paste0("GSEA_ridgeplot_",cluster,".svg")), plot = ridgeplot, 
           width = 10, height = 8)
    
    
    
}



###
#matrix_norm_counts <- vst(dds_SRP042620)
#matrix_norm_counts
#assay(matrix_norm_counts)
assay(dds_SRP042620)

####
resultsNames(deseq_result)
# "condition_normal_vs_cancer"
####
# LFC shrinkage uses information from all genes to generate more accurate estimates. 
# Specifically, the distribution of LFC estimates for all genes is used (as a prior) to shrink the 
# LFC estimates of genes with little information or high dispersion toward more likely (lower) LFC estimates.
#res <- lfcShrink(dds, coef="condition_normal_vs_cancer", type="apeglm")
