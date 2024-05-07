# Dataset is SRP042620


library(recount3)
library(DESeq2)
library(tidyverse)
library(STRINGdb)
library(igraph)
library(biomaRt)
library(org.Hs.eg.db)
BiocManager::install("CorLevelPlot")
library(CorLevelPlot)
library(gridExtra)
library(WGCNA)
library(recount3)
library(DESeq2)

getwd()


setwd(file.path(getwd(),"Scripts"))

source("Functions.R")
getwd()


# Retrieving Dataset "SRP042620" using recount 3 with designed condition (normal, cancer) in form of DeSeq2 object
dds_SRP042620 <- retrieve_dds("SRP042620")
View(dds_SRP042620)
#dds_SRP042620@colData
#colData(dds_SRP042620)

# Subsetting it to 20 sample because it's computationally heavy to take more(even 20), but 20 is enough for nice WGCNA
#subset_dds <- dds_SRP042620[,1:20]

# The below is demonstrated the way how to get expression matrix from DeSeq2 object
# assay(dds_SRP042620)


# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes
gsg <- goodSamplesGenes(t(assay(dds_SRP042620)))
summary(gsg)
gsg$allOK
# if false, there are outliers in genes or in samples
table(gsg$goodGenes)
# there are 6072 genes that are outliers
table(gsg$goodSamples)
# all samles are nice, no outliers

# Exclude Outliers 
dds_SRP042620 <- dds_SRP042620[gsg$goodGenes == TRUE,]
#dds_SRP042620
# Detect outliers samples
htree <- hclust(dist(t(assay(dds_SRP042620)), method = "manhattan"))

plot(htree)



# pca - method 2
pca <- prcomp(t(assay(dds_SRP042620)))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)
    
# remove.packages("gtable")
# install.packages("gtable")

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))


### NOTE: If there are batch effects observed, correct for them before moving ahead

# exclude outlier samples
# samples.to.be.excluded <- c('')
# data.subset <- dds_SRP042620[,!(colnames(assay(dds_SRP042620)) %in% samples.to.be.excluded)]


## Remove all genes with counts < 15 in more than 75% of samples (168*0.75=126)
## Suggested by WGCNA on RNAseq FAQ

# dds_SRP042620 <- dds_SRP042620[rowSums(counts(dds_SRP042620) >= 15) >= 126,]
# rowSums(counts(dds_SRP042620))
# nrow(dds_SRP042620) # 20244 genes
# dds_SRP042620

# without removing
subset_dds <- dds_SRP042620

# Doing DeSeq function (finding DEGs) on DDS object
deseq_result <- DESeq(subset_dds)
#View(deseq_result)

# getting results of DeSeq
res <- results(deseq_result)
res
#head(assay(deseq_result))

vsd <- varianceStabilizingTransformation(deseq_result)
#head(assay(vsd))

wpn_vsd <- getVarianceStabilizedData(deseq_result)
#wpn_vsd

rv_wpn <- rowVars(wpn_vsd)
#rv_wpn
summary(rv_wpn)

library(genefilter)

q75_wpn <- quantile( rowVars(wpn_vsd), .75)  # <= original
q75_wpn
q95_wpn <- quantile( rowVars(wpn_vsd), .95)  # <= changed to 95 quantile to reduce dataset
q95_wpn
expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]
expr_normalized[1:5,1:10]

names50samples <- colnames(expr_normalized[,1:50])

dim(expr_normalized)

expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(
    Gene_id = row.names(expr_normalized)
  ) %>%
  pivot_longer(-Gene_id)


expr_normalized_df
names50samples
expr_normalized_df[expr_normalized_df$name %in% names50samples,] %>% ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "Normalized and 95 quantile Expression",
    x = "treatment",
    y = "normalized expression"
  )

# Preparing and transposing out expr matrix for WGCNA
input_mat = t(expr_normalized)
# We can see now that the rows = treatments and columns = gene probes
input_mat[1:5,1:10]           # Look at first 5 rows and 10 columns
nrow(input_mat)
ncol(input_mat)

######  WGCNA !!!
#library(WGCNA)
allowWGCNAThreads()          # allow multi-threading (optional)
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)
sft.data <- sft$fitIndices
sft.data
par(mfrow = c(1,2));
cex1 = 0.9;

# Plotting Scale Independence and Mean Connectivity  
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

# Second Way

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

# picked_power is optional depending on plots 
picked_power = 14
temp_cor <- cor       

cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

View(netwk)

cor <- temp_cor     # Return cor function to original namespace

# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- netwk$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
table(netwk$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(netwk$dendrograms[[1]], cbind(netwk$unmergedColors, netwk$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# grey module = all genes that doesn't fall into other modules were assigned to the grey module




# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)

# Plot the dendrogram and the module colors underneath
# Plot for merged modules
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

#netwk$colors[netwk$blockGenes[[1]]]
#table(netwk$colors)

unique(colData(deseq_result)$condition)


module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

module_df[1:5,]

table(module_df$colors)


write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes
MEs0


# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")



# 6B. Intramodular analysis: Identifying driver genes ---------------

# Calculate the module membership and the associated p-values
# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.
module.membership.measure <- cor(module_eigengenes, input_mat, use = 'p')
nSamples <- nrow(input_mat)
nGenes <- ncol(input_mat)
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples) # Calculates Student asymptotic p-value for given correlations.
write.xlsx(module.membership.measure.pvals, file.path("C:/Users/RTIntelektFBT/Desktop/Roman_Project/Networks-in-BC/Networks-in-BC/Files/module.membership.measure.pvals.xlsx"))

View(module.membership.measure.pvals)
# Using the module membership measures you can identify genes with high module membership in interesting modules.


ncol(input_mat)
length(netwk$colors)
top.hub_genes <- chooseTopHubInEachModule(
      datExpr = input_mat,     # Gene expression data with rows as samples and columns as genes
      colorh = netwk$colors,   # The module assignments (color vectors) corresponding to the rows in datExpr
      omitColors = "grey",
      power = 2,
      type = "signed", 
      )
top.hub_genes
  



# pick out a few modules of interest here
# OPTIONAL< DEPENDING ON OUR RESULTS !!!!
modules_of_interest = c("turquoise", "blue")
#modules_of_interest = c("turquoise")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes
expr_normalized[1:5,1:10]

subexpr = expr_normalized[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")



### GSEA with clusterprofiler
library(clusterProfiler)

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors == "turquoise")

row.names(module_df) = module_df$gene_id

module_df

# Get normalized expression for those genes
expr_normalized[1:5,1:10]
expr_normalized

# subset genes from expression matrix that are in our module 
subexpr = expr_normalized[submod$gene_id,]

# res
# View(res)
res_1module <- res[submod$gene_id,]
res_1module
view(res_1module)


df_2_5_6 <- res_1module[,c(2,5,6)]
df_2_5_6$Symbol <- rownames(df_2_5_6)
df_2_5_6

df_1_4_7_8 <- df_2_5_6[,c(4,1,2,3)]
rownames(df_1_4_7_8) <- NULL
df_1_4_7_8



rankings <- sign(df_1_4_7_8$log2FoldChange)*(-log10(df_1_4_7_8$pvalue)) # we will use the signed p values from spatial DGE as ranking
names(rankings) <- df_1_4_7_8$Symbol # genes as names
rankings
#View(rankings)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
head(rankings)

plot(rankings)



#rankings_gseGO <- df_1_4_7_8$log2FoldChange
#names(rankings_gseGO) <- df_1_4_7_8$Symbol
#rankings_gseGO <- na.omit(rankings_gseGO)
#rankings_gseGO = sort(rankings_gseGO, decreasing = TRUE)
#rankings_gseGO
organism <- "org.Hs.eg.db"
gse <- gseGO(geneList = rankings, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

View(gse)
View(gse@result)

getwd()
gse_result <- arrange(gse@result, p.adjust, NES)
gse_result
library(xlsx)
write.xlsx(gse_result, file.path("C:/Users/RTIntelektFBT/Desktop/Roman_Project/Networks-in-BC/Networks-in-BC/Files/WGCNA_GSEA_Turq_Module.xlsx"))
# gse_pairwise <- pairwise_termsim(enrichment_TNBC)
require(DOSE)
library(enrichplot)
enriched_terms <- pairwise_termsim(gse)

plot_path <-"C:/Users/RTIntelektFBT/Desktop/Roman_Project/Networks-in-BC/Networks-in-BC/Plots"


#Showing 10 categories as more is unreadable
Dotlot <- dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
Dotlot
ggsave(file = file.path(plot_path,"WGCNA_GSEA_Turq_Module_Dotlot.svg"), plot = Dotlot, width=10, height=8)
dev.off()

#Showing 30 categories as it looks good
Emapplot <- emapplot(enriched_terms, showCategory = 20)
Emapplot
ggsave(file = file.path(plot_path,"WGCNA_GSEA_Turq_Module_Emapplot.svg"), plot = Emapplot, width=10, height=8)


# categorySize can be either 'pvalue' or 'geneNum'. Nice thing, but should be polished
# as it looks ugly now
Cnetplot <- cnetplot(gse, categorySize="pvalue", foldChange=rankings, showCategory = 20, 
         cex_label_gene=0, cex_label_category=0.5, cex_category=0.5,
         cex_gene=0.5, layout = "kk")
Cnetplot
ggsave(file = file.path(plot_path,"WGCNA_GSEA_Turq_Module_Cnetplot.svg"), plot = Emapplot, width=10, height=8)

# No clue what it shows or could it be useful. The distribution is quite simialr
# A Ridgeplot can be used to visualize these enrichment scores across multiple gene sets 
# and conditions simultaneously. Each ridge in the plot represents a gene set, and the width 
# and height of the ridge represent the density of enrichment scores for that gene set across 
# different conditions. This visualization allows for the comparison of enrichment patterns across 
# conditions and helps identify which gene sets are consistently enriched or depleted across 
# experimental groups.

Ridgeplot <- ridgeplot(gse, core_enrichment = TRUE, label_format=40, orderBy = "NES",
          decreasing = TRUE, showCategory = 10) + labs(x = "enrichment distribution")
Ridgeplot
ggsave(file = file.path(plot_path,"WGCNA_GSEA_Turq_Module_Ridgeplot.svg"), plot = Ridgeplot, width=10, height=8)

# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
Gseaplot <- gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)
Gseaplot
ggsave(file = file.path(plot_path,"WGCNA_GSEA_Turq_Module_Gseaplot.svg"), plot = Gseaplot, width=10, height=8)


terms <- gse$Description[1:20]
terms_with_TNBC <- paste0(terms, " TNBC")
# Nice graph, maybe will show what terms are interesting 
Pmcplot <- pmcplot(terms_with_TNBC, 2010:2023, proportion=TRUE)
Pmcplot
ggsave(file = file.path(plot_path,"WGCNA_GSEA_Turq_Module_Pmcplot.svg"), plot = Pmcplot, width=10, height=8)






#
###### GSEA  with fgsea
library("fgsea")
library(RColorBrewer) # for a colourful plot
# Functions ===================================================
## Function: Adjacency matrix to list -------------------------
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

## Function: prepare_gmt --------------------------------------
prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE){
  # for debug
  #file <- gmt_files[1]
  #genes_in_data <- df$gene_symbol
  
  # Read in gmt file
  gmt <- gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  #Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(genes_in_data, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again
  final_list <- matrix_to_list(mat) # for this we use the function we previously defined
  
  if(savefile){
    saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  
  print('Wohoo! .gmt conversion successfull!:)')
  return(final_list)
}


# pick out a few modules of interest here
# OPTIONAL< DEPENDING ON OUR RESULTS !!!!
#modules_of_interest = c("turquoise", "blue")

modules_of_interest = c("turquoise")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors == "turquoise")

row.names(module_df) = module_df$gene_id

module_df

# Get normalized expression for those genes
expr_normalized[1:5,1:10]
expr_normalized

subexpr = expr_normalized[submod$gene_id,]

# res
# View(res)
res_1module <- res[submod$gene_id,]
res_1module
view(res_1module)




gmt_file <- "BP_subset_GO_Hs.symbols.gmt"


my_genes <- rownames(res_1module)

bg_genes <- prepare_gmt(gmt_file, my_genes, savefile = FALSE)
View(bg_genes)

# df_1_4_7_8 <- data.frame( Symbol = character(), 
#                           Log2FC = numeric(), 
#                           P_value = numeric(), 
#                           Adj.P_value = numeric())
# View(df_1_4_7_8)
# df_1_4_7_8$Symbol <- rownames(res_1module)
# rm(df_1_4_7_8)
# 

#rm(df_1_4_7_8)
df_2_5_6 <- res_1module[,c(2,5,6)]
df_2_5_6$Symbol <- rownames(df_2_5_6)
df_2_5_6

df_1_4_7_8 <- df_2_5_6[,c(4,1,2,3)]
rownames(df_1_4_7_8) <- NULL
df_1_4_7_8

# Prepare your ranked list of genes
rankings <- sign(df_1_4_7_8$log2FoldChange)*(-log10(df_1_4_7_8$pvalue)) # we will use the signed p values from spatial DGE as ranking
names(rankings) <- df_1_4_7_8$Symbol # genes as names
#View(rankings)
#names(rankings)


rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
head(rankings)
plot(rankings)


ggplot(data.frame(gene_symbol = names(rankings)[1:50], ranks = rankings[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## 4. Run GSEA ---------------------------------------------------------------
# Easy peasy! Run fgsea with the pathways 
GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

head(GSEAres)
View(GSEAres)
GSEAres_ordered <- GSEAres[order(padj), ]
View(GSEAres_ordered)


ggplot(GSEAres_ordered[1:15,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj < 0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()


library(pathview)
require(DOSE)

BiocManager::install("clusterProfiler", version = "3.8")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)


#Showing 10 categories as more is unreadable
dotplot(GSEAres_ordered, showCategory=10, split=".sign") + facet_grid(.~.sign)


## 6. Check results ------------------------------------------------------
# Top 6 enriched pathways (ordered by p-val)
head(GSEAres[order(pval), ])

sum(GSEAres[, padj < 0.01])
sum(GSEAres[, pval < 0.01])


number_of_top_pathways_up = 10
number_of_top_pathways_down = 10
topPathwaysUp <- GSEAres[ES > 0][head(order(padj), n = number_of_top_pathways_up), pathway]
topPathwaysDown <- GSEAres[ES < 0][head(order(padj), n = number_of_top_pathways_down), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

#pdf(file = paste0(filename, '_gsea_top30pathways.pdf'), width = 20, height = 15)
plotGseaTable(bg_genes[topPathways], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)
#dev.off()

# Select only independent pathways, removing redundancies/similar pathways
collapsedPathways <- collapsePathways(GSEAres[order(padj)][padj < 0.05], bg_genes, rankings)
mainPathways <- GSEAres[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
#pdf(file = paste0('GSEA/Selected_pathways/', paste0(filename, background_genes, '_gsea_mainpathways.pdf')), width = 20, height = 15)
plotGseaTable(bg_genes[mainPathways], rankings, GSEAres, gseaParam = 0.5)
#dev.off()

#If you’d like to export the tables, just uncomment the 2 lines above. You can also export to .png, or other formats:
#png(file = paste0(filename, ‘_gsea_mainpathways.png’), width = 1500, height = 800)
#plotGseaTable(bg_genes[mainPathways], rankings, GSEAres, gseaParam = 0.5)





# plot the most significantly enriched pathway
plotEnrichment(bg_genes[[head(GSEAres[order(padj), ], 1)$pathway]],
               rankings) + 
  labs(title = head(GSEAres[order(padj), ], 1)$pathway)


bg_genes[[head(GSEAres[order(padj), ], 1)$pathway]]

#View(GSEAres_ordered)
#head(GSEAres_ordered,1)
#head(GSEAres[order(padj), ], 2)$pathway
#GSEAres_ordered[2,]

bg_genes[[GSEAres_ordered[2,]$pathway]]
View(bg_genes)
# plot the 2 most significantly enriched pathway
plotEnrichment(bg_genes[[GSEAres_ordered[2,]$pathway]],
               rankings) + 
  labs(title = GSEAres_ordered[2,]$pathway)

# plot the 2 most significantly enriched pathway

plotGseaTable(
  pathways = bg_genes[GSEAres_ordered$pathway[1:5]],
  stats = rankings,
  fgseaRes = GSEAres_ordered,
  gseaParam = 1,
  colwidths = c(5, 3, 0.8, 1.2, 1.2),
  pathwayLabelStyle = NULL,
  headerLabelStyle = NULL,
  valueStyle = NULL,
  axisLabelStyle = NULL,
  render = NULL
)





