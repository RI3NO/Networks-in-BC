library(recount3)
library(DESeq2)


setwd(file.path(getwd(),"Scripts"))

source("Functions.R")

getwd()

# Retrieving Dataset "SRP042620" using recount 3 with designed condition (normal, cancer) in form of DeSeq2 object
dds_SRP042620 <- retrieve_dds("SRP042620")
dds_SRP042620

!any(rownames(dds_SRP042620) == rowData(dds_SRP042620)$gene_name)
!any(rownames(assay(dds_SRP042620)) == rowData(dds_SRP042620)$gene_name)

# View(dds_SRP042620)

length(rownames(dds_SRP042620))
length(unique(rownames(dds_SRP042620)))

rowData(dds_SRP042620)

assay(dds_SRP042620)

!any(rownames(rowData(dds_SRP042620)) == rowData(dds_SRP042620)$gene_name)

annotation <- rowData(dds_SRP042620)[,c("gene_id","gene_type","gene_name")]
annotation$gene_name_row <- rownames(annotation)
View(annotation)

!any(rownames(annotation) == annotation$gene_name)
sum(rownames(annotation) == annotation$gene_name)
length(annotation$gene_name)
length(unique(annotation$gene_name))

!any(rownames(annotation) == annotation$gene_id)
!any(rownames(annotation) == annotation$gene_name_row)
sum(rownames(annotation) == annotation$gene_name_row)
sum(rownames(annotation) == annotation$gene_id)

str(dds_SRP042620)
annotation
View(annotation)

unique(annotation$gene_type)


gene_names_All <- rownames(annotation)
gene_names_All
length(gene_names_All)
length(unique(gene_names_All))


ensemble_IDs <- annotation$gene_id
ensemble_IDs

length(ensemble_IDs)
length(unique(ensemble_IDs))


dds_SRP042620
row.names(dds_SRP042620[gene_names_All,]) <- ensemble_IDs
rownames(dds_SRP042620)
View(assay(dds_SRP042620))




# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes
library(WGCNA)
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

length(rownames(assay(dds_SRP042620)))
length(unique(rownames(assay(dds_SRP042620))))


View(colData(dds_SRP042620))
# annotation$
View(annotation)
rownames(dds_SRP042620)
# Exclude outliers from annotation
annotation <- annotation[annotation$gene_id %in% rownames(dds_SRP042620),]
View(annotation)

# colData(dds_SRP042620)$sra.library_selection

# res_NAomited <- na.omit(res)
# res_NAomited

BiocManager::install("GDCRNATools")
library(GDCRNATools)

deseq_result <- DESeq(dds_SRP042620)
assay(deseq_result)
assay(dds_SRP042620)

# !any()
# dataFrameProbe <- data.frame( a = c(1,3,5),b = c(2,4,6))
# row.names(dataFrameProbe) <- c("A","B","C")
# FrameProbe <- dataFrameProbe[c("C","A","B"),]
# !any(rownames(dataFrameProbe[rownames(FrameProbe),]) == rownames(FrameProbe))
# dataFrameProbe["A",]


res <- results(deseq_result)
rownames(res)

# !any(rownames(res) == rownames(deseq_result))
# !any(rownames(res) == rownames(dds_SRP042620))


# row.names(res[gene_names_All,]) <- ensemble_IDs
# res
# res[gene_names_All,]
# c("BX004987.1","AC145212.2 ","Y_RNA")
# res[c("BX004987.1","AC145212.2 ","Y_RNA"),]


annotation
protein.coding.genes_id <- annotation[annotation$gene_type == "protein_coding",]$gene_id
protein.coding.genes_id
is.vector(protein.coding.genes_id)
dim(protein.coding.genes_id)
# protein.coding.genes_names <- annotation[annotation$gene_type == "protein_coding",]$gene_name
# protein.coding.genes_names

linc.RNA.genes_id <- annotation[annotation$gene_type == "lincRNA",]$gene_id
linc.RNA.genes_id
is.vector(linc.RNA.genes_id)
# linc.RNA.genes_names <- annotation[annotation$gene_type == "lincRNA",]$gene_name
# linc.RNA.genes_names

miRNA.genes_id <- annotation[annotation$gene_type == "miRNA",]$gene_id
miRNA.genes_id
is.vector(miRNA.genes_id)
# miRNA.genes_name <- annotation[annotation$gene_type == "miRNA",]$gene_name
# miRNA.genes_name

View(res)
diff_expr_all_miRNA <- res[miRNA.genes_id,]
diff_expr_all_miRNA

View(diff_expr_all_miRNA@listData)
diff_expr_miRNA <- na.omit(diff_expr_all_miRNA)
diff_expr_miRNA

miRNA.diff_expr <- rownames(diff_expr_miRNA)
miRNA.diff_expr
length(miRNA.diff_expr)

# assay(dds_SRP042620)  
# annotation[annotation$gene_id == "ENSG00000277400.1",]

# dds_SRP042620["ENSG00000277400.1",]
# assay(dds_SRP042620[protein.coding.genes_id,])
# counts(dds_SRP042620)
# assay(dds_SRP042620)




####### Normalization of RNAseq data #######


# assay(dds_SRP042620[c(protein.coding.genes_id,linc.RNA.genes_id)])
# assay(deseq_result[c(protein.coding.genes_id,linc.RNA.genes_id)])

# counts(dds_SRP042620[c(protein.coding.genes_id,linc.RNA.genes_id)])
# counts(deseq_result[c(protein.coding.genes_id,linc.RNA.genes_id)])
# counts(dds_SRP042620[c("ENSG00000215717.5")])
# assay(deseq_result[c("ENSG00000215717.5")])

rnaExpr <- gdcVoomNormalization(counts = assay(dds_SRP042620[c(protein.coding.genes_id,linc.RNA.genes_id),]), filter = FALSE)
rnaExpr

is.matrix(rnaExpr)
dim(rnaExpr)

library(biomaRt)
####### Normalization of miRNAs data #######
# Connect to miRBase database in biomaRt

# Initializes a connection to the Ensembl database
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
listEnsembl()
attributes = listAttributes(ensembl)
attributes
# 
# 
# # Annotate genes from dds_SRP042620 
dds_SRP042620_miRNAs <- dds_SRP042620[miRNA.diff_expr,]
dds_SRP042620_miRNAs

biomart_list <- getBM(filter = "ensembl_gene_id_version",  # is a vector of filters that one wil use as input to the query.
                      attributes = c("mirbase_accession","mirbase_id","ensembl_gene_id_version"), # is a vector of attributes that one wants to retrieve (= the output of the query).
                      values = rownames(assay(dds_SRP042620_miRNAs)),  # values is a vector of values for the filters.
                      mart = ensembl) # is an object of class Mart, which is created by the useEnsembl() function.
biomart_list
# biomart_list_clean <- biomart_list[complete.cases(biomart_list), ]
biomart_list_clean <- biomart_list[!nchar(biomart_list$mirbase_id) < 3 , ]
biomart_list_clean
# nchar('fdfd')

dds_SRP042620_miRNAs_ANN <- dds_SRP042620_miRNAs[biomart_list_clean$ensembl_gene_id_version,]
dds_SRP042620_miRNAs
dds_SRP042620_miRNAs_ANN
View(assay(dds_SRP042620_miRNAs_ANN))

nrow(biomart_list)
nrow(biomart_list_clean)

length(biomart_list$mirbase_id)
nrow(na.omit(biomart_list))

nrow(dds_SRP042620_miRNAs)
# length(unique(biomart_list$hgnc_symbol))
c("mirbase_accession","mirbase_id","mirbase_trans_name")

length(biomart_list_clean$mirbase_id)

!any(rownames(dds_SRP042620_miRNAs_ANN) == biomart_list_clean$ensembl_gene_id_version)

rownames(dds_SRP042620_miRNAs_ANN) <- biomart_list_clean$mirbase_id 
dds_SRP042620_miRNAs_ANN
View(assay(dds_SRP042620_miRNAs_ANN))
rownames(dds_SRP042620_miRNAs_ANN)


mirExpr <- gdcVoomNormalization(counts = assay(dds_SRP042620_miRNAs_ANN), filter = FALSE)
is.matrix(mirExpr)
dim(mirExpr)
mirExpr
linc.RNA.genes_id
gdcCEAnalysis(lnc = linc.RNA.genes_id, pc = protein.coding.genes_id, deMIR = NULL, lnc.targets = "miRcode",
              pc.targets = "miRcode", rna.expr = rnaExpr , mir.expr = mirExpr )

linc.RNA.genes_id

# Arguments
# lnc a vector of Ensembl long non-coding gene ids
# pc a vector of Ensembl protein coding gene ids
# deMIR a vector of differentially expressed miRNAs. Default is NULL
# lnc.targets a character string specifying the database of miRNA-lncRNA interactions. Should
# be one of 'spongeScan', 'starBase', and 'miRcode'. Default is 'starBase'.
# Or a list of miRNA-lncRNA interactions generated by users
# pc.targets a character string specifying the database of miRNA-lncRNA interactions. Should
# be one of 'spongeScan', 'starBase', and 'miRcode'. Default is 'starBase'.
# Or a list of miRNA-lncRNA interactions generated by users
# rna.expr voom transformed gene expression data
# mir.expr voom transformed mature miRNA expression data









# # annotating genes
# library(org.Hs.eg.db)
# keytypes(org.Hs.eg.db)
# 
# 
# # Initializes a connection to the Ensembl database
# ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# #ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
# listEnsembl()
# attributes = listAttributes(ensembl)
attributes
# 
# 
# # Annotate genes from dds_SRP042620 
biomart_list <- getBM(filter = "ensembl_gene_id",  # is a vector of filters that one wil use as input to the query.
                      attributes = c("ensembl_gene_id", "hgnc_symbol", "transcript_biotype",  #
                                     "description"), # is a vector of attributes that one wants to retrieve (= the output of the query).
                      values = rownames(assay(dds_SRP042620)),  # values is a vector of values for the filters.
                      mart = ensembl) # is an object of class Mart, which is created by the useEnsembl() function.
biomart_list
# length(unique(biomart_list$hgnc_symbol))
# 
# 
# unique(biomart_list$transcript_biotype) #Виводить всі унікальні типи РНК, які є в зразках







# rownames(assay(dds_SRP042620))
# ### Annotation of genes for diff_table
# sel = AnnotationDbi::select(org.Hs.eg.db, rownames(assay(dds_SRP042620)) , 
#                             columns = c("ENTREZID","GENENAME","ENSEMBL"), keytype="SYMBOL")
# sel
row(sel)
nrow(na.omit(sel))
# 
# # Filtering out Genes, that were not annotated
# sel = sel[!is.na(sel$ENTREZID),]
# 







project <- 'TCGA-CHOL'
rnadir <- paste(project, 'RNAseq', sep='/')
mirdir <- paste(project, 'miRNAs', sep='/')

####### Download RNAseq data #######
gdcRNADownload(project.id     = 'TCGA-CHOL', 
               data.type      = 'RNAseq', 
               write.manifest = FALSE,
               method         = 'gdc-client',
               directory      = rnadir)

####### Download mature miRNA data #######
gdcRNADownload(project.id     = 'TCGA-CHOL', 
               data.type      = 'miRNAs', 
               write.manifest = FALSE,
               method         = 'gdc-client',
               directory      = mirdir)



####### Parse RNAseq metadata #######
metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-CHOL',
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)

metaMatrix.RNA
####### Filter duplicated samples in RNAseq metadata #######
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)

####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in RNAseq metadata #######
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)

####### Parse miRNAs metadata #######
metaMatrix.MIR <- gdcParseMetadata(project.id = 'TCGA-CHOL',
                                   data.type  = 'miRNAs', 
                                   write.meta = FALSE)

####### Filter duplicated samples in miRNAs metadata #######
metaMatrix.MIR <- gdcFilterDuplicate(metaMatrix.MIR)

####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in miRNAs metadata #######
metaMatrix.MIR <- gdcFilterSampleType(metaMatrix.MIR)







library(GDCRNATools)
deLNC <- c('ENSG00000260920.2','ENSG00000242125','ENSG00000261211')
dePC <- c('ENSG00000043355.1','ENSG00000109586','ENSG00000144355')
genes <- c(deLNC, dePC)
samples <- c('TCGA-2F-A9KO-01', 'TCGA-2F-A9KP-01', 
             'TCGA-2F-A9KQ-01', 'TCGA-2F-A9KR-01', 
             'TCGA-2F-A9KT-01', 'TCGA-2F-A9KW-01')
rnaExpr <- data.frame(matrix(c(2.7,7.0,4.9,6.9,4.6,2.5,
                               0.5,2.5,5.7,6.5,4.9,3.8,
                               2.1,2.9,5.9,5.7,4.5,3.5,
                               2.7,5.9,4.5,5.8,5.2,3.0,
                               2.5,2.2,5.3,4.4,4.4,2.9,
                               2.4,3.8,6.2,3.8,3.8,4.2),6,6), 
                      stringsAsFactors=FALSE)
rownames(rnaExpr) <- genes
colnames(rnaExpr) <- samples

mirExpr <- data.frame(matrix(c(7.7,7.4,7.9,8.9,8.6,9.5,
                               5.1,4.4,5.5,8.5,4.4,3.5,
                               4.9,5.5,6.9,6.1,5.5,4.1,
                               12.4,13.5,15.1,15.4,13.0,12.8,
                               2.5,2.2,5.3,4.4,4.4,2.9,
                               2.4,2.7,6.2,1.5,4.4,4.2),6,6),
                      stringsAsFactors=FALSE)
colnames(mirExpr) <- samples
rownames(mirExpr) <- c('hsa-miR-340-5p','hsa-miR-181b-5p',
                       'hsa-miR-181a-5p', 'hsa-miR-181c-5p',
                       'hsa-miR-199b-5p','hsa-miR-182-5p')

ceOutput <- gdcCEAnalysis(lnc       = deLNC, 
                          pc          = dePC, 
                          lnc.targets = 'starBase', 
                          pc.targets  = 'starBase', 
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)
ceOutput




getwd()
clusters_path <- file.path("C:/Users/RTIntelektFBT/Desktop/Roman_Project/Networks-in-BC/Networks-in-BC/","Files","clusters")

all_clusters <- list.files(clusters_path)
all_clusters


res

dds_SRP042620

# Regular expression pattern
pattern <- "\\..*2$"

# Count the number of strings that match the pattern
count <- sum(grepl(pattern, rownames(dds_SRP042620)))
count

rownames(dds_SRP042620)
unique(rownames(dds_SRP042620))

for (cluster in all_clusters) {
    cluster.vec <- as.vector(read.table(file.path(clusters_path,cluster)))
    print(cluster.vec)
    
    #print(res[cluster.vec,])
    ENSG00000223972.5
    print(res[cluster.vec,])
    stop()
}


strings <- c("BX004987.1", "BX004987.2", "BX004987.3", "BX004987.12", "BX004987.21", "BX004987.22")

# Regular expression to match strings with a "2" after the dot
pattern <- "\\d+\\.2\\b"

# Count how many strings match the pattern
count <- sum(grepl(pattern, strings))
count


strings <- c("BX004987.1", "BX004987.2", "BX004987.3", "BX004987.12", "BX004987.21", "BX004987.22")
library(tidyverse)
strings_detected <- str_detect(strings, pattern = ".2$")
strings_detected



dds_SRP042620["Metazoa_SRP.1",]
dds_SRP042620
rowData(dds_SRP042620)$gene_id
strings_detected <- str_detect(rownames(dds_SRP042620), pattern = ".2$")
strings_detected
sum(strings_detected)

rowData(dds_SRP042620)$gene_id
strings_detected <- str_detect(rowData(dds_SRP042620)$gene_id, pattern = ".5$")
strings_detected
rownames(dds_SRP042620[strings_detected,])

sum(strings_detected)

Symbol_ENSGid <- data.frame(SYMBOL = rownames(dds_SRP042620), 
                            gene_id = rowData(dds_SRP042620)$gene_id, 
                            gene_type =rowData(dds_SRP042620)$gene_type )

Symbol_ENSGid <- Symbol_ENSGid[Symbol_ENSGid$gene_type == "protein_coding",]
# Regular expression pattern to match up to the first dot
pattern <- "\\..*$"

# Use sub to extract the substring
Symbol_ENSGid$sub_gene_id <- sub(pattern, "", Symbol_ENSGid$gene_id)
sum(duplicated(Symbol_ENSGid$sub_gene_id))

Symbol_ENSGid[duplicated(Symbol_ENSGid$sub_gene_id),]



deseq_result <- DESeq(dds_SRP042620)
res <- results(deseq_result)
res
