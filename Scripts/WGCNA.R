install.packages("WGCNA")
BiocManager::install("impute")
BiocManager::install("preprocessCore")
# WGCNA is available on CRAN

library(WGCNA)


data <- readr::read_delim("GSE61333_ligule_count.txt",     # <= path to the data file
                          delim = "\t")

col_sel = names(data)[-1]      # Get all but first column name
mdata <- data %>%
  tidyr::pivot_longer(
    .,                        # The dot is the the input data, magrittr tutorial
    col = all_of(col_sel)
  ) %>%  
  mutate(
    group = gsub("-.*","", name) %>% gsub("[.].*","", .)   # Get the shorter treatment names
  )
p <- mdata %>%
  ggplot(., aes(x = name, y = value)) +             # x = treatment, y = RNA Seq count
  geom_violin() +                                   # violin plot, show distribution
  geom_point(alpha = 0.2) +                         # scatter plot
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)          # Rotate treatment text
  ) +
  labs(x = "Treatment Groups", y = "RNA Seq Counts") +
  facet_grid(cols = vars(group), drop = TRUE, scales = "free_x")      # Facet by hour


p

de_input = as.matrix(data[,-1])
row.names(de_input) = data$GeneId

meta_df <- data.frame( Sample = names(data[-1])) %>%
  mutate(
    Type = gsub("-.*","", Sample) %>% gsub("[.].*","", .)
  )

dds <- DESeqDataSetFromMatrix(round(de_input),
                              meta_df,
                              design = ~Type)

dds <- DESeq(dds)

vsd <- varianceStabilizingTransformation(dds)
library(genefilter)  

wpn_vsd <- getVarianceStabilizedData(dds)
rv_wpn <- rowVars(wpn_vsd)
summary(rv_wpn)

q75_wpn <- quantile( rowVars(wpn_vsd), .75)  # <= original
q95_wpn <- quantile( rowVars(wpn_vsd), .95)  # <= changed to 95 quantile to reduce dataset
expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]

expr_normalized[1:5,1:10]

expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(
    Gene_id = row.names(expr_normalized)
  ) %>%
  pivot_longer(-Gene_id)

expr_normalized_df %>% ggplot(., aes(x = name, y = value)) +
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

install.packages("WGCNA")
BiocManager::install("impute")
BiocManager::install("preprocessCore")
# WGCNA is available on CRAN

library(WGCNA)


counts_TNBC <- counts(dds_TNBC) %>%
  as_tibble(.,rownames="gene_name")

annotation <- tibble(sample_name=colData(dds_TNBC)$external_id, condition=colData(dds_TNBC)$condition)





col_sel <- names(counts_TNBC)[-1]      # Get all but first column name
mdata <- counts_TNBC %>%
  tidyr::pivot_longer(
    .,                        # The dot is the the input data, magrittr tutorial
    col = all_of(col_sel)
  )

# Assuming the second tibble is named 'annotation'
mdata <- left_join(mdata, annotation, by = c("name" = "sample_name"))

p <- mdata %>%
  ggplot(., aes(x = name, y = value)) +             # x = treatment, y = RNA Seq count
  geom_violin() +                                   # violin plot, show distribution
  geom_point(alpha = 0.2) +                         # scatter plot
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)          # Rotate treatment text
  ) +
  labs(x = "Treatment Groups", y = "RNA Seq Counts") +
  facet_grid(cols = vars(condition), drop = TRUE, scales = "free_x")      # Facet by hour


p

wde_input = as.matrix(data[,-1])
row.names(de_input) = data$GeneId

meta_df <- data.frame( Sample = names(data[-1])) %>%
  mutate(
    Type = gsub("-.*","", Sample) %>% gsub("[.].*","", .)
  )


DEGs_dds_TNBC <- dds_TNBC[rowData(dds_TNBC)$gene_name  %in%  rownames(DEGs_TNBC),]
DEGs_normCounts_TNBC <- counts(DEGs_dds_TNBC, normalized = TRUE)

input_mat = t(DEGs_normCounts_TNBC)

input_mat[1:5,1:10]

library(WGCNA)
allowWGCNAThreads() 

powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)

par(mfrow = c(1,2));
cex1 = 0.9;

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

picked_power = 9
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

cor <- temp_cor

mergedColors = labels2colors(netwk$colors)

plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

module_df[1:5,]

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

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
annotation

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

modules_of_interest = c("green", "turquoise", "tan")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes
DEGs_normCounts_TNBC[1:5,1:10]

subexpr = DEGs_normCounts_TNBC[submod$gene_id,]

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

