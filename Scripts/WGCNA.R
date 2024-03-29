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
)

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
