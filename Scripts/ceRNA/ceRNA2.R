library(miRBaseConverter)
result3$Mature1

lncrna <- sub("\\..*", "", lncrna)
proteing_coding
proteing_coding <- sub("\\..*", "", proteing_coding)
# Extract entries with "-5p" from the vector
demiRNA_names <- result3$Mature1[grep("-5p", result3$Mature1, ignore.case = TRUE)] %>%
  unique()

result3$Mature1

library(GDCRNATools)


samples <- c('TCGA-2F-A9KO-01', 'TCGA-2F-A9KP-01', 
             'TCGA-2F-A9KQ-01', 'TCGA-2F-A9KR-01', 
             'TCGA-2F-A9KT-01', 'TCGA-2F-A9KW-01')
             
num_rows <- length(demiRNAs)
num_columns <- length(samples)

mir.expr <- as.data.frame(matrix(runif(num_rows * num_columns, min = -1, max = 1), nrow = num_rows, ncol = num_columns))

rownames(mir.expr) <- demiRNA_names
colnames(mir.expr) <- samples

num_rows <- length(de_proteing_coding_names + length(delncRNA_names)
num_columns <- length(samples)

# Create a random dataframe with specified dimensions
rna.expr <- as.data.frame(matrix(runif(num_rows * num_columns, min = -1, max = 1), nrow = num_rows, ncol = num_columns))

rownames(expression) <- c(de_proteing_coding_names, delncRNA_names)
colnames(expression) <- samples

ceOutput <- gdcCEAnalysis(lnc       = delncRNA_names, 
                          pc          = de_proteing_coding_names, 
                          deMIR = rownames(demiRNA_names),
                          lnc.targets = 'starBase', 
                          pc.targets  = 'starBase', 
                          rna.expr    = rna.expr, 
                          mir.expr    = mir.expr)
