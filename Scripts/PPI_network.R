# Translating genes to proteins using BiomaRt
ensembl <- useEnsembl(biomart = "genes")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
protein_DEGs <- getBM(attributes = c("ensembl_peptide_id"),
                      filters = "ensembl_gene_id_version",
                      values = rownames(DEGs),
                      mart = ensembl)

# Making Protein-Protein Interation (PPI) network using STRING
string_db <- STRINGdb$new(version = "11.0", species = 9606, score_threshold=200, input_directory="")
example1_mapped <- string_db$map(protein_DEGs, "ensembl_peptide_id", removeUnmappedRows = TRUE )
ppi_data <- string_db$plot_network(example1_mapped)
