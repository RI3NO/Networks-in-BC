

setwd("/Users/matviimykhailichenko/Documents/GitHub/Networks-in-BC/Scripts")

library(recount3)
library(DESeq2)
library(tidyverse)
library(STRINGdb)
library(igraph)
library(biomaRt)
library(org.Hs.eg.db)


# Loading dataset from recount3
human_projects <- available_projects()
proj_info <- subset(
  human_projects,
  project == "SRP042620" & project_type == "data_sources"
)
SRP042620 <- create_rse(proj_info)

colData(SRP042620) %>%
  as.data.frame() %>%
  .$sra.sample_attributes



