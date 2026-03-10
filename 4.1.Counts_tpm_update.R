library(tidyverse)
library(openxlsx)
library(here)

species <- c("Arabidopsis thaliana", "Triticum aestivum", "Solanum lycopersicum")
for(sp in species) {
  coldata <- read.xlsx(here(paste0("Prova/", sp, "/Data/metadata_filtered_2.xlsx")))
  counts <- read.table(here(paste0("Prova/", sp, "/Data/counts_gene_level.txt")))
  tpm <- read.table(here(paste0("Prova/", sp, "/Data/tpm_gene_level.txt")))
  
  counts2 <- read.delim(here(paste0("Prova/", sp, "/Data/compiled_counts_table.txt")))
  
  tpm2 <- read.delim(here(paste0("Prova/", sp, "/Data/compiled_tpm_table.txt"))) 
  
  metadata <- coldata %>%
    filter(Run %in% colnames(counts2))
  
  
  # transform transcript level counts to gene level counts
  counts_gene_level <-
    counts2 %>% 
    select(target_id, metadata$Run) %>%
    mutate(target_id = str_to_upper(target_id)) %>%
    mutate(GENEID = str_remove(target_id, "\\..")) %>%
    relocate(target_id, GENEID) %>%
    group_by(GENEID) %>%
    summarise(across(all_of(metadata$Run), sum, .names = "{.col}")) %>%
    column_to_rownames(var = "GENEID")
  
  counts_tot <- add_column(counts, counts_gene_level) 
  coldata <- coldata %>%
    filter(Run %in% colnames(counts_tot))
  counts_tot <- counts_tot %>%
    select(coldata$Run)
  
  tpm_gene_level <-
    tpm2 %>% 
    select(target_id, metadata$Run) %>%
    mutate(target_id = str_to_upper(target_id)) %>%
    mutate(GENEID = str_remove(target_id, "\\..")) %>%
    relocate(target_id, GENEID) %>%
    group_by(GENEID) %>%
    summarise(across(all_of(metadata$Run), sum, .names = "{.col}")) %>%
    column_to_rownames(var = "GENEID")
  
  tpm_tot <- add_column(tpm, tpm_gene_level) %>%
    select(coldata$Run)
  
  write.table(counts_tot, here(paste0("Prova/", sp, "/Data/counts_gene_level.txt")))
  write.table(tpm_tot, here(paste0("Prova/", sp, "/Data/tpm_gene_level.txt")))
  write.xlsx(coldata, here(paste0("Prova/", sp, "/Data/metadata_filtered_2.xlsx")))
  
}
