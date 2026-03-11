library(tidyverse)
library(openxlsx)
library(here)

species <- c("Arabidopsis thaliana", "Triticum aestivum", "Solanum lycopersicum")
for(sp in species) {
  coldata <- read.xlsx(here(paste0("Prova/", sp, "/Data/metadata_filtered_2.xlsx")))
  
  counts <- read.delim(here(paste0("Prova/", sp, "/Data/compiled_counts_table.txt"))) 
  
  # transform transcript level counts to gene level counts
  counts_gene_level <-
    counts %>% 
    select(target_id, coldata$Run) %>%
    mutate(target_id = str_to_upper(target_id)) %>%
    mutate(GENEID = str_remove(target_id, "\\.[0-9]+$")) %>%
    relocate(target_id, GENEID) %>%
    group_by(GENEID) %>%
    summarise(across(coldata$Run, sum, .names = "{.col}")) %>%
    column_to_rownames(var = "GENEID")
  
  
  
  write.table(counts_gene_level, here(paste0("Prova/", sp, "/Data/counts_gene_level.txt")))
  
}

