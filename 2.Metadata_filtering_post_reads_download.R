library(tidyverse)
library(openxlsx)
library(here)
library(DESeq2)
library(sva)

species <- c("Arabidopsis thaliana", "Triticum aestivum", "Solanum lycopersicum")
for(sp in species) {
  coldata <- read.xlsx(here(paste0("Project/", sp, "/Data/metadata_filtered.xlsx")))
  counts <- read.table(here(paste0("Project/", sp, "/Data/compiled_counts_table.txt")), header = T)
  coldata <- coldata %>%
    filter(Run %in% colnames(counts))
  
  general_metadata <- read.xlsx(here("Project/general_metadata.xlsx"))
  
  control_numbers <- coldata %>%
    group_by(BioProject, Type.of.sample) %>%
    summarise(number_control = n()) %>% ungroup() %>% filter(Type.of.sample == "CONTROL")
  
  stressed_numbers <- coldata %>%
    group_by(BioProject, Type.of.sample) %>%
    summarise(number_stressed = n()) %>% ungroup() %>% filter(Type.of.sample == "STRESSED")
  
  treatment_numbers <- control_numbers %>%
    inner_join(stressed_numbers, by = "BioProject") %>% select(!c(Type.of.sample.x, Type.of.sample.y))
  
  # removing studies on non-drought stress
  coldata <- coldata %>%
    left_join(treatment_numbers) %>%
    filter(number_control >= 1, number_stressed >= 1) %>%
    inner_join(general_metadata, by = join_by(BioProject == `Repository.Identifier.BioProject.(SRA,.ArrayExpress)`), keep = FALSE) %>% distinct(Run, .keep_all = TRUE)  %>%
    select(!c(Species, Study.short.summary, Stress.type, number_sample, `repeats.(biological.replicates)`, total_sample, `Platform(s)`, Publication, Pubmed.ID, DOI, `(article.title)`, Repository.Identifier.as.indicated.in.the.original.publication, Major_Category, number_control, number_stressed))
  
  write.xlsx(coldata, here(paste0("Project/", sp, "/Data/metadata_filtered_2.xlsx")))
  
}

