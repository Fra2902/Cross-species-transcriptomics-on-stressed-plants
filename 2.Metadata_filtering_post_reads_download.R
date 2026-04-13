library(tidyverse)
library(openxlsx)
library(here)
library(DESeq2)
library(sva)

# List of species to process
species <- c("Arabidopsis thaliana", "Triticum aestivum", "Solanum lycopersicum")
for(sp in species) {
  
  # --- Load filtered metadata and count matrix ---
  coldata <- read.xlsx(here(paste0("Project/", sp, "/Data/metadata_filtered.xlsx")))
  counts <- read.table(here(paste0("Project/", sp, "/Data/compiled_counts_table.txt")), header = T)
  
  # --- Keep only samples present in the count matrix ---
    coldata <- coldata %>%
    filter(Run %in% colnames(counts))
  
    # --- Count CONTROL and STRESSED samples per BioProject ---
    control_numbers <- coldata %>%
    group_by(BioProject, Type.of.sample) %>%
    summarise(number_control = n()) %>% ungroup() %>% filter(Type.of.sample == "CONTROL")
  
  stressed_numbers <- coldata %>%
    group_by(BioProject, Type.of.sample) %>%
    summarise(number_stressed = n()) %>% ungroup() %>% filter(Type.of.sample == "STRESSED")
  
  treatment_numbers <- control_numbers %>%
    inner_join(stressed_numbers, by = "BioProject") %>% select(!c(Type.of.sample.x, Type.of.sample.y))
  
  # --- Keep only BioProjects with both conditions represented ---
  coldata <- coldata %>%
    left_join(treatment_numbers) %>%
    filter(number_control >= 1, number_stressed >= 1) %>%
    select(!c(number_control, number_stressed))
  
  # --- Save updated metadata ---
  write.xlsx(coldata, here(paste0("Project/", sp, "/Data/metadata_filtered_2.xlsx")))
  
}

