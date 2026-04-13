library(tidyverse)
library(openxlsx)
library(here)
library(DESeq2)
library(sva)

# List of species to process
species <- c("Solanum lycopersicum", "Arabidopsis thaliana", "Triticum aestivum")
for(sp in species) {
  coldata <- read.xlsx(here(paste0("Project/", sp, "/Data/metadata_filtered_2.xlsx")))
  counts_gene_level <- read.table(here(paste0("Project/", sp, "/Data/counts_gene_level.txt")))
  
  # --- Process each tissue separately ---
  tissues <- c("leaf", "root")
  for (tissue in tissues) {
    print(paste0(sp, " - ", tissue))
    
    metadata <- coldata %>%
      filter(Major.Category == tissue)
    
    if (length(unique(metadata$BioProject)) > 1) {
      metadata_numbers <- metadata %>%
        group_by(BioProject) %>%
        summarise(count = n()) %>%
        ungroup()
      
      metadata <- metadata %>%
        left_join(metadata_numbers) %>%
        filter(count > 2)
      
      counts <- counts_gene_level %>%
        select(metadata$Run)
      
      print(identical(colnames(counts), metadata$Run))
      dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = metadata, design = ~ Type.of.sample)
      
      # filtraggio dei dati, tenendo solo le righe che hanno almeno un count di 10 in almeno un certo numero di colonne
      metadataS <- metadata %>% filter(Type.of.sample == "STRESSED")
      metadataC <- metadata %>% filter(Type.of.sample == "CONTROL")
      sizeC <- length(metadataC$Run) / length(unique(metadataC$BioProject))
      sizeS <- length(metadataS$Run) / length(unique(metadataS$BioProject))
      
      keep <- rowSums(counts(dds) >= 10) >= min(sizeC, sizeS)
      dds <- dds[keep,]
      counts_not_normalised <- assay(dds)
      dds <- estimateSizeFactors(dds)
      counts_fil <- counts(dds, normalized = TRUE)
      
      if (sp == "Arabidopsis thaliana" & tissue == "leaf") {
        max_counts <- apply(counts_fil, 1, max)
        threshold <- quantile(max_counts, 0.995)
        counts_fil <- counts_fil[max_counts < threshold, ]
        
        batch <- metadata$BioProject
        group <- metadata$Type.of.sample
        adjusted <- as.data.frame(ComBat_seq(counts = counts_fil, batch = batch, group = group)) %>% select(-c("SRR19846610", "SRR19846699", "SRR19846630", "SRR19846652", "SRR19846653", "SRR19846659", "SRR19846602", "SRR19846644"))
      } else {
        batch <- metadata$BioProject
        group <- metadata$Type.of.sample
        adjusted <- as.data.frame(ComBat_seq(counts = counts_fil, batch = batch, group = group))
        if (sp == "Solanum lycopersicum" & tissue == "leaf") {
          adjusted <- adjusted %>% select(-c("SRR10192967", "SRR10192966", "SRR10192965", "SRR10192972", "SRR10192973", "SRR10192964", "SRR23959355", "SRR23959357", "SRR23959356", "SRR23959352", "SRR23959354", "SRR23959353", "SRR7218473", "SRR7206505", "SRR7218472", "SRR7206503", "SRR7206501", "SRR7206499", "SRR7206497", "SRR7206495"))
        }
      }
      
      adjusted <- adjusted %>% filter(if_all(everything(), ~ . < 10^9))
      
      write.table(adjusted, here(paste0("Project/", sp, "/Data/counts_gene_level_", tissue, "_normalised_adjusted.txt")), quote = FALSE, sep = "\t")
      
      
    }
    
    
  }
  
  if (sp == "Arabidopsis thaliana") {
    counts_gene_level <- counts_gene_level %>% select(-c("SRR19846610", "SRR19846699", "SRR19846630", "SRR19846652", "SRR19846653", "SRR19846659", "SRR19846602", "SRR19846644"))
    
  } else if (sp == "Solanum lycopersicum") {
    counts_gene_level <- counts_gene_level %>% select(-c("SRR10192967", "SRR10192966", "SRR10192965", "SRR10192972", "SRR10192973", "SRR10192964", "SRR23959355", "SRR23959357", "SRR23959356", "SRR23959352", "SRR23959354", "SRR23959353", "SRR7218473", "SRR7206505", "SRR7218472", "SRR7206503", "SRR7206501", "SRR7206499", "SRR7206497", "SRR7206495"))
  }
  write.table(counts_gene_level, here(paste0("Project/", sp, "/Data/counts_gene_level_2.txt")))
  coldata <- coldata %>% filter(Run %in% colnames(counts_gene_level))
  write.xlsx(coldata, here(paste0("Project/", sp, "/Data/metadata_filtered_3.xlsx")))
}