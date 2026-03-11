library(tidyverse)
library(openxlsx)
library(here)
library(DESeq2)
library(sva)

species <- c("Solanum lycopersicum", "Arabidopsis thaliana", "Triticum aestivum")
for(sp in species) {
  coldata <- read.xlsx(here(paste0("Project/", sp, "/Data/metadata_filtered_2.xlsx")))
  counts_gene_level <- read.table(here(paste0("Project/", sp, "/Data/counts_gene_level.txt")))
  
  dds <- DESeqDataSetFromMatrix(countData = round(counts_gene_level), colData = coldata, design = ~ Type.of.sample)

  metadataS <- coldata %>% filter(Type.of.sample == "STRESSED")
  metadataC <- coldata %>% filter(Type.of.sample == "CONTROL")
  sizeC <- length(metadataC$Run) / length(unique(metadataC$BioProject))
  sizeS <- length(metadataS$Run) / length(unique(metadataS$BioProject))

  # filtraggio dei dati, tenendo solo le righe che hanno almeno un count di 10 in almeno un certo numero di colonne
  keep <- rowSums(counts(dds) >= 10) >= min(sizeC, sizeS)
  dds <- dds[keep,]
  dds <- estimateSizeFactors(dds)
  counts <- counts(dds, normalized = TRUE)

  if (sp == "Arabidopsis thaliana") {
    max_counts <- apply(counts, 1, max)
    threshold <- quantile(max_counts, 0.995)
    counts_fil <- counts[max_counts < threshold, ]

    batch <- coldata$BioProject
    group <- coldata$Type.of.sample
    adjusted <- ComBat_seq(counts = counts_fil, batch = batch, group = group)
  } else if (sp == "Solanum lycopersicum") {
    max_counts <- apply(counts, 1, max)
    threshold <- quantile(max_counts, 0.80)
    counts_fil <- counts[max_counts < threshold, ]

    batch <- coldata$BioProject
    group <- coldata$Type.of.sample
    adjusted <- ComBat_seq(counts = counts_fil, batch = batch, group = group)
  } else {
    batch <- coldata$BioProject
    group <- coldata$Type.of.sample
    adjusted <- ComBat_seq(counts = counts, batch = batch, group = group)
  }

  adjusted <- as.data.frame(adjusted) %>% filter(across(everything(), ~ . < 10^9))
  write.table(adjusted, here(paste0("Project/", sp, "/Data/counts_gene_level_normalised_adjusted_group_parameter.txt")), quote = FALSE, sep = "\t")

  # single tissues

  tissues <- c("leaf", "root")
  for (tissue in tissues) {
    print(paste0(sp, " - ", tissue))

    metadata <- coldata %>%
      filter(Major.Category == tissue)

    if (length(unique(metadata$BioProject)) == 1) {
      next
    }

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
      counts <- counts(dds, normalized = TRUE)

      if (sp == "Arabidopsis thaliana") {
        max_counts <- apply(counts, 1, max)
        threshold <- quantile(max_counts, 0.995)
        counts_fil <- counts[max_counts < threshold, ]

        batch <- metadata$BioProject
        group <- metadata$Type.of.sample
        adjusted <- as.data.frame(ComBat_seq(counts = counts_fil, batch = batch, group = group)) %>% select(-c("SRR19846610", "SRR19846699", "SRR19846630", "SRR19846652", "SRR19846653", "SRR19846659", "SRR19846602", "SRR19846644"))
      } else {
        batch <- metadata$BioProject
        group <- metadata$Type.of.sample
        adjusted <- as.data.frame(ComBat_seq(counts = counts, batch = batch, group = group))
        if (sp == "Solanum lycopersicum" & tissue == "leaf") {
          adjusted <- adjusted %>% select(-c("SRR10192967", "SRR10192966", "SRR10192965", "SRR10192972", "SRR10192973", "SRR10192964", "SRR23959355", "SRR23959357", "SRR23959356", "SRR23959352", "SRR23959354", "SRR23959353", "SRR7218473", "SRR7206505", "SRR7218472", "SRR7206503", "SRR7206501", "SRR7206499", "SRR7206497", "SRR7206495"))
        }
      }

      adjusted <- adjusted %>% filter(across(everything(), ~ . < 10^9))

      write.table(adjusted, here(paste0("Project/", sp, "/Data/counts_gene_level_", tissue, "_normalised_adjusted_group_parameter.txt")), quote = FALSE, sep = "\t")


      if (sp == "Arabidopsis thaliana") {
        max_counts <- apply(counts_not_normalised, 1, max)
        threshold <- quantile(max_counts, 0.90)
        counts_not_normalised_fil <- counts_not_normalised[max_counts < threshold, ]

        batch <- metadata$BioProject
        group <- metadata$Type.of.sample
        adjusted <- as.data.frame(ComBat_seq(counts = counts_not_normalised_fil, batch = batch, group = group)) %>% select(-c("SRR19846610", "SRR19846699", "SRR19846630", "SRR19846652", "SRR19846653", "SRR19846659", "SRR19846602", "SRR19846644"))
      } else {
        batch <- metadata$BioProject
        group <- metadata$Type.of.sample
        adjusted <- as.data.frame(ComBat_seq(counts = counts_not_normalised, batch = batch, group = group))
        if (sp == "Solanum lycopersicum" & tissue == "leaf") {
          adjusted <- adjusted %>% select(-c("SRR10192967", "SRR10192966", "SRR10192965", "SRR10192972", "SRR10192973", "SRR10192964", "SRR23959355", "SRR23959357", "SRR23959356", "SRR23959352", "SRR23959354", "SRR23959353", "SRR7218473", "SRR7206505", "SRR7218472", "SRR7206503", "SRR7206501", "SRR7206499", "SRR7206497", "SRR7206495"))
        }
      }

      adjusted <- adjusted %>% filter(across(everything(), ~ . < 10^9))

      write.table(adjusted, here(paste0("Project/", sp, "/Data/counts_gene_level_", tissue, "_notnormalised_adjusted_group_parameter.txt")), quote = FALSE, sep = "\t")


    }


  }
  
  if (sp == "Arabidopsis thaliana") {
    counts_gene_level <- counts_gene_level %>% select(-c("SRR19846610", "SRR19846699", "SRR19846630", "SRR19846652", "SRR19846653", "SRR19846659", "SRR19846602", "SRR19846644"))
    write.table(counts_gene_level, here(paste0("Prova/", sp, "/Data/counts_gene_level.txt")))
    coldata <- coldata %>% filter(Run %in% colnames(counts_gene_level))
    write.xlsx(coldata, here(paste0("Project/", sp, "/Data/metadata_filtered_2.xlsx")))
  }
  else if (sp == "Solanum lycopersicum") {
    coldata <- coldata %>% filter(!Run %in% c("SRR10192967", "SRR10192966", "SRR10192965", "SRR10192972", "SRR10192973", "SRR10192964", "SRR23959355", "SRR23959357", "SRR23959356", "SRR23959352", "SRR23959354", "SRR23959353", "SRR7218473", "SRR7206505", "SRR7218472", "SRR7206503", "SRR7206501", "SRR7206499", "SRR7206497", "SRR7206495"))
    counts_gene_level <- counts_gene_level %>% select(coldata$Run)
    write.table(counts_gene_level, here(paste0("Prova/", sp, "/Data/counts_gene_level.txt")))
    write.xlsx(coldata, here(paste0("Project/", sp, "/Data/metadata_filtered_2.xlsx")))
  }
}