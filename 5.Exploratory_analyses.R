library(tidyverse)
library(openxlsx)
library(here)
library(DESeq2)
library(sva)
library(PCAtools)
library(patchwork)
library(cowplot)
library(magick)

species <- c("Arabidopsis thaliana", "Triticum aestivum", "Solanum lycopersicum")
for(sp in species) {
  print(sp)
   
  coldata <- read.xlsx(here(paste0("Project/", sp, "/Data/metadata_filtered_2.xlsx")))
  counts_gene_level <- read.table(here(paste0("Project/", sp, "/Data/counts_gene_level.txt")))
  
  rownames(coldata) <- coldata$Run
  
  dds <- DESeqDataSetFromMatrix(countData = round(counts_gene_level), colData = coldata, design = ~ Type.of.sample)
  
  metadataS <- coldata %>% filter(Type.of.sample == "STRESSED")
  metadataC <- coldata %>% filter(Type.of.sample == "CONTROL")
  sizeC <- length(metadataC$Run) / length(unique(metadataC$BioProject))
  sizeS <- length(metadataS$Run) / length(unique(metadataS$BioProject))
  
  # filtraggio dei dati, tenendo solo le righe che hanno almeno un count di 10 in almeno un certo numero di colonne
  keep <- rowSums(counts(dds) >= 10) >= min(sizeC, sizeS)
  dds <- dds[keep,]
  counts <- counts(dds)
  
  dds <- estimateSizeFactors(dds)
  counts <- counts(dds, normalized = TRUE)
  
  counts_vst <- vst(round(counts), blind = TRUE)
  
  
  #########################################
  
  pca <- PCAtools::pca(counts_vst, metadata = coldata)
  
  plot1 <- biplot(pca, colby = "Type.of.sample", shape = "Major.Category", max.overlaps = 0, colkey = c(CONTROL = 'lightblue', STRESSED = 'lightcoral'), colLegendTitle = "Treatment", legendPosition = "right", pointSize = 3, axisLabSize = 8, subtitle = sp) & labs(shape = "Organ")
  pdf(here(paste0("Project/", sp, "/Products/PCA/norm_counts_vst.pdf")), width = 8)
  print(plot1)
  dev.off()
  
  plot2 <- biplot(pca, colby = "Type.of.sample", shape = "BioProject", max.overlaps = 0, colkey = c(CONTROL = 'lightblue', STRESSED = 'lightcoral'), colLegendTitle = "Treatment", legendPosition = "right", pointSize = 3, axisLabSize = 8) + scale_shape_manual(values = 0:24)
  pdf(here(paste0("Project/", sp, "/Products/PCA/norm_counts_vst_bioproject.pdf")), width = 8)
  print(plot2)
  dev.off()
  
  plot3 <- eigencorplot(pcaobj = pca, components = getComponents(pca, 1:6), metavars = c("BioProject", "Major.Category", "Type.of.sample", "Method", "Days.of.stress", "Sampling.time", "Genotype"), cexCorval = 0.6, cexLabY = 1, cexLabX = 1, fontLabX = 1, fontLabY = 1, fontCorval = 1, plotRsquared = T, col = c("white", "red1", "red2", "red3", "red4"), scale = FALSE)
  plot3$y.limits <- c("BioProject", "Treatment", "Organ", "Method", "Days of stress", "Sampling time", "Genotype")
  pdf(here(paste0("Project/", sp, "/Products/PCA/corr_matrix_norm_counts_vst.pdf")), width = 4.9, height = 3.5)
  print(plot3)
  dev.off()
  
  plot3 <- grid::grid.grabExpr(print(plot3))
  if (sp == "Arabidopsis thaliana") {
    plot3 <- ggdraw() + draw_grob(plot3, x = 0, y = 0.06, width = 1, height = 0.8)
    plot1_2 <- plot1 + plot2 + plot_annotation(tag_levels = list(c('A'))) + plot_layout(ncol = 2, guides = 'collect', axes = 'collect') & theme(plot.tag = element_text(size = 16, face = "bold"), legend.key.size = unit(0.3, "cm"), legend.key.spacing = unit(0.2, "cm"), legend.spacing.y = unit(0.25, "cm"), legend.margin = margin(t = 0, b = 0, unit = "cm"), legend.text = element_text(size = 9), legend.title = element_text(size = 10)) & scale_colour_discrete(name = "Treatment", labels = c("CONTROL" = "control", "STRESSED" = "stressed"))
  } else if (sp == "Triticum aestivum") {
    plot3 <- ggdraw() + draw_grob(plot3, x = 0, y = 0.06, width = 1, height = 0.8)
    plot1_2 <- plot1 + plot2 + plot_annotation(tag_levels = list(c('B'))) + plot_layout(ncol = 2, guides = 'collect', axes = 'collect') & theme(plot.tag = element_text(size = 16, face = "bold"), legend.key.size = unit(0.3, "cm"), legend.key.spacing = unit(0.1, "cm"), legend.spacing.y = unit(0.2, "cm"), legend.key.spacing.y = unit(-0.05, "cm"), legend.margin = margin(t = 0, b = 0, unit = "cm"), legend.text = element_text(size = 6), legend.title = element_text(size = 8), legend.box.margin = margin(0, 25, 0, 0)) & guides(shape = guide_legend(ncol = 1)) & scale_colour_discrete(name = "Treatment", labels = c("CONTROL" = "control", "STRESSED" = "stressed"))
  } else {
    plot3 <- ggdraw() + draw_grob(plot3, x = 0, y = 0.06, width = 1, height = 0.8)
    plot1_2 <- plot1 + plot2 + plot_annotation(tag_levels = list(c('C'))) + plot_layout(ncol = 2, guides = 'collect', axes = 'collect') & theme(plot.tag = element_text(size = 16, face = "bold"), legend.key.size = unit(0.3, "cm"), legend.key.spacing = unit(0.2, "cm"), legend.spacing.y = unit(0.25, "cm"), legend.margin = margin(t = 0, b = 0, unit = "cm"), legend.text = element_text(size = 9), legend.title = element_text(size = 10)) & scale_colour_discrete(name = "Treatment", labels = c("CONTROL" = "control", "STRESSED" = "stressed"))
  }


  figure <- plot_grid(plot1_2, plot3, ncol = 2, rel_widths = c(1.4, 0.8))


  pdf(here(paste0("Project/", sp, "/Products/PCA/Figura_2.pdf")), width = 15, height = 4.5)
  print(figure)
  dev.off()

  #######################################
  tissues <- c("leaf", "root")
  for (tissue in tissues) {
    print(paste0(sp, " - ", tissue))
    metadata <- coldata %>%
      filter(Major.Category == tissue)


    counts <- counts_gene_level %>%
      select(metadata$Run)

    dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = metadata, design = ~ Type.of.sample)

    metadataS <- metadata %>% filter(Type.of.sample == "STRESSED")
    metadataC <- metadata %>% filter(Type.of.sample == "CONTROL")
    sizeC <- length(metadataC$Run) / length(unique(metadataC$BioProject))
    sizeS <- length(metadataS$Run) / length(unique(metadataS$BioProject))

    # data filtering
    keep <- rowSums(counts(dds) >= 10) >= min(sizeC, sizeS)
    dds <- dds[keep,]
    counts <- counts(dds)

    dds <- estimateSizeFactors(dds)
    counts <- counts(dds, normalized = TRUE)
    counts_vst <- vst(round(counts), blind = TRUE)


    norm_adj <- read.table(here(paste0("Project/", sp, "/Data/counts_gene_level_", tissue, "_normalised_adjusted_group_parameter.txt")))

    norm_adj_vst <- vst(round(as.matrix(norm_adj)), blind = TRUE)

    dfs <- list(norm_counts_vst = counts_vst, norm_adj_vst = norm_adj_vst)
   

    pos <- 1
    for(object in dfs) {
      print(names(dfs)[pos])
      data <- metadata %>%
        filter(Run %in% colnames(object))
      object <- as.data.frame(object) %>%
        select(data$Run)
      pca <- PCAtools::pca(as.matrix(object), metadata = data)
      
      plot1 <- PCAtools::biplot(pca, colby = "Type.of.sample", max.overlaps = 0, colLegendTitle = "Treatment", colkey = c(CONTROL = 'lightblue', STRESSED = 'lightcoral'), legendPosition = "right", pointSize = 5)
      pdf(here(paste0("Project/", sp, "/Products/PCA/", tissue, "/", names(dfs)[pos], ".pdf")), width = 8)
      print(plot1)
      dev.off()

      plot2 <- eigencorplot(pcaobj = pca, components = getComponents(pca, 1:6), metavars = c("BioProject", "Type.of.sample", "Method", "Days.of.stress", "Sampling.time", "Genotype"), cexCorval = 0.6, cexLabY = 1, cexLabX = 1, fontLabX = 1, fontLabY = 1, fontCorval = 1, plotRsquared = T, col = c("white", "red1", "red2", "red3", "red4"), scale = FALSE)
      pdf(here(paste0("Project/", sp, "/Products/PCA/", tissue, "/corr_matrix_", names(dfs)[pos], ".pdf")), width = 4.9, height = 3.5)
      print(plot2)
      dev.off()
      
      pos <- pos + 1
    }
    
    
  }
}


