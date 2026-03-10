library(tidyverse)
library(openxlsx)
library(here)
library(DESeq2)
library(gplots)

species <- c("Arabidopsis thaliana", "Triticum aestivum", "Solanum lycopersicum")
for(sp in species) {
  
  coldata <- read.xlsx(here(paste0("Prova/", sp, "/Data/metadata_filtered_2.xlsx")))
  
  for (tissue in c("leaf", "root")) {
    counts <- read.table(here(paste0("Prova/", sp, "/Data/counts_gene_level_", tissue, "_normalised_adjusted_group_parameter.txt")))
    
    metadata <- coldata %>%
      filter(Major.Category == tissue, Run %in% colnames(counts))
    counts <- counts %>%
      select(metadata$Run)
    
    if (length(unique(metadata$Type.of.sample)) == 2 & length(metadata$Run) > 2) {
      
      dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = metadata, design = ~ Type.of.sample)
      metadataS <- metadata %>% filter(Type.of.sample == "STRESSED")
      metadataC <- metadata %>% filter(Type.of.sample == "CONTROL")
      sizeC <- length(metadataC$Run) / length(unique(metadataC$BioProject))
      sizeS <- length(metadataS$Run) / length(unique(metadataS$BioProject))
      keep <- rowSums(counts(dds) >= 10) >= min(sizeC, sizeS)
      dds <- dds[keep,]
      write(rownames(assay(dds)), here(paste0("Prova/", sp, "/Products/DEGs analysis/expressed_genes_", tissue, ".txt")))
      dds <- DESeq(dds)
      res <- results(dds, contrast = c("Type.of.sample", "STRESSED", "CONTROL"), alpha = 0.05)
      res <- lfcShrink(dds, contrast= c("Type.of.sample", "STRESSED", "CONTROL"), type="normal", res=res)
      res <- as.data.frame(res) %>%
        filter(abs(log2FoldChange) >= log2(1.5) & padj <= 0.05) %>%
        rownames_to_column(var = "Gene") %>%
        mutate(expression = case_when(log2FoldChange < 0 ~ 'UNDER', log2FoldChange > 0 ~ 'OVER'), dataset = tissue)
      write.xlsx(res, here(paste0("Prova/", sp, "/Products/DEGs analysis/GDE_", tissue, ".xlsx", sep = "")), quote = FALSE)
      
      
      studies <- unique(metadata$BioProject)
      
      if (length(studies) > 1) {
        for (study in studies) {
          data <- metadata %>%
            filter(BioProject == study)
          micro_counts <- counts %>%
            select(data$Run)
          
          
          if (length(unique(data$Type.of.sample)) == 2 & length(data$Run) > 2) {
            dds <- DESeqDataSetFromMatrix(countData = round(micro_counts), colData = data, design = ~ Type.of.sample)
            metadataS <- data %>% filter(Type.of.sample == "STRESSED")
            metadataC <- data %>% filter(Type.of.sample == "CONTROL")
            sizeC <- length(metadataC$Run) / length(unique(metadataC$BioProject))
            sizeS <- length(metadataS$Run) / length(unique(metadataS$BioProject))
            keep <- rowSums(counts(dds) >= 10) >= min(sizeC, sizeS)
            dds <- dds[keep,]
            write(rownames(assay(dds)), here(paste0("Prova/", sp, "/Products/DEGs analysis/expressed_genes_", study, "_", tissue, ".txt")))
            dds <- DESeq(dds)
            res <- results(dds, contrast = c("Type.of.sample", "STRESSED", "CONTROL"), alpha = 0.05)
            res <- lfcShrink(dds, contrast= c("Type.of.sample", "STRESSED", "CONTROL"), type="normal", res=res)
            res <- as.data.frame(res) %>%
              filter(abs(log2FoldChange) >= log2(1.5) & padj <= 0.05) %>%
              rownames_to_column(var = "Gene") %>%
              mutate(expression = case_when(log2FoldChange < 0 ~ 'UNDER', log2FoldChange > 0 ~ 'OVER'))
            write.xlsx(res, here(paste0("Prova/", sp, "/Products/DEGs analysis/GDE_", study, "_", tissue, ".xlsx", sep = "")), quote = FALSE)
            
            res <- res %>%
              mutate(dataset = paste0(study, "_", tissue, "_"))
          }  
          
        }
      }
    }
    
  }
  
  ## barplot
  projects <- list.files(path = here(paste0("Prova/", sp, "/Products/DEGs analysis/")), pattern = "\\.xlsx$", full.names = TRUE)
  
  exclude_pattern <- "jaccard"
  
  # Exclude the specified file(s) from the list
  projects <- projects[!grepl(exclude_pattern, projects)]
  
  df <- projects %>%
    keep(~ file.info(.x)$size > 6450) %>%
    map_dfr(~ read.xlsx(.x) %>% mutate(Gene = as.character(Gene), expression = as.character(expression), study = substring(basename(.x), 5, nchar(basename(.x))-5)))
  
  barplot <- ggplot() +
    geom_bar(data=df, mapping=aes(x=study, fill=expression), color = "black", position = position_dodge()) +
    scale_fill_manual(values = c("OVER" = "salmon", "UNDER" = "deepskyblue")) +
    scale_y_continuous(trans = "log10") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  pdf(here(paste0("Prova/", sp, "/Products/DEGs analysis/barplot.pdf")), width = 4.9, height = 4.9)
  print(barplot)
  dev.off()
  
  
  
  ## jaccard index heatmap
  
  projects <- projects %>%
    map(~ substring(basename(.x), 5, nchar(basename(.x))-5))
  projects <- unlist(projects)
  
  for (i in c("OVER", "UNDER")) {
    jaccard_index_df <- data.frame(matrix(nrow=length(projects),ncol=length(projects)))
    rownames(jaccard_index_df) <- projects
    colnames(jaccard_index_df) <- projects
    for (row in 1:length(projects)) {
      for (col in 1:length(projects)) {
        if (col == row) {
          jaccard_index_df[row, col] = NA
        }
        else {
          row_name <- rownames(jaccard_index_df)[row]
          col_name <- colnames(jaccard_index_df)[col]
          
          
          df1 <- read.xlsx(here(paste0("Prova/", sp, "/Products/DEGs analysis/GDE_", row_name, ".xlsx", sep = "")), rowNames = TRUE) %>%
            filter(expression == i)
          
          df2 <- read.xlsx(here(paste0("Prova/", sp, "/Products/DEGs analysis/GDE_", col_name, ".xlsx", sep = "")), rowNames = TRUE) %>%
            filter(expression == i)
          
          genes1 <- rownames(df1)
          genes2 <- rownames(df2)
          jaccard_index <- length(intersect(genes1, genes2)) / length(union(genes1, genes2))
          jaccard_index_df[row, col] <- jaccard_index
        }
      }
      rownames(jaccard_index_df)[row] <- paste0(projects[row], " - ", length(rownames(df1)), " genes")
    }
    
    jaccard_index_df[upper.tri(jaccard_index_df)] <- NA
    my_palette <- colorpanel(10, "white", "blue")
    pdf(here(paste0("Prova/", sp, "/Products/DEGs analysis/heatmap_jaccard_index_", i, ".pdf")), width = 8.6, height = 8.6)
    heatmap.2(as.matrix(jaccard_index_df), trace="none", density.info="none", dendrogram="none", sepcolor = "black", colsep = 1:ncol(jaccard_index_df), rowsep = 1:ncol(jaccard_index_df), sepwidth = c(0.0005, 0.0005), col=my_palette, margins=c(12,12), cexRow=1, cexCol=1, Rowv = NA, Colv = NA, offsetRow = -40)
    dev.off()
    
    jaccard_index_df <- jaccard_index_df %>%
      rownames_to_column(var = "Subset")
    write.xlsx(jaccard_index_df, here(paste0("Prova/", sp, "/Products/DEGs analysis/jaccard_index_", i, ".xlsx")), quote = FALSE)
  }
  
  
}



#########################################

# boxplot JIs

jaccard_indexes <- list()
expression <- list()
organisms <- list()
tissues <- list()
for (i in c("OVER", "UNDER")) {
  
  expr <- ifelse(i == "OVER", "UP", "DOWN")
  
  for (sp in species) {
    index_table <- read.xlsx(here(paste0("Prova/", sp, "/Products/DEGs analysis/jaccard_index_", i, ".xlsx")))
    index_table <- index_table %>% column_to_rownames("Subset")
    
    for (tissue in c("leaf", "root")) {
      
      tissue_table <- index_table %>% filter(str_starts(rownames(index_table), "P"), str_detect(rownames(index_table), tissue)) %>% select(!c("leaf", "root") & ends_with(tissue))
      JIs <- tissue_table %>% pivot_longer(everything(), names_to = "Set", values_to = "JI") %>% drop_na(JI) %>% pull(JI)
    
      expression <- append(expression, rep(expr, length(JIs)))
      organisms <- append(organisms, rep(sp, length(JIs)))
      tissues <- append(tissues, rep(tissue, length(JIs)))
      jaccard_indexes <- append(jaccard_indexes, JIs)
      
    }
    
  }
  
}

expression <- unlist(expression)
organisms <- unlist(organisms)
tissues <- unlist(tissues)
jaccard_indexes <- unlist(jaccard_indexes)
df_boxplot <- data_frame(Expression = expression, Species = organisms, Organ = tissues, Index = jaccard_indexes)
df_boxplot$Expression <- factor(df_boxplot$Expression, levels = c("UP", "DOWN"))
write.xlsx(df_boxplot, here("Prova/general_jaccard_index_sets対sets.xlsx"))
boxplot <- ggplot(data=df_boxplot, mapping=aes(x = interaction(Expression, Organ, Species), y = Index, fill = Organ)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  scale_x_discrete(labels = c(
    "UP.leaf.Arabidopsis thaliana" = "Arabidopsis leaf UP",
    "UP.root.Arabidopsis thaliana" = "Arabidopsis root UP",
    "UP.leaf.Solanum lycopersicum" = "Tomato leaf UP",
    "UP.root.Solanum lycopersicum" = "Tomato root UP",
    "UP.leaf.Triticum aestivum" = "Wheat leaf UP",
    "UP.root.Triticum aestivum" = "Wheat root UP", 
    "DOWN.leaf.Arabidopsis thaliana" = "Arabidopsis leaf DOWN",
    "DOWN.root.Arabidopsis thaliana" = "Arabidopsis root DOWN",
    "DOWN.leaf.Solanum lycopersicum" = "Tomato leaf DOWN",
    "DOWN.root.Solanum lycopersicum" = "Tomato root DOWN",
    "DOWN.leaf.Triticum aestivum" = "Wheat leaf DOWN",
    "DOWN.root.Triticum aestivum" = "Wheat root DOWN"
  )) +
  scale_fill_manual(values = c("leaf" = "darkgreen", "root" = "brown")) +
  labs(x = "Set", y = "Jaccard index") +
  facet_wrap(~ Expression, scale = "free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 


pdf(here("Prova/general_jaccard_index_sets対sets.pdf"))
print(boxplot)
dev.off()

# expression type and organ DEG set similarity difference significance calculation

for (sp in species) {
  sh <- shapiro.test((df_boxplot %>% filter(Species == sp))$Index)
  pdf(here(paste0("Prova/", sp, "/Products/DEGs analysis/JI_distribution.pdf")))
  print(hist((df_boxplot %>% filter(Species == sp))$Index))
  dev.off()
  
  ## ttest requirements and test calculation
  if (sh[[2]] > 0.05) {
    ba <- bartlett.test(Index ~ Expression, data = df_boxplot)
    
    if (ba[[3]] > 0.05) {
      mod <- t.test(Index ~ Expression, data = df_boxplot, var.equal = TRUE)
      write_lines(mod, here(paste0("Prova/", sp, "/Products/DEGs analysis/ttest_expression.txt")))
    } else {
      mod <- t.test(Index ~ Expression, data = df_boxplot)
      write_lines(mod, here(paste0("Prova/", sp, "/Products/DEGs analysis/welch_expression.txt")))
    }
    
    ba <- bartlett.test(Index ~ Organ, data = df_boxplot)
    if (ba[[3]] > 0.05) {
      mod <- t.test(Index ~ Organ, data = df_boxplot, var.equal = TRUE)
      write_lines(mod, here(paste0("Prova/", sp, "/Products/DEGs analysis/ttest_organ.txt")))
    } else {
      mod <- t.test(Index ~ Expression, data = df_boxplot)
      write_lines(mod, here(paste0("Prova/", sp, "/Products/DEGs analysis/welch_expression.txt")))
    }
    
  } else {
    ## wilcoxon-mann-whitney calculation
    
    wi <- wilcox.test((df_boxplot %>% filter(Species == sp, Expression == "UP"))$Index, (df_boxplot %>% filter(Species == sp, Expression == "DOWN"))$Index)
    write_lines(wi, here(paste0("Prova/", sp, "/Products/DEGs analysis/wilcoxon_expression.txt")))
    
    wi <- wilcox.test((df_boxplot %>% filter(Species == sp, Organ == "leaf"))$Index, (df_boxplot %>% filter(Species == sp, Organ == "root"))$Index)
    write_lines(wi, here(paste0("Prova/", sp, "/Products/DEGs analysis/wilcoxon_organ.txt")))
  }
  
  
  
}



####################################
# boxplot DEG number

degs_number <- list()
expression <- list()
organisms <- list()
tissues <- list()
integrated_sets <- list()

for (sp in species) {
  projects <- list.files(path = here(paste0("Prova/", sp, "/Products/DEGs analysis/")), pattern = "\\.xlsx$", full.names = TRUE)
  exclude_pattern <- "jaccard"
  projects <- projects[!grepl(exclude_pattern, projects)]
  
  projects <- projects %>%
    map(~ substring(basename(.x), 5, nchar(basename(.x))-5))
  projects <- unlist(projects)
  
  for (name in projects) {
    df_degs <- read.xlsx(here(paste0("Prova/", sp, "/Products/DEGs analysis/GDE_", name, ".xlsx", sep = "")), rowNames = TRUE)
    for (i in c("OVER", "UNDER")) {
      degs <- df_degs %>% filter(expression == i)
      
      expr <- ifelse(i == "OVER", "UP", "DOWN")
      
      if (!name %in% c("leaf", "root")) {
        tissue <- str_sub(name, -4, -1)
        
        expression <- append(expression, expr)
        organisms <- append(organisms, sp)
        tissues <- append(tissues, tissue)
        degs_number <- append(degs_number, length(rownames(degs)))
      } else {
        integrated_sets <- append(integrated_sets, length(rownames(degs)))
      }
      
    }
      
  }
}
expression <- unlist(expression)
organisms <- unlist(organisms)
tissues <- unlist(tissues)
degs_number <- unlist(degs_number)
df_boxplot <- data_frame(Expression = expression, Species = organisms, Organ = tissues, DEGs_number = degs_number)
df_boxplot$Expression <- factor(df_boxplot$Expression, levels = c("UP", "DOWN"))
df_boxplot <- df_boxplot %>%
  mutate(Integrated = case_when(
    Species == "Arabidopsis thaliana" & Organ == "leaf" & Expression == "UP" ~ integrated_sets[[1]],
    Species == "Arabidopsis thaliana" & Organ == "leaf" & Expression == "DOWN" ~ integrated_sets[[2]],
    Species == "Arabidopsis thaliana" & Organ == "root" & Expression == "UP" ~ integrated_sets[[3]],
    Species == "Arabidopsis thaliana" & Organ == "root" & Expression == "DOWN" ~ integrated_sets[[4]],
    Species == "Triticum aestivum" & Organ == "leaf" & Expression == "UP" ~ integrated_sets[[5]],
    Species == "Triticum aestivum" & Organ == "leaf" & Expression == "DOWN" ~ integrated_sets[[6]],
    Species == "Triticum aestivum" & Organ == "root" & Expression == "UP" ~ integrated_sets[[7]],
    Species == "Triticum aestivum" & Organ == "root" & Expression == "DOWN" ~ integrated_sets[[8]],
    Species == "Solanum lycopersicum" & Organ == "leaf" & Expression == "UP" ~ integrated_sets[[9]],
    Species == "Solanum lycopersicum" & Organ == "leaf" & Expression == "DOWN" ~ integrated_sets[[10]],
    Species == "Solanum lycopersicum" & Organ == "root" & Expression == "UP" ~ integrated_sets[[11]],
    Species == "Solanum lycopersicum" & Organ == "root" & Expression == "DOWN" ~ integrated_sets[[12]]
  ))

boxplot <- ggplot(data=df_boxplot, mapping=aes(x = interaction(Expression, Organ, Species), y = DEGs_number, fill = Organ)) +
  geom_boxplot(width = 0.7) +
  geom_crossbar(aes(x = interaction(Expression, Organ, Species), y = Integrated, ymin = Integrated, ymax = Integrated), colour = "salmon", width = 0.7, linewidth = 0.2) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  scale_x_discrete(labels = c(
    "UP.leaf.Arabidopsis thaliana" = "Arabidopsis leaf UP",
    "UP.root.Arabidopsis thaliana" = "Arabidopsis root UP",
    "UP.leaf.Solanum lycopersicum" = "Tomato leaf UP",
    "UP.root.Solanum lycopersicum" = "Tomato root UP",
    "UP.leaf.Triticum aestivum" = "Wheat leaf UP",
    "UP.root.Triticum aestivum" = "Wheat root UP", 
    "DOWN.leaf.Arabidopsis thaliana" = "Arabidopsis leaf DOWN",
    "DOWN.root.Arabidopsis thaliana" = "Arabidopsis root DOWN",
    "DOWN.leaf.Solanum lycopersicum" = "Tomato leaf DOWN",
    "DOWN.root.Solanum lycopersicum" = "Tomato root DOWN",
    "DOWN.leaf.Triticum aestivum" = "Wheat leaf DOWN",
    "DOWN.root.Triticum aestivum" = "Wheat root DOWN"
  )) +
  scale_fill_manual(values = c("leaf" = "darkgreen", "root" = "brown")) +
  labs(x = "Set", y = "DEGs number") +
  facet_wrap(~ Expression, scale = "free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

pdf(here("Prova/general_number_of_degs.pdf"))
print(boxplot)
dev.off()

