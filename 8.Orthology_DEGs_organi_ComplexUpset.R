library(tidyverse)
library(openxlsx)
library(here)
library(DESeq2)
library(networkD3)
library(circlize)
library(shadowtext)
library(ComplexUpset)


species <- c("Arabidopsis thaliana", "Triticum aestivum", "Solanum lycopersicum")

# Import the Orthofinder output
orthogroups <- read.csv(here("Project/N0.tsv"), sep = "\t")

colnames(orthogroups) <- c("Orthogroup", "Arabidopsis thaliana", "Solanum lycopersicum", "Triticum aestivum")

orthogroups <- orthogroups %>%
  group_by(Orthogroup) %>%
  summarise(across(everything(), ~ paste(unique(.), collapse = ", ")))

mart <- biomaRt::useMart("plants_mart", dataset = "athaliana_eg_gene", host = "https://plants.ensembl.org")

for (tissue in c("leaf", "root")) {
  
  deg_list <- lapply(species, function(sp) {
    read.xlsx(here(paste0("Project/", sp, "/Products/DEGs analysis/GDE_", tissue, ".xlsx")))
  })
  
  names(deg_list) <- species
  
  orthogroups_2 <- orthogroups %>% column_to_rownames("Orthogroup")
  
  for (i in 1:length(orthogroups$`Arabidopsis thaliana`)) {
    
    up <- 0
    down <- 0
    
    for (sp in species) {
      
      
      ort_genes <- str_split_1(str_to_upper(orthogroups[[sp]][i]), ", ")
      if ("TRUE" %in% (ort_genes %in% deg_list[[sp]]$Gene)) {
        
        orthogroups_2[[sp]][i] <- 1
        
        for (gene in ort_genes) {
          if (gene %in% deg_list[[sp]]$Gene) {
            if ((deg_list[[sp]] %>% filter(Gene == gene))$expression == "OVER") {
              up <- up + 1
              
            } else {
              down <- down + 1
              
            }
          }
        }
      } else {
        orthogroups_2[[sp]][i] <- 0
      }
      
      
    }
    
    if (up > down) {
      orthogroups_2$expression[i] <- "UP"
    } else if (up < down) {
      orthogroups_2$expression[i] <- "DOWN"
    } else {
      orthogroups_2$expression[i] <- "EVEN"
    }
    
    if (i == length(orthogroups_2$`Arabidopsis thaliana`)) {
      break
    }
  }
  orthogroups_2 <- orthogroups_2 %>% filter(!(`Arabidopsis thaliana` == 0 & `Solanum lycopersicum` == 0 & `Triticum aestivum` == 0))
  orthogroups_2$`Arabidopsis thaliana` <- as.numeric(orthogroups_2$`Arabidopsis thaliana`)
  orthogroups_2$`Solanum lycopersicum` <- as.numeric(orthogroups_2$`Solanum lycopersicum`)
  orthogroups_2$`Triticum aestivum` <- as.numeric(orthogroups_2$`Triticum aestivum`)
  
  plot <- upset(
    orthogroups_2, 
    intersect = species,
    annotations = list(
      "Expression" = ggplot(mapping=aes(x=intersection, fill = expression)) + 
        geom_bar(position = "fill") +
        scale_fill_manual(values = c(
          "UP" = "salmon",
          "DOWN" = "lightblue",
          "EVEN" = "orange"
        )) 
    )
  )
  
  pdf(here(paste0("Project/Orthology analysis/orthology_ComplexUpset_", tissue, ".pdf")))
  print(plot)
  dev.off()
  
  
  results <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "description"),
    filters = "ensembl_gene_id",
    values = deg_list[[1]]$Gene,
    mart = mart
  )
  
  
  orthogroups_2 <- orthogroups_2 %>% filter(`Arabidopsis thaliana` == 1, `Solanum lycopersicum` == 1, `Triticum aestivum` == 1) %>% 
    rownames_to_column("Orthogroup") %>% 
    left_join(orthogroups, by = "Orthogroup") %>%
    select(-c(`Arabidopsis thaliana.x`, `Solanum lycopersicum.x`, `Triticum aestivum.x`)) %>%
    mutate(`Arabidopsis thaliana.y` = str_split(`Arabidopsis thaliana.y`, ", ")) %>%
    unnest(`Arabidopsis thaliana.y`) %>%
    inner_join(deg_list[["Arabidopsis thaliana"]], by = join_by(`Arabidopsis thaliana.y` == Gene)) %>%
    select(-c(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, expression.y, dataset)) %>%
    inner_join(results, by = join_by(`Arabidopsis thaliana.y` == ensembl_gene_id)) %>%
    group_by(Orthogroup) %>%
    summarise(across(everything(), ~ paste(unique(.), collapse = ", "))) %>%
    mutate(`Triticum aestivum.y` = str_split(`Triticum aestivum.y`, ", ")) %>%
    unnest(`Triticum aestivum.y`) %>%
    mutate(`Triticum aestivum.y` = str_to_upper(`Triticum aestivum.y`), 
           `Solanum lycopersicum.y` = str_to_upper(`Solanum lycopersicum.y`)) %>%
    inner_join(deg_list[["Triticum aestivum"]], by = join_by(`Triticum aestivum.y` == Gene)) %>%
    select(-c(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, expression, dataset)) %>%
    group_by(Orthogroup) %>%
    summarise(across(everything(), ~ paste(unique(.), collapse = ", "))) %>%
    mutate(`Solanum lycopersicum.y` = str_split(`Solanum lycopersicum.y`, ", ")) %>%
    unnest(`Solanum lycopersicum.y`) %>%
    inner_join(deg_list[["Solanum lycopersicum"]], by = join_by(`Solanum lycopersicum.y` == Gene)) %>%
    select(-c(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, expression, dataset)) %>%
    group_by(Orthogroup) %>%
    summarise(across(everything(), ~ paste(unique(.), collapse = ", ")))
  
  colnames(orthogroups_2) <- c("Orthogroup", "Expression", "Arabidopsis thaliana", "Solanum lycopersicum", "Triticum aestivum", "Description")
    
  orthogroups_2 <- orthogroups_2 %>% select(Orthogroup, `Arabidopsis thaliana`, `Triticum aestivum`, `Solanum lycopersicum`, Description, Expression)
  write.xlsx(orthogroups_2, here(paste0("Project/Orthology analysis/conserved_degs_", tissue, ".xlsx")))

  degs <- orthogroups_2 %>%
    select(-Expression) %>%
    mutate(`Arabidopsis thaliana` = str_split(`Arabidopsis thaliana`, ", ")) %>% 
    unnest(`Arabidopsis thaliana`) %>% 
    inner_join(read.xlsx(here(paste0("Project/Arabidopsis thaliana/Products/DEGs analysis/GDE_", tissue, ".xlsx"))) %>% 
                select(Gene, log2FoldChange), 
              by = join_by(`Arabidopsis thaliana` == Gene)) %>%
    group_by(Orthogroup) %>% 
    summarise(log2FoldChange = mean(log2FoldChange),
              `Arabidopsis thaliana` = paste(`Arabidopsis thaliana`, collapse = ", "),
              across(-c(log2FoldChange, `Arabidopsis thaliana`), dplyr::first),
              .groups = "drop") %>%
    mutate(`Triticum aestivum` = str_split(`Triticum aestivum`, ", ")) %>% 
    unnest(`Triticum aestivum`) %>% 
    inner_join(read.xlsx(here(paste0("Project/Triticum aestivum/Products/DEGs analysis/GDE_", tissue, ".xlsx"))) %>% 
                select(Gene, log2FoldChange),
              by = join_by(`Triticum aestivum` == Gene),
              suffix = c(".ath", ".tae")) %>% 
    group_by(Orthogroup) %>% 
    summarise(log2FoldChange.tae = mean(log2FoldChange.tae),
              `Triticum aestivum` = paste(`Triticum aestivum`, collapse = ", "),
              across(-c(log2FoldChange.tae, `Triticum aestivum`), dplyr::first),
              .groups = "drop") %>%
    mutate(`Solanum lycopersicum` = str_split(`Solanum lycopersicum`, ", ")) %>%
    unnest(`Solanum lycopersicum`) %>%
    inner_join(read.xlsx(here(paste0("Project/Solanum lycopersicum/Products/DEGs analysis/GDE_", tissue, ".xlsx"))) %>%
                select(Gene, log2FoldChange),
              by = join_by(`Solanum lycopersicum` == Gene)) %>% 
    group_by(Orthogroup) %>% 
    summarise(log2FoldChange.sly = mean(log2FoldChange),
              `Solanum lycopersicum` = paste(`Solanum lycopersicum`, collapse = ", "),
              across(-c(log2FoldChange.sly, `Solanum lycopersicum`), dplyr::first),
              .groups = "drop") %>% 
    mutate(meanLF = rowMeans(across(c(log2FoldChange.sly, log2FoldChange.tae, log2FoldChange.ath))))
  
  write.xlsx(degs, here(paste0("Project/Orthology analysis/conserved_degs_", tissue, "_lfc.xlsx")))
  
}


