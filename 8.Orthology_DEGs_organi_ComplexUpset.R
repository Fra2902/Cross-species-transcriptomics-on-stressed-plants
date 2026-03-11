library(tidyverse)
library(openxlsx)
library(here)
library(DESeq2)
library(networkD3)
library(circlize)
library(shadowtext)
library(UpSetR)
library(biomaRt)
library(ComplexUpset)


species <- c("Arabidopsis thaliana", "Triticum aestivum", "Solanum lycopersicum")

orthogroups <- read.csv(here("Project/N0.tsv"), sep = "\t")

colnames(orthogroups) <- c("Orthogroup", "Arabidopsis thaliana", "Solanum lycopersicum", "Triticum aestivum")

orthogroups <- orthogroups %>%
  group_by(Orthogroup) %>%
  summarise(across(everything(), ~ paste(unique(.), collapse = ", ")))

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
  
  plot <- ComplexUpset::upset(
    orthogroups_2, 
    intersect = species,
    annotations = list(
      "Expression" = ggplot(mapping=aes(x=intersection, fill = expression)) + geom_bar(position = "fill") +
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
  
  mart <- useMart("plants_mart", dataset = "athaliana_eg_gene", host = "https://plants.ensembl.org")
  results <- getBM(
    attributes = c("ensembl_gene_id", "description"),
    filters = "ensembl_gene_id",
    values = deg_list[[1]]$Gene,
    mart = mart
  )
  
  detach("package:biomaRt", unload = TRUE)
  
  orthogroups_2 <- orthogroups_2 %>% filter(`Arabidopsis thaliana` == 1, `Solanum lycopersicum` == 1, `Triticum aestivum` == 1) %>% 
    rownames_to_column("Orthogroup") %>% 
    left_join(orthogroups, by = "Orthogroup") %>%
    select(-c(`Arabidopsis thaliana.x`, `Solanum lycopersicum.x`, `Triticum aestivum.x`)) %>%
    mutate(`Arabidopsis thaliana.y` = str_split(`Arabidopsis thaliana.y`, ", ")) %>%
    unnest(`Arabidopsis thaliana.y`) %>%
    inner_join(deg_list[[1]], by = join_by(`Arabidopsis thaliana.y` == Gene)) %>%
    select(-c(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, expression.y, dataset)) %>%
    inner_join(results, by = join_by(`Arabidopsis thaliana.y` == ensembl_gene_id)) %>%
    group_by(Orthogroup) %>%
    summarise(across(everything(), ~ paste(unique(.), collapse = ", "))) %>%
    mutate(`Triticum aestivum.y` = str_split(`Triticum aestivum.y`, ", ")) %>%
    unnest(`Triticum aestivum.y`) %>%
    mutate(`Triticum aestivum.y` = str_to_upper(`Triticum aestivum.y`), 
           `Solanum lycopersicum.y` = str_to_upper(`Solanum lycopersicum.y`)) %>%
    inner_join(deg_list[[2]], by = join_by(`Triticum aestivum.y` == Gene)) %>%
    select(-c(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, expression, dataset)) %>%
    group_by(Orthogroup) %>%
    summarise(across(everything(), ~ paste(unique(.), collapse = ", "))) %>%
    mutate(`Solanum lycopersicum.y` = str_split(`Solanum lycopersicum.y`, ", ")) %>%
    unnest(`Solanum lycopersicum.y`) %>%
    inner_join(deg_list[[3]], by = join_by(`Solanum lycopersicum.y` == Gene)) %>%
    select(-c(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, expression, dataset)) %>%
    group_by(Orthogroup) %>%
    summarise(across(everything(), ~ paste(unique(.), collapse = ", ")))
  
  colnames(orthogroups_2) <- c("Orthogroup", "Expression", "Arabidopsis thaliana", "Solanum lycopersicum", "Triticum aestivum", "Description")
    
  orthogroups_2 <- orthogroups_2 %>% select(Orthogroup, `Arabidopsis thaliana`, `Triticum aestivum`, `Solanum lycopersicum`, Description, Expression)
  write.xlsx(orthogroups_2, here(paste0("Project/Orthology analysis/conserved_degs_", tissue, ".xlsx")))
}


