library(tidyverse)
library(openxlsx)
library(here)
library(DESeq2)
library(sva)
library(PCAtools)
library(gplots)
library(pheatmap)
library(enrichplot)
library(DOSE)

species <- c("Arabidopsis thaliana", "Triticum aestivum", "Solanum lycopersicum")
generic_BP <- read_file(here("generic_BP.txt"))
generic_BP <- str_split_1(generic_BP, "\n")

conserved_terms <- list()
i <- 1
for(sp in species) {
  organism <- paste0(str_to_lower(substr(sp, 1, 1)), str_split_i(sp, " ", 2))
  gmt_token <- gprofiler2::upload_GMT_file(gmtfile = here(paste0("Prova/", sp, "/Data/", organism, ".GO_BP.ENSG_PlantConnectomeDroughtAdded.gmt")))
  
  coldata <- read.xlsx(here(paste0("Prova/", sp, "/Data/metadata_filtered_2.xlsx")))
  counts_gene_level <- read.table(here(paste0("Prova/", sp, "/Data/counts_gene_level.txt")))
  
  terms <- list()
  tissues <- c("leaf", "root")
  
  for (tissue in tissues) {
    
    expressed_genes_adj <- read_file(here(paste0("Prova/", sp, "/Products/DEGs analysis/expressed_genes_", tissue, ".txt")))
    expressed_genes_adj <- str_split_1(expressed_genes_adj, "\n")
    
    degs <- read.xlsx(here(paste0("Prova/", sp, "/Products/DEGs analysis/GDE_", tissue, ".xlsx")))
    
    for (expr in c("OVER", "UNDER")) {
      genes <- (degs %>% filter(expression == expr))$Gene
      
      gost <- gprofiler2::gost(query = genes, organism = organism, highlight = F, custom_bg = expressed_genes_adj, evcodes = T)
      df <- data.frame(gost$result$term_id, gost$result$source, gost$result$term_name, gost$result$p_value, gost$result$query, gost$result$term_size, gost$result$query_size, gost$result$intersection_size, gost$result$intersection)
      if (length(df) > 0) {
        colnames(df) <- c("ID", "Source", "Description", "p.adjust", "Query", "Term_size", "Query_size", "Intersection_size", "geneID")
        df <- df %>%
          filter(!(ID %in% generic_BP))
        
        write.xlsx(df, here(paste0("Prova/", sp, "/Products/Enrichment analysis/enrichment_analysis_degs_", tissue, "_", expr, ".xlsx")))
        terms <- append(terms, df$Description)
        
        if (length(df$ID) > 0) {
          df$Count <- sapply(df$geneID, function(x) length(strsplit(x, ",")[[1]]))
          df$GeneRatio <- paste0(df$Count, "/", df$Term_size)
          df$pvalue <- df$p.adjust
          df$qvalue <- df$p.adjust
          
          er <- new("enrichResult",
                    result = df,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "gSCS",
                    qvalueCutoff = 0.2,
                    organism = substr(organism, 1, 3),  
                    ontology = "BP",
                    gene = unique(unlist(strsplit(df$geneID, ","))),
                    universe = rownames(assay(dds)),
                    geneSets = list(),  # Not strictly required unless doing downstream analysis
                    readable = FALSE)
          
          dotplot <- dotplot(er,
                             showCategory = 25,
                             color = "p.adjust",
                             size = "Count",
                             font.size = 12,
                             orderBy = "p.adjust",
                             label_format = 50) + 
            scale_fill_gradientn(
              colors = c("red", "blue"),
              trans = "log10",  # log scale to spread values more evenly
              name = "p.adjust"
            ) 
          
          ggsave(here(paste0("Prova/", sp, "/Products/Enrichment analysis/enrichment_dotplot_degs_", tissue, "_", expr, ".pdf")), dotplot, height = 6, width = 8, units = "in", dpi = 600)
        }
        
        
      }
      
    }
    
    
    
  }
  terms <- unique(terms)
  conserved_terms[[i]] <- terms
  i <- i + 1
}
conserved_terms <- union(union(intersect(conserved_terms[[1]], conserved_terms[[2]]), intersect(conserved_terms[[1]], conserved_terms[[3]])), intersect(conserved_terms[[2]], conserved_terms[[3]]))
conserved_terms <- unlist(unique(conserved_terms))
write(conserved_terms, here("Prova/conserved_enriched_terms_organism.txt"))

heatmap_df <- data.frame(row.names = conserved_terms)
for (sp in species) {
  for (tissue in tissues) {
    expression <- list()
    terms_list <- lapply(c("OVER", "UNDER"), function(expr) {
      enr <- read.xlsx(here(paste0("Prova/", sp, "/Products/Enrichment analysis/enrichment_analysis_degs_", tissue, "_", expr, ".xlsx")))
      enr <- enr %>% filter(Description %in% conserved_terms)
      
      return(enr)
    })
    
    if (length(terms_list[[1]]$Description) > 0 & length(terms_list[[2]]$Description) > 0) {
      terms_list <- terms_list[[1]] %>% full_join(terms_list[[2]], by = "Description")
      for (term in conserved_terms) {
        temp_list <- terms_list %>% filter(Description == term)
        if (length(temp_list$Description) == 0) {
          expression <- append(expression, "NO")
        } else {
          if (!is.na(temp_list$p.adjust.x) & !is.na(temp_list$p.adjust.y)) {
            expression <- append(expression, 0)
          } else if (!is.na(temp_list$p.adjust.x) & is.na(temp_list$p.adjust.y)) {
            expression <- append(expression, 1)
          } else {
            expression <- append(expression, -1)
          }
        }
      }
    } else {
      if (length(terms_list[[1]]$Description) > 0) {
        terms_list <- terms_list[[1]]
        for (term in conserved_terms) {
          temp_list <- terms_list %>% filter(Description == term)
          if (length(temp_list$Description) == 0) {
            expression <- append(expression, "NO")
          } else {
            expression <- append(expression, 1)
          }
        }
      } else {
        terms_list <- terms_list[[2]]
        for (term in conserved_terms) {
          temp_list <- terms_list %>% filter(Description == term)
          if (length(temp_list$Description) == 0) {
            expression <- append(expression, "NO")
          } else {
            expression <- append(expression, -1)
          }
        }
      }
    }
    
    
    
    expression <- as.numeric(unlist(expression))
    heatmap_df <- add_column(heatmap_df, !!paste(sp, tissue, sep = " ") := expression)
    
    
    
  }
}

annotation <- data_frame(Column = c("Arabidopsis thaliana leaf", "Arabidopsis thaliana root", "Triticum aestivum leaf", "Triticum aestivum root", "Solanum lycopersicum leaf", "Solanum lycopersicum root"), Organ = rep(c("Leaf", "Root"), 3)) %>% column_to_rownames("Column")
palette <- c("steelblue1", "navy", "lightsalmon")

heatmap_df[heatmap_df == "NO"] <- NA
heatmap <- pheatmap(heatmap_df, color = palette, na_col = "grey", legend_breaks = c(1, 0, -1), legend_labels = c("UP", "UP & DOWN", "DOWN"), border_color = "white", annotation_col = annotation, annotation_colors = list(Organ = c(Leaf = "darkgreen", Root = "brown")), cellheight = 10, cellwidth = 10, angle_col = "45", gaps_col = c(2, 4), cluster_rows = F, cluster_cols = F, labels_col = c("Arabidopsis", "Arabidopsis", "Wheat", "Wheat", "Tomato", "Tomato"))
pdf(here("Prova/heatmap_enrichment_organism.pdf"), width = 10, height = 25)
heatmap
dev.off()

