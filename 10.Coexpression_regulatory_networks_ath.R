library(tidyverse)
library(openxlsx)
library(here)
library(DESeq2)
library(sva)
library(BioNERO)
library(WGCNA)
library(enrichplot)
library(biomaRt)

orthology <- read.csv(here("Orthology analysis/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"), sep = "\t") %>%
  group_by(Orthogroup) %>%
  summarise(across(everything(), ~ paste(unique(.), collapse = ", "))) %>%
  mutate(Arabidopsis.thaliana = str_split(Arabidopsis.thaliana, ", ")) %>%
  rowwise() %>%
  unnest(Arabidopsis.thaliana) %>%
  filter(Arabidopsis.thaliana != "") %>%
  group_by(Orthogroup) %>%
  summarise(across(everything(), ~ paste(unique(.), collapse = ", ")))

mart <- useMart("plants_mart", dataset = "athaliana_eg_gene", host = "https://plants.ensembl.org")
connectome <- str_to_upper(str_split_1(read_file(here(paste0("Prova/Arabidopsis thaliana/Data/genes_connectome.txt"))), "\t"))
tf <- read_tsv(here("Prova/Arabidopsis thaliana/Data/Ath_TF_list.txt"))

orthogroups <- parse_orthofinder(here("Orthology analysis/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")) %>%
  mutate(Gene = str_to_upper(Gene))
species <- c("Arabidopsis thaliana", "Triticum aestivum", "Solanum lycopersicum")
generic_BP <- read_file(here("generic_BP.txt"))
generic_BP <- str_split_1(generic_BP, "\n")
gmt_token <- gprofiler2::upload_GMT_file(gmtfile = here("Prova/Arabidopsis thaliana/Data/athaliana.GO_BP.ENSG_PlantConnectomeDroughtAdded.gmt"))
plaza <- read.csv(here("Prova/Arabidopsis thaliana/Data/id_conversion.ath.csv"), skip = 8, row.names = NULL, sep = "\t")
gene_descr_1 <- plaza %>% filter(id_type == "symbol")
gene_descr <- plaza %>% 
  filter(id_type == "Alias") %>%
  group_by(X.gene_id) %>%
  summarise(id = paste(unique(id), collapse = ", ")) %>%
  add_row(gene_descr_1 %>% dplyr::select(-id_type)) %>%
  mutate(id = str_extract(id, "^[^,]+"))
gene_descr <- gene_descr %>%
  add_row(plaza %>% filter(!(X.gene_id %in% gene_descr$X.gene_id), id_type == "uniprot") %>% dplyr::select(-id_type)) %>%
  mutate(id = str_remove(id, regex("^at", ignore_case = TRUE)))
go <- clusterProfiler::read.gmt(here("Prova/Arabidopsis thaliana/Data/athaliana.GO_BP.ENSG.gmt")) %>%
  filter(term == "GO:0009415") %>%
  pull(gene)

set.seed(27000)
for (tissue in c("leaf", "root")) {
  
  
  degs_ath <- read.xlsx(here(paste0("Prova/Arabidopsis thaliana/Products/DEGs analysis/GDE_", tissue, ".xlsx")))
  
  conserved_degs <- read.xlsx(here(paste0("Prova/Orthology analysis/conserved_degs_", tissue, ".xlsx"))) %>%
    mutate(Arabidopsis.thaliana = str_split(Arabidopsis.thaliana, ", ")) %>%
    unnest(Arabidopsis.thaliana) %>%
    pull(Arabidopsis.thaliana)
  
  coldata <- read.xlsx(here(paste0("Prova/Arabidopsis thaliana/Data/metadata_filtered_2.xlsx"))) %>% filter(Major.Category == tissue)
  
  
  tpm <- read.table(here(paste0("Prova/Arabidopsis thaliana/Data/tpm_gene_level.txt"))) %>% 
    select(coldata$Run) %>% 
    filter(rownames(.) %in% degs_ath$Gene)
  
  
  final_exp <- exp_preprocess(
    tpm
  )
  
  power_ortho <- SFT_fit(final_exp)
  pdf(here(paste0("Prova/WGCNA/spearman/ath/", tissue, "_power_plot.pdf")), height = 4)
  print(power_ortho$plot)
  dev.off()
  
  net <- exp2gcn(final_exp, SFTpower = power_ortho$power)
  
  pdf(here(paste0("Prova/WGCNA/spearman/ath/", tissue, "_genes_per_module.pdf")))
  print(plot_ngenes_per_module(net))
  dev.off()
  
  
  hubs <- get_hubs_gcn(final_exp, net)
  
  results <- getBM(
    attributes = c("ensembl_gene_id", "description"),
    filters = "ensembl_gene_id",
    values = hubs$Gene,
    mart = mart
  )
  
  hubs <- hubs %>% left_join(results, by = join_by(Gene == ensembl_gene_id))
  write.xlsx(hubs, here(paste0("Prova/WGCNA/spearman/ath/hub_genes_", tissue, ".xlsx")))
  
  for (colour in unique(net$genes_and_modules$Modules)) {
    
    genes <- net$genes_and_modules %>% filter(Modules == colour) %>% pull(Genes)
    
    gost <- gprofiler2::gost(query = genes, organism = "athaliana", custom_bg = rownames(filt_tpm), highlight = F, evcodes = T)
    df <- data.frame(gost$result$term_id, gost$result$source, gost$result$term_name, gost$result$p_value, gost$result$query, gost$result$term_size, gost$result$query_size, gost$result$intersection_size, gost$result$intersection)
    if (length(df) > 0) {
      colnames(df) <- c("ID", "Source", "Description", "p.adjust", "Query", "Term_size", "Query_size", "Intersection_size", "geneID")
      df <- df %>%
        filter(Source == "GO:BP", !(ID %in% generic_BP))
      
      write.xlsx(df, here(paste0("Prova/WGCNA/spearman/ath/enrichment_analysis_", tissue, "_", colour, ".xlsx")))
      
      # vedere i moduli stabili, fare arricchimento su tutti i moduli (conservati e non), fare tabella unica dell'arricchimento di tutti i moduli
      
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
                           font.size = 9,
                           orderBy = "GeneRatio",
                           label_format = 50) + 
          scale_fill_gradientn(
            colors = c("red", "blue"),
            trans = "log10",  # log scale to spread values more evenly
            name = "p.adjust"
          ) 
        
        ggsave(here(paste0("Prova/WGCNA/spearman/ath/enrichment_", tissue, "_", colour, ".pdf")), dotplot, height = 7, width = 8, units = "in", dpi = 600)
      }
      
      
    }
  }
  
  if (tissue == "leaf") {
    conserved_modules <- c("black")
  } else {
    conserved_modules <- c("darkturquoise")
  }
  
  for (colour in conserved_modules) {
  
    edges_filtered_sft <- get_edge_list(net, module = colour, filter = T)
    
    write.xlsx(edges_filtered_sft, here(paste0("Prova/WGCNA/spearman/ath/", colour, "_module_", tissue, ".xlsx")))
    
    
    hub_genes <- unique(hubs$Gene)
    
    edges_hubs <- edges_filtered_sft %>%
      filter(map_lgl(Var1, ~ length(intersect(.x, hub_genes)) > 0))
    
    if (length(edges_hubs$Var1) > 0) {
      results <- getBM(
        attributes = c("ensembl_gene_id", "description"),
        filters = "ensembl_gene_id",
        values = edges_hubs %>% pull(Var2),
        mart = mart
      )
      
      edges_hubs <- edges_hubs %>% 
        left_join(results, by = join_by(Var2 == ensembl_gene_id))
      
      write.xlsx(edges_hubs, here(paste0("Prova/WGCNA/spearman/ath/edges_hubs_", colour, "_", tissue, ".xlsx")))
      
      
      edges_connectome <- edges_filtered_sft %>%
        filter(map_lgl(Var1, ~ length(intersect(.x, connectome)) > 0))
      
      if (length(edges_connectome$Var1) > 0) {
        results <- getBM(
          attributes = c("ensembl_gene_id", "description"),
          filters = "ensembl_gene_id",
          values = edges_connectome %>% pull(Var2),
          mart = mart
        )
        
        edges_connectome <- edges_connectome %>% 
          left_join(results, by = join_by(Var2 == ensembl_gene_id))
        
        write.xlsx(edges_connectome, here(paste0("Prova/WGCNA/spearman/ath/edges_connectome_", colour, "_", tissue, ".xlsx")))
        
        tf_set        <- unique(tf$Gene_ID)
        conn_set      <- unique(connectome)
        hub_set       <- unique(hub_genes)
        degs_over_set <- degs_ath %>% filter(expression == "OVER") %>% pull(Gene)
        degs_under_set<- degs_ath %>% filter(expression == "UNDER") %>% pull(Gene)
        
        edges_filtered_sft <- edges_filtered_sft %>%
          pivot_longer(
            cols = c(Var1, Var2),
            names_to = ".value",
            names_pattern = "([A-Za-z]+)\\d"
          ) %>%
          select(-Freq) %>% 
          mutate(gene_flags = map(
                   Var,
                   ~ list(
                     tf   = any(.x %in% tf_set),
                     conn = any(.x %in% conn_set),
                     hub  = any(.x %in% hub_set),
                     over = any(.x %in% degs_over_set),
                     under = any(.x %in% degs_under_set),
                     conserved = (any(.x %in% conserved_degs)),
                     go = (any(.x %in% go))
                   )
                 ),
                 transcription.factor = if_else(map_lgl(gene_flags, "tf"), "Yes", "No"),
                 plant.connectome     = if_else(map_lgl(gene_flags, "conn"), "Yes", "No"),
                 hub                  = if_else(map_lgl(gene_flags, "hub"), "Yes", "No"),
                 conserved            = if_else(map_lgl(gene_flags, "conserved"), "Yes", "No"),
                 go                   = if_else(map_lgl(gene_flags, "go"), "Yes", "No"),
                 degs = case_when(
                   map_lgl(gene_flags, "over")  ~ "OVER",
                   map_lgl(gene_flags, "under") ~ "UNDER",
                   TRUE                         ~ NA_character_
                 ),
                 connectome_degs = paste(plant.connectome, degs)
          ) %>%
          select(-gene_flags) %>% 
          distinct(Var, .keep_all = TRUE) %>%
          left_join(hubs, by = join_by(Var == Gene)) %>%
          select(-c("Module", "description")) %>% 
          mutate(kWithin = tidyr::replace_na(kWithin, 0)) %>%
          left_join(gene_descr, by = join_by(Var == X.gene_id))
        
        write.xlsx(edges_filtered_sft, here(paste0("Prova/WGCNA/spearman/ath/genes_annotation_", colour, "_", tissue, ".xlsx")))
        
      }
    }
    
  }
}
