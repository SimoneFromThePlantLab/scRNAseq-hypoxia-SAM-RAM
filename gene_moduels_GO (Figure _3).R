
pseudotime_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 6)

pr_deg_ids <- row.names(subset(pseudotime_genes, q_value < 0.05))

gene_module_df <- find_gene_modules(cds[pr_deg_ids, ], resolution = 1e-3)

valid_cells <- colnames(cds)[is.finite(pseudotime(cds))]
cds <- cds[, valid_cells]

cell_group_df <- tibble::tibble(
  cell = row.names(colData(cds)), 
  cell_group = cut(pseudotime(cds), breaks = 10, labels = paste0("Bin_", 1:10)) 
)

agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)

row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Pseudotime ", colnames(agg_mat))

tiff("heatmap_pseudotime_new.tiff", units="in", width=8, height=5, res=450)
pheatmap::pheatmap(
  agg_mat, 
  cluster_rows = TRUE, 
  cluster_cols = FALSE,          # Nessun clustering delle colonne per rispettare l'ordine dello pseudotime
  scale = "column",              # Normalizzazione per colonna
  clustering_method = "ward.D2", # Metodo di clustering per righe
  fontsize = 6,
  main = "Gene Modules Ordered by Pseudotime"
)
dev.off()

modules_list <- split(gene_module_df$id, gene_module_df$module)

output_dir <- "gene_modules"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

for (module in names(modules_list)) {
  module_genes <- modules_list[[module]]
  write.csv(module_genes, file = file.path(output_dir, paste0("module_", module, ".csv")), row.names = FALSE)
}

str(gene_module_df)
head(gene_module_df)

##################################################
# Specify "early phase" and "latephase" modules
##################################################

early_modules <- c(x1...x2)
late_modules <- c(x1...x2)

early_genes <- gene_module_df %>%
  filter(module %in% early_modules) %>%
  pull(id)  

late_genes <- gene_module_df %>%
  filter(module %in% late_modules) %>%
  pull(id)

early_genes <- unique(early_genes)

write.csv(early_genes, "early_phase_genes.csv", row.names = FALSE)
write.csv(late_genes, "late_phase_genes.csv", row.names = FALSE)

library(clusterProfiler)
library(org.At.tair.db)

early_go_results <- enrichGO(
  gene         = early_genes,
  OrgDb        = org.At.tair.db,
  keyType      = "TAIR",
  ont          = "BP",   # Ontologia biologica
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

tiff("barplot_GO.tiff", units="in", width=10, height=7, res=450)
barplot(early_go_results, showCategory = 15)
dev.off()

enrich_result <- enrichGO(
  gene = early_genes,
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,    
  minGSSize = 5,          
  maxGSSize = 200         
)

barplot(enrich_result, showCategory = 10)

#Enrichment for shoot and root
early_genes_shoot <- early_phase_genes$x
early_genes_root <- early_phase_genes_root$x

late_genes_shoot <- late_phase_genes$x
late_genes_root <- late_phase_genes_root$x

get_go_enrichment <- function(gene_list) {
  enrichGO(
    gene         = gene_list,
    OrgDb        = org.At.tair.db,
    keyType      = "TAIR",
    ont          = "BP",   
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable    = TRUE
  ) %>% as.data.frame()
}

go_shoot <- get_go_enrichment(early_genes_shoot) %>%
  arrange(p.adjust) %>%  
  head(20)               
go_shoot$Tissue <- "Shoot Apex"

go_root <- get_go_enrichment(early_genes_root) %>%
  arrange(p.adjust) %>%
  head(20)
go_root$Tissue <- "Root Tip"

common_go_terms <- intersect(go_shoot$Description, go_root$Description)

go_shoot_common <- go_shoot %>% filter(Description %in% common_go_terms)
go_root_common <- go_root %>% filter(Description %in% common_go_terms)

go_combined <- bind_rows(go_shoot_common, go_root_common)

go_combined$Tissue <- factor(go_combined$Tissue, levels = c("Shoot Apex", "Root Tip"))

CairoSVG("GO_plot.svg", width = 10, height = 8)
ggplot(go_combined, aes(x = Tissue, y = reorder(Description, desc(Description)), size = Count, color = Tissue)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(5, 20)) + 
  labs(
    title = "Common GO Terms in Early Modules",
    x = "Organ",
    y = "GO Term",
    size = "Gene Count",
    color = "Tissue"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12)) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    panel.background = element_blank(),                                
    axis.line = element_line(color = "black")                          
  ) 
dev.off()

CairoSVG("GO_plot.svg", width = 10, height = 8)
ggplot(go_combined, aes(x = Tissue, y = fct_rev(factor(Description)), size = Count, color = Tissue)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(5, 20)) +
  labs(
    title = "Common GO Terms in Early Modules",
    x = "Organ",
    y = "GO Term",
    size = "Gene Count",
    color = "Tissue"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12)) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")
  )
dev.off()

go_output <- go_combined %>%
  dplyr::select(Description, Count, p.adjust, Tissue)

write.csv(go_output, "GO_enrichment_results.csv", row.names = FALSE)

go_terms_of_interest <- c(
  "cellular response to decreased oxygen levels",
  "cellular response to hypoxia",
  "cellular response to oxygen levels",
  "response to decreased oxygen levels",
  "response to hypoxia",
  "response to oxidative stress",
  "response to oxygen levels",
  "response to water",
  "response to water deprivation",
  "response to wounding",
  "secondary metabolic process",
  "sulfur compound metabolic process"
)

extract_genes_for_term <- function(go_df, tissue_label) {
  go_df %>%
    dplyr::filter(Description %in% go_terms_of_interest) %>%
    dplyr::select(Description, geneID) %>%
    dplyr::mutate(
      Tissue = tissue_label,
      Genes = strsplit(geneID, "/")
    ) %>%
    dplyr::select(-geneID) %>%
    tidyr::unnest(Genes)
}

go_shoot_genes <- extract_genes_for_term(go_shoot, "Shoot Apex")
go_root_genes <- extract_genes_for_term(go_root, "Root Tip")

all_go_genes <- bind_rows(go_shoot_genes, go_root_genes)

write.csv(all_go_genes, "GO_genes_common_terms_gene_list.csv", row.names = FALSE)