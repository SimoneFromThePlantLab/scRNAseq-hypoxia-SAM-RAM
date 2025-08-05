
##########################
##FIGURE 2 HEATMAP
##########################


genes_of_interest <- gene_list$GeneID

for (gene in genes_of_interest) {
  # Identifica le cellule che esprimono il gene
  cells_expressing <- colnames(GetAssayData(seurat_obj, slot = "counts"))[
    GetAssayData(seurat_obj, slot = "counts")[gene, ] > 0
  ]
  
  # Salva il file CSV con le cellule trovate
  write.csv(cells_expressing, paste0("cells_expressing_", gene, ".csv"), row.names = FALSE)
}

# Directory dei file marker e ipossici
marker_dir <- "~/path/to/Cells_expressing/"   # Modifica con il percorso corretto # QUI PER GRUPPI DI GENI
hypoxia_dir <- "~path/to/Cells_expressing/" # Modifica con il percorso corretto

# Ottieni i file marker e ipossici
marker_files <- list.files(marker_dir, pattern = "*.csv", full.names = TRUE)
hypoxia_files <- list.files(hypoxia_dir, pattern = "*.csv", full.names = TRUE)

# Inizializza un data frame per i risultati
results <- data.frame(
  Marker_Gene = character(),        # Nome del gene marker
  Hypoxia_Gene = character(),       # Nome del gene ipossico
  Common_Cells = integer(),         # Numero di cellule in comune
  Total_Cells = integer(),          # Totale delle cellule nel file ipossico
  Percentage_Common = numeric()     # Percentuale di cellule in comune
)

# Loop per confrontare ogni marker con ogni ipossico
for (marker_file in marker_files) {
  # Leggi il file marker
  marker_cells <- read.csv(marker_file, header = TRUE, stringsAsFactors = FALSE)$x
  
  # Loop sui file ipossici
  for (hypoxia_file in hypoxia_files) {
    # Leggi il file ipossico
    hypoxia_cells <- read.csv(hypoxia_file, header = TRUE, stringsAsFactors = FALSE)$x
    
    # Trova le cellule in comune
    common_cells <- intersect(marker_cells, hypoxia_cells)
    
    # Calcola il totale e la percentuale
    total_cells <- length(hypoxia_cells)
    percentage_common <- if (total_cells > 0) {
      (length(common_cells) / total_cells) * 100
    } else {
      NA
    }
    
    # Aggiungi i risultati al data frame
    results <- rbind(results, data.frame(
      Marker_Gene = tools::file_path_sans_ext(basename(marker_file)),   # Nome del file marker
      Hypoxia_Gene = tools::file_path_sans_ext(basename(hypoxia_file)), # Nome del file ipossico
      Common_Cells = length(common_cells),
      Total_Cells = total_cells,
      Percentage_Common = percentage_common
    ))
  }
}

library(reshape2)

# Converti i risultati in una matrice per la heatmap
heatmap_data <- reshape2::acast(results, Marker_Gene ~ Hypoxia_Gene, value.var = "Percentage_Common", fill = 0)

# Filtra i dati per includere solo le combinazioni con almeno il 10% di overlap
filtered_results <- subset(results, Percentage_Common >= 20)

# Converti i dati filtrati in una matrice per la heatmap
filtered_heatmap_data <- reshape2::acast(filtered_results, Marker_Gene ~ Hypoxia_Gene, value.var = "Percentage_Common", fill = 0)

# Heatmap con miglioramenti grafici
tiff("Heatmap_core_xylem.tiff", units = "in", width = 8, height = 4, res = 450)
ggplot(melt(heatmap_data), aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "black") +  # Cornice nera per ogni cella
  scale_fill_gradient(low = "white", high = "red", na.value = "gray") +
  theme_minimal(base_size = 12) +  # Cambia lo stile base e aumenta la dimensione del testo
  labs(x = "Hypoxia Genes", y = "Marker Genes", fill = "Percentage") +
  theme(
    panel.grid.major = element_line(color = "gray80"),  # Linee di griglia
    panel.grid.minor = element_blank(),  # Rimuove griglia secondaria
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Cornice nera
    axis.text.x = element_text(angle = 90, hjust = 1),  # Ruota le etichette delle colonne
    axis.text.y = element_text(hjust = 1)  # Migliora leggibilità delle righe
  )
dev.off()


####################################
##SUPPLEMENTAL FIGURE 2 HEATMAP
####################################

genes_of_interest <- gene_list$GeneID  

counts_mat <- GetAssayData(seurat_obj, slot = "counts")

genes_present <- genes_of_interest[genes_of_interest %in% rownames(counts_mat)]
counts_sub <- counts_mat[genes_present, ]

#Create binary matrix (1 = expressed, 0 = not expressed)
binary_mat <- counts_sub > 0

marker_count_per_cell <- colSums(binary_mat)

cells_expressing_10 <- names(marker_count_per_cell[marker_count_per_cell >= 10])

write.csv(cells_expressing_10, "cells_expressing_10_markers.csv", row.names = FALSE)

library(pheatmap)

HRG_genes <- gene_list$GeneID  # prima colonna del tuo file
HRG_genes <- HRG_genes[HRG_genes %in% rownames(seurat_obj)]  # tieni solo geni presenti

counts_mat <- GetAssayData(seurat_obj, slot = "counts")

counts_HRG <- counts_mat[HRG_genes, cells_expressing_10, drop = FALSE]

binary_HRG <- counts_HRG > 0

percent_expressing <- rowSums(binary_HRG) / ncol(binary_HRG) * 100

heatmap_matrix <- matrix(percent_expressing, nrow = length(percent_expressing))
rownames(heatmap_matrix) <- HRG_genes
colnames(heatmap_matrix) <- "Percent expressing"

heatmap_matrix_t <- t(heatmap_matrix)

gene_labels <- setNames(gene_list$Gene_name, gene_list$GeneID)

colnames(heatmap_matrix_t) <- gene_labels[colnames(heatmap_matrix_t)]

pheatmap(
  heatmap_matrix_t,
  color = colorRampPalette(c("white", "red"))(100),
  cluster_rows = FALSE,   # ora le righe sono solo 1 (percentuali)
  cluster_cols = TRUE,    # puoi clusterizzare i geni se vuoi
  main = "HRG expression in stem-like cells (≥10 SCM)"
)
