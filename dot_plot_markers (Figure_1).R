genes_of_interest <- marker_genes

for (gene in genes_of_interest) {
  plot <- FeaturePlot(ROOT, features = gene, label = FALSE) + NoLegend() + theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Bordo nero
    panel.background = element_blank(),                                # Rimuovi sfondo
    axis.line = element_line(color = "black")                          # Linee degli assi
  ) 
  
  ggsave(filename = paste0(gene, "_FeaturePlot.png"), plot = plot, width = 6, height = 6)
}

dot_data <- DotPlot(ROOT, features = marker_genes)$data


exp_mat <- dot_data %>%
  select(features.plot, id, avg.exp.scaled) %>%
  tidyr::pivot_wider(names_from = id, values_from = avg.exp.scaled, values_fill = 0) %>%
  column_to_rownames("features.plot") %>%
  as.matrix()


row_clust <- hclust(dist(exp_mat))
col_clust <- hclust(dist(t(exp_mat)))

a
gene_order <- rownames(exp_mat)[row_clust$order]
cluster_order <- colnames(exp_mat)[col_clust$order]


dot_data$features.plot <- factor(dot_data$features.plot, levels = gene_order)
dot_data$id <- factor(dot_data$id, levels = cluster_order)

filtered_data <- dot_data %>%
  filter(pct.exp >= 25)

tiff("Marker_genes_SAM_marker_order.tiff", units = "in", width = 15, height = 3.5, res = 450)
CairoSVG("Marker_genes_SAM_marker_order.svg", width = 15, height = 3.5)
postscript("aiuto.eps", width = 15, height = 3.5, horizontal = FALSE, onefile = FALSE, paper = "special")
ggplot(filtered_data, aes(x = id, y = features.plot, size = pct.exp, color = avg.exp.scaled)) +
  geom_point() +
  scale_color_gradient(low = "lightblue", high = "red4") +
  scale_size(range = c(1, 8)) +
  labs(
    title = "Bubble Plot",
    x = "Clusters",
    y = "Marker Genes",
    size = "% Cells Expressing",
    color = "Scaled Expression"
  )  +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),       # Rimuove la griglia
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Aggiunge una cornice nera
  )
dev.off()
