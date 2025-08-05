
library(monocle3)
library(SeuratWrappers)
library(dplyr)
library(tibble)


seurat_obj <- JoinLayers(seurat_obj)

cds <- as.cell_data_set(seurat_obj)
cds <- estimate_size_factors(cds)

cds <- preprocess_cds(cds, num_dim = 50)  
cds <- reduce_dimension(cds, reduction_method="UMAP") 
cds <- cluster_cells(cds, resolution = 1e-5) 
cds <- learn_graph(cds) 
plot_cells(cds, color_cells_by = "cluster", cell_size = 0.5, group_label_size = 10)

rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

tiff("seurat_obj_STM.tiff", units="in", width=6, height=5, res=450)
plot_cells(cds, genes=c("ATXGXXXX"), cell_size = 1) + theme(
  panel.border = element_rect(color = "black", fill = NA, size = 1),  
  panel.background = element_blank(),                               
  axis.line = element_line(color = "black")                          
) 
dev.off()

#PSEUDOTIME
cds <- order_cells(cds) 

tiff("seurat_obj_pesudotime.tiff", units="in", width=6, height=5, res=450)
plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE, cell_size = 0.5) + theme(
  panel.border = element_rect(color = "black", fill = NA, size = 1),  
  panel.background = element_blank(),                                
  axis.line = element_line(color = "black")                          
) 
dev.off()

my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("ATXGXXXX")))
cds_subset <- cds[my_genes,]

tiff("seurat_obj_pesudotime_STM.tiff", units="in", width=5, height=5, res=450)
plot_genes_in_pseudotime(cds_subset) + theme(
  panel.border = element_rect(color = "black", fill = NA, size = 1), 
  panel.background = element_blank(),                                
  axis.line = element_line(color = "black")                          
) + NoLegend() 
dev.off()