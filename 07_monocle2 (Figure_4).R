
library(monocle)

seurat_obj <- readRDS("seurat_obj_subset.rds")

expression_matrix <- as(as.matrix(GetAssayData(seurat_obj, slot = "counts")), "sparseMatrix")

cell_metadata <- seurat_obj@meta.data

gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix), row.names = rownames(expression_matrix))

fd <- new("AnnotatedDataFrame", data = gene_annotation)

pd <- new("AnnotatedDataFrame", data = cell_metadata)

cds <- newCellDataSet(expression_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.1,
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

expr_matrix <- as.matrix(exprs(cds))

variable_genes <- order(matrixStats::rowVars(expr_matrix), decreasing = TRUE)[1:1500]
cds_subset <- cds[variable_genes, ]

cds_subset <- reduceDimension(cds_subset, method = "DDRTree", max_components = 3)

cds_subset <- orderCells(cds_subset)

plot_cell_trajectory(cds_subset, color_by = "State") + 
  scale_color_viridis_c(option = "plasma") + 
  theme_void()

p1

ggsave(file="test.svg", plot=p1, width=10, height=8, dpi = 600)

plot_genes_in_pseudotime(cds_subset["ATXGXXXX"])

pData(cds_subset)$gene_expression <- exprs(cds_subset)["ATXGXXXX", ]

plot_cell_trajectory(cds_subset, color_by = "gene_expression") +
  scale_color_viridis_c(option = "plasma") +  
  theme_minimal()


gene_name <- "ATXGXXXX"  
pData(cds_state)$gene_expression <- exprs(cds_state)[gene_name, ]

CairoSVG("gene.svg", width = 6, height = 10)
plot_cell_trajectory(cds_state, color_by = "gene_expression", theta = 270, cell_size = 3) +
  scale_color_gradientn(colors = c("grey88", "yellow", "orange", "red", "red4"), 
                        values = c(0, 0.2, 0.3, 0.7, 1), guide = "none") + 
  theme_void()
dev.off()
