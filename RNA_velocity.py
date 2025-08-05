
#questo in R
# 1. Salva i barcodes delle cellule da tenere
write.csv(Cells(ROOT_subset), file = "barcodes_subset.csv", row.names = FALSE)

# 2. Salva le coordinate UMAP
write.csv(Embeddings(ROOT_subset, "umap"), file = "umap_subset.csv")

# 3. Salva i cluster
clusters_df <- data.frame(cluster = Idents(ROOT_subset))
write.csv(clusters_df, file = "clusters_subset.csv")

#poi Python
import scanpy as sc
import scvelo as scv
import pandas as pd

adata = sc.read('root_2.loom', cache=True)

barcodes = pd.read_csv("barcodes_subset.csv")
selected_barcodes = barcodes.iloc[:, 0].values
adata.obs_names = [bc.replace('-1', '') for bc in adata.obs_names]

print(selected_barcodes[:5])
print(adata.obs_names[:5])

cleaned_obs_names = [bc.split(':')[1].replace('x', '').replace('t', '') for bc in adata.obs_names]
valid_barcodes = [bc for bc in selected_barcodes if bc.replace('-1', '') in cleaned_obs_names]
mask = [bc.split(':')[1].replace('x','').replace('t','') in [v.replace('-1','') for v in selected_barcodes]
         for bc in adata.obs_names]

adata_subset = adata[mask]

print(adata_subset.shape)

# Pulisci i nomi in adata per confrontarli
adata_cleaned_barcodes = [x.split(":")[1].replace('x','').replace('t','') for x in adata.obs_names]

# Pulisci anche i selected barcodes
selected_barcodes_clean = [x.replace("-1", "") for x in selected_barcodes]

# Trova gli indici validi
mask = [bc in selected_barcodes_clean for bc in adata_cleaned_barcodes]

# Subset
adata_subset = adata[mask].copy()

umap = pd.read_csv("umap_subset.csv", index_col=0)
#adata_subset = adata_subset[umap.index]  # riordina in base all'UMAP <----------------------------------------

print(umap.index[:5])
print(adata_subset.obs_names[:5])
adata_cleaned_barcodes = [x.split(":")[1].replace('x','').replace('t','') for x in adata_subset.obs_names]
umap_cleaned_index = umap.index.str.replace("-1", "")

common = list(set(adata_cleaned_barcodes) & set(umap_cleaned_index))
print(f"Matched cells: {len(common)}")

barcode_to_index = {bc.split(":")[1].replace('x','').replace('t',''): i for i, bc in enumerate(adata_subset.obs_names)}
ordered_indices = [barcode_to_index[bc] for bc in umap_cleaned_index if bc in barcode_to_index]
adata_subset = adata_subset[ordered_indices].copy()
adata_subset.obsm["X_umap"] = umap.values
print(adata_subset)

# Carica le annotazioni di cluster da file o da un altro oggetto
cluster_df = pd.read_csv("clusters_subset.csv", index_col=0)  # deve avere index = cell barcodes

# Rimuovi '-1' o aggiustali per farli combaciare
cluster_df.index = cluster_df.index.str.replace("-1", "")

print(cluster_df.index[:5])
print(adata_subset.obs_names[:5])

# Aggiusta i barcodes di cluster_df per farli combaciare con quelli di adata_subset
cluster_df.index = ["root_2:" + bc + "x" for bc in cluster_df.index]
adata_subset.obs["clusters"] = cluster_df.loc[adata_subset.obs_names, "x"]

scv.tl.velocity(adata_subset)
scv.tl.velocity_graph(adata_subset, n_jobs=16, show_progress_bar=False)

scv.pl.velocity_embedding_stream(adata_subset, basis="umap", color="clusters", dpi=600)

######################
## FIN QUI PERFETTO
######################

###################################################
## IMPORTANTISSIMO
## SE VUOI UTILIZZARE I CODICI AGI E NON I NOMI DEI GENI USA LA LINA QUI SOTTO

adata_subset.var_names = adata_subset.var["Accession"]
adata_subset.var_names_make_unique()

###################################################

# 1. Recupero dei parametri dinamici e calcolo del latent time
scv.tl.recover_dynamics(adata_subset)
scv.tl.latent_time(adata_subset)

# Visualizzazione del latent time
scv.pl.scatter(adata_subset, color='latent_time', color_map='viridis', save='_latent_time.png')

# 2. Calcolo della confidenza nei vettori di velocity
scv.tl.velocity_confidence(adata_subset)

# Visualizzazione della confidenza
scv.pl.scatter(adata_subset, color=['velocity_length', 'velocity_confidence'], save='_velocity_confidence.png')

# 3. Identificazione dei geni dinamici più rilevanti
scv.tl.rank_dynamical_genes(adata_subset, groupby='clusters')
top_genes = adata_subset.uns['rank_dynamical_genes']['names']
top_genes_df = pd.DataFrame(top_genes)

# Selezione dei primi 5 geni dinamici nel primo cluster
genes_to_plot = top_genes_df.iloc[:5, 0].tolist()

# 4. Visualizzazione dei profili dinamici
scv.pl.dynamic(adata_subset, var_names=genes_to_plot, ncols=2, save='_dynamic_genes.png')

#################
## OK ANCHE QUI
#################

scv.settings.verbosity = 3  # per vedere i log
scv.settings.presenter_view = True

# 1. Calcolo della pseudotime velocity-based
scv.tl.velocity_pseudotime(adata_subset)

# 2. Identificazione degli stati terminali
scv.tl.terminal_states(adata_subset)

# 3. Visualizza i terminal states sulla UMAP
scv.pl.scatter(adata_subset, color=['velocity_pseudotime', 'terminal_states'], cmap='viridis')

scv.tl.velocity_confidence(adata_subset)
scv.pl.scatter(adata_subset, color='velocity_confidence', cmap='coolwarm')

scv.tl.recover_dynamics(adata_subset)
scv.tl.velocity(adata_subset, mode='dynamical')
scv.tl.velocity_graph(adata_subset)
scv.tl.latent_time(adata_subset)
scv.pl.scatter(adata_subset, color='latent_time', cmap='viridis')

scv.tl.rank_velocity_genes(adata_subset, groupby='clusters', min_corr=.3)
df = pd.DataFrame(adata_subset.uns['rank_velocity_genes']['names'])
df.head()

http://localhost:8889/tree?token=34be8082868e69fb4a652b4584b6a7c17adc5d5fb3d6dd92