# scRNAseq-hypoxia-SAM-RAM

This repository contains all scripts used for the analyses presented in the manuscript:  
**"Inferring hypoxia-responsive regulators of cell fate in plant meristems through single-cell transcriptomics"**.

The repository is organized into **R** and **Python** scripts that reproduce, step-by-step, the analyses performed on single-cell RNA-seq datasets of *Arabidopsis thaliana* shoot and root apices, with a focus on hypoxia responses and stem cell regulation.

---

## Repository Structure

The scripts are organized according to the main analyses performed in the study:

- **SAM (Figure 1)**
  - Preprocessing and clustering of shoot apex datasets.
  - Identification of major cell populations and marker genes.

- **RAM (Figure 1)**
  - Preprocessing and clustering of root apex datasets.
  - Identification of major cell populations and marker genes.

- **dot_plot_markers (Figure 1)**
  - Dot plot visualization of selected cell-type marker genes.

- **gene_modules (Figure 3)**
  - Pseudotime analysis and identification of gene modules in shoot and root apices.
  - Functional enrichment analysis of gene modules.

- **heatmap_overlap (Figure 2 and Supplemental Figures)**
  - Overlap analysis between stem cell markers (SCMs) and hypoxia-responsive genes (HRGs).
  - Generation of heatmaps comparing co-expression patterns across tissues.

- **monocle2 / monocle3 (Figure 3, 4)**
  - Pseudotime trajectory inference and visualization of gene dynamics.

- **RNA_velocity (Figure 4)**
  - RNA velocity analysis of stem cell and differentiated populations to study transcriptional dynamics.

---

## How to Use

Each script is self-contained and includes the commands necessary to reproduce the corresponding figures and analyses.  
To execute the analyses:

1. Clone the repository:  
   ```bash
   git clone https://github.com/SimoneFromThePlantLab/scRNAseq-hypoxia-SAM-RAM.git
   ```
2. Open the desired script in **R** or **Python**.
3. Follow the script sequentially to reproduce preprocessing, analysis, and figure generation.

---

## Requirements

- **R (≥ 4.0)** with packages:
  - `Seurat`, `ggplot2`, `dplyr`, `ComplexHeatmap`, `Cairo`, `tidyr`
- **Python (≥ 3.8)** with packages:
  - `scanpy`, `scvelo`, `numpy`, `pandas`, `matplotlib`

---

## Notes

- The scripts are provided for transparency and reproducibility.
- Input data are not included in the repository but can be accessed from [link to data repository, if available] or requested from the corresponding author.
- File names and figure references correspond to those reported in the manuscript for easy navigation.
