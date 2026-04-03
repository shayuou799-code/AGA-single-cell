# ============================================================================
# NicheNet Ligand-Target Prediction Analysis
# Project: AGA single-cell RNA-seq (Nature Communications)
# Input:  Annotated Seurat object with level3_identity and Type (V/O)
# Output: Ligand activity predictions, ligand-target heatmaps
# Example: Mast cell -> HFSC intercellular signaling
# ============================================================================

set.seed(42)

library(nichenetr)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggrepel)

# ============================================================================
# Part 1: Load NicheNet Pre-trained Networks (Human)
# ============================================================================

# --- 1. Load human ligand-target matrix ---
ligand_target_matrix <- readRDS(url(
  "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))

# MODIFIED: using human lr_network database
lr_network_human <- readRDS(url(
  "https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
lr_network <- lr_network_human %>%
  dplyr::rename(ligand = from, receptor = to) %>%
  distinct(ligand, receptor)

weighted_networks <- readRDS(url(
  "https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
weighted_networks_lr <- weighted_networks$lr_sig %>%
  inner_join(lr_network %>% distinct(ligand, receptor),
             by = c("from" = "ligand", "to" = "receptor"))

# ============================================================================
# Part 2: Prepare Sender and Receiver Gene Sets
# ============================================================================

# --- 2. Set cell type identities ---
# seuratObj: AGA Seurat object with level3_identity and Type metadata
Idents(seuratObj) <- seuratObj$level3_identity

# Helper: extract expressed genes (>5% of cells)
get_expressed_genes <- function(celltype, seurat_obj, pct = 0.05) {
  cells <- WhichCells(seurat_obj, idents = celltype)
  expr_mat <- GetAssayData(seurat_obj, slot = "data")[, cells]
  n_cells <- ncol(expr_mat)
  rownames(expr_mat)[rowSums(expr_mat > 0) > (n_cells * pct)]
}

# --- 3. Define sender (Mast) and receiver (HFSC) ---
sender_celltypes <- c("MAST_1", "MAST_2")
receiver_celltype <- c("HFSC1", "HFSC2")

expressed_genes_sender <- sender_celltypes %>%
  lapply(get_expressed_genes, seuratObj, 0.05) %>%
  unlist() %>% unique()

expressed_genes_receiver <- receiver_celltype %>%
  lapply(get_expressed_genes, seuratObj, 0.05) %>%
  unlist() %>% unique()

# --- 4. Define Vertex-upregulated genes in receiver (HFSC) ---
seurat_hfsc <- subset(seuratObj, idents = receiver_celltype)
Idents(seurat_hfsc) <- seurat_hfsc$Type

de_hfsc <- FindMarkers(seurat_hfsc,
                       ident.1 = "Vertex",
                       ident.2 = "Occipital",
                       min.pct = 0.1) %>%
  rownames_to_column("gene")

# Filter to significant upregulated genes as target set
genes_of_interest_hfsc <- de_hfsc %>%
  filter(avg_log2FC > 0.25 & p_val_adj < 0.05) %>%
  pull(gene) %>%
  .[. %in% rownames(ligand_target_matrix)]

# Background: all expressed genes in receiver
background_expressed_genes <- expressed_genes_receiver %>%
  .[. %in% rownames(ligand_target_matrix)]

# ============================================================================
# Part 3: Identify Potential Ligands and Predict Activity
# ============================================================================

# --- 5. Filter ligands with cognate receptors ---
ligands <- lr_network %>% pull(ligand) %>% unique()
receptors <- lr_network %>% pull(receptor) %>% unique()

expressed_ligands <- intersect(ligands, expressed_genes_sender)
expressed_receptors <- intersect(receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>%
  filter(ligand %in% expressed_ligands & receptor %in% expressed_receptors) %>%
  pull(ligand) %>% unique()

# --- 6. Predict ligand activity (Mast -> HFSC) ---
ligand_activities <- predict_ligand_activities(
  geneset = genes_of_interest_hfsc,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)

# Top 20 ligands by Pearson correlation
best_upstream_ligands <- ligand_activities %>%
  top_n(20, pearson) %>%
  arrange(-pearson) %>%
  pull(test_ligand)

print("Top 20 Mast -> HFSC ligands:")
print(best_upstream_ligands)

# ============================================================================
# Part 4: Ligand-Target Visualization
# ============================================================================

# --- 7. Extract top target genes for best ligands ---
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = genes_of_interest_hfsc,
         ligand_target_matrix = ligand_target_matrix,
         n = 200) %>%
  bind_rows()

# --- 8. Ligand-target heatmap (top 5 ligands, top 20 targets) ---
top_ligands_show <- best_upstream_ligands[1:5]
targets_to_show <- active_ligand_target_links_df %>%
  filter(ligand %in% top_ligands_show) %>%
  top_n(20, weight) %>%
  pull(target) %>% unique()

p_ligand_target <- make_heatmap_ggplot(
  ligand_target_matrix[targets_to_show, top_ligands_show],
  y_name = "Prioritized Ligands",
  x_name = "Predicted Target Genes",
  legend_title = "Regulatory Potential"
)
print(p_ligand_target)

# --- 9. Ligand activity bar plot ---
ligand_aupr_matrix <- ligand_activities %>%
  filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>%
  select(aupr_corrected) %>%
  arrange(aupr_corrected) %>%
  as.matrix()

p_ligand_aupr <- make_heatmap_ggplot(
  ligand_aupr_matrix,
  y_name = "Prioritized Ligands",
  x_name = "",
  legend_title = "AUPR"
)
print(p_ligand_aupr)

# --- 10. Validate ligand expression (DotPlot) ---
DotPlot(seuratObj,
        features = best_upstream_ligands[1:10],
        idents = c(sender_celltypes, receiver_celltype),
        cols = c("lightgrey", "red")) +
  RotatedAxis() +
  ggtitle("Top Mast -> HFSC Ligand Expression")
