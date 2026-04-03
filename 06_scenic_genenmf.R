# ============================================================================
# Gene Regulatory Network (pySCENIC) and GeneNMF Meta-program Analysis
# Project: AGA single-cell RNA-seq (Nature Communications)
# Input:  Seurat/AnnData object, pySCENIC outputs
# Output: Regulon AUCell scores, GeneNMF meta-programs, heatmaps
# ============================================================================

set.seed(42)

# ============================================================================
# Part 1: pySCENIC Command-line Workflow (run in terminal)
# ============================================================================

# --- Step 1 (Python): Export expression matrix to loom format ---
# import scanpy as sc
# import numpy as np
# import loompy
#
# sc.pp.filter_genes(adata, min_cells=3)
# row_attrs = {"Gene": np.array(adata.var_names)}
# col_attrs = {"CellID": np.array(adata.obs_names)}
# lp = adata.to_df().T
# loompy.create("scenic_input.loom", lp.values, row_attrs, col_attrs)

# --- Step 2 (CLI): GRN inference ---
# pyscenic grn scenic_input.loom \
#     hs_hgnc_tfs.txt \
#     -o adjacencies.tsv \
#     --num_workers 16

# --- Step 3 (CLI): Regulon prediction (cisTarget) ---
# pyscenic ctx adjacencies.tsv \
#     hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
#     --annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
#     --expression_mtx_fname scenic_input.loom \
#     --mode "dask_multiprocessing" \
#     --output regulons.csv \
#     --num_workers 16

# --- Step 4 (CLI): AUCell scoring ---
# pyscenic aucell scenic_input.loom \
#     regulons.csv \
#     -o auc_output.loom \
#     --num_workers 16

# ============================================================================
# Part 2: Import AUCell Scores into R/Seurat
# ============================================================================

library(Seurat)
library(SCopeLoomR)
library(dplyr)
library(ggplot2)

# --- 5. Read AUCell matrix from pySCENIC output loom ---
loom_file <- "auc_output.loom"  # [USER: set path]
loom_conn <- open_loom(loom_file, mode = "r")

# Extract regulon AUC matrix (cells x regulons)
regulon_auc <- get_regulons_AUC(loom_conn)
close_loom(loom_conn)

# Align cell barcodes between regulon_auc and Seurat object
common_cells <- intersect(colnames(regulon_auc), Cells(seurat_obj))
regulon_auc <- regulon_auc[, common_cells]

# Add AUC scores to Seurat metadata
for (reg in rownames(regulon_auc)) {
  seurat_obj[[reg]] <- regulon_auc[reg, Cells(seurat_obj)]
}

# --- 6. Visualize regulon activity on UMAP ---
FeaturePlot(seurat_obj,
            features = c("FOS(+)", "JUN(+)"),
            cols = c("lightgrey", "red")) +
  ggtitle("AP-1 Regulon Activity (AUCell)")

# ============================================================================
# Part 3: AUCell Scoring with msigdbr Gene Sets
# ============================================================================

library(AUCell)
library(msigdbr)
library(rstatix)

# --- 7. Retrieve gene set from MSigDB ---
chronic_genes <- msigdbr(species = "Homo sapiens", category = "C5") %>%
  dplyr::filter(gs_name == "GOBP_CHRONIC_INFLAMMATORY_RESPONSE") %>%
  dplyr::pull(gene_symbol)

# Prepare gene set list for AUCell
geneSets <- list(Chronic_Inflammation = chronic_genes)

# --- 8. Build AUCell rankings and calculate scores ---
expr_matrix <- GetAssayData(seurat_obj, slot = "data")
cells_rankings <- AUCell_buildRankings(expr_matrix, plotStats = FALSE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)

# Add scores to Seurat metadata
seurat_obj$Chronic_Inflammation_AUC <- getAUC(cells_AUC)["Chronic_Inflammation", ]

# --- 9. Compare scores across cell types and regions ---
plot_data <- data.frame(
  Score    = seurat_obj$Chronic_Inflammation_AUC,
  CellType = seurat_obj$level1_identity,
  Group    = seurat_obj$Type
)

# Filter to immune cells
immune_cells <- c("Mast", "T cells", "Myeloid")
plot_data_clean <- plot_data %>%
  filter(CellType %in% immune_cells) %>%
  filter(Group %in% c("Vertex", "Occipital"))

# Statistical testing per cell type
stat_res <- plot_data_clean %>%
  group_by(CellType) %>%
  wilcox_test(Score ~ Group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

print(stat_res)

# ============================================================================
# Part 4: GeneNMF Meta-program Analysis
# ============================================================================

library(GeneNMF)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# --- 10. Run GeneNMF on immune cell subset ---
# immune_subset: Seurat object containing immune cells only
geneNMF_results <- geneNMF.run(immune_subset,
                               assay = "SCT",
                               n.programs = 9,
                               niter = 50)

# --- 11. Extract meta-program genes ---
mp_genes <- geneNMF.metaprograms(geneNMF_results,
                                  n.genes = 50,
                                  min.confidence = 0.8)

# --- 12. Heatmap of meta-program activity ---
mp_matrix <- geneNMF.score(immune_subset, mp_genes)

ph <- Heatmap(
  mp_matrix,
  name = "MP Score",
  col = colorRamp2(c(0, 0.5, 1), c("lightyellow", "#E64B35", "black")),
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE
)

pdf("Figure_GeneNMF_Clustered_Heatmap.pdf", width = 8, height = 7)
draw(ph)
dev.off()

# ============================================================================
# Part 5: Regulon Heatmap (pheatmap)
# ============================================================================

library(pheatmap)

# --- 13. Average expression heatmap (Occipital vs Vertex) ---
# Example: Wnt pathway genes in DP cells
pathways_map <- list(
  "Wnt_Canonical"    = c("WNT5A", "WNT10B", "CTNNB1", "LEF1", "TCF7"),
  "Wnt_NonCanonical" = c("WNT5A", "ROR2", "RHOA", "RAC1", "JNK1")
)

# Build average expression matrix
avg_object <- AverageExpression(AGA,
                                 group.by = c("level3_identity", "Type"),
                                 assays = "SCT",
                                 return.seurat = TRUE)

pathway.matrix.list <- list()
for (pathway.x in names(pathways_map)) {
  genes <- pathways_map[[pathway.x]]
  genes_found <- genes[genes %in% rownames(avg_object)]
  if (length(genes_found) > 0) {
    sub_mat <- GetAssayData(avg_object, slot = "data")[genes_found, , drop = FALSE]
    pathway.matrix.list[[pathway.x]] <- as.matrix(sub_mat)
  }
}

pathway.matrix <- do.call(rbind, pathway.matrix.list)

pheatmap(pathway.matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale = "row",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         border_color = "white",
         main = "Wnt Signaling (Z-scaled)")
