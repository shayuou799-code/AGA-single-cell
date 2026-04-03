# ============================================================================
# Pseudotime Trajectory (Monocle3) and Developmental Potential (CytoTRACE2)
# Project: AGA single-cell RNA-seq (Nature Communications)
# Input:  Annotated Seurat object (AGA) with level3_identity
# Output: Pseudotime trajectories, CytoTRACE2 stemness scores
# ============================================================================

set.seed(42)

library(Seurat)
library(monocle3)       # v1.2.9
library(SeuratWrappers)
library(dplyr)
library(ggplot2)

# ============================================================================
# Part 1: Monocle3 Trajectory Analysis
# ============================================================================

# --- 1. Subset cells of interest (e.g., BARX2+ ORS/HFSC lineage) ---
obj <- subset(AGA, subset = level3_identity %in%
                c("HFSC1", "HFSC2", "ORS Basal", "ORS Suprabasal",
                  "Upper HF", "IFE Basal"))

# --- 2. Convert Seurat object to Monocle3 cell_data_set ---
cds <- as.cell_data_set(obj)

# Transfer UMAP from Seurat to Monocle3
cds <- cluster_cells(cds, reduction_method = "UMAP")

# --- 3. Learn principal graph ---
cds <- learn_graph(cds, use_partition = FALSE)

# --- 4. Order cells in pseudotime ---
# Set root node at HFSC1 population (most quiescent/stem-like)
# get_earliest_principal_node() selects the node closest to the specified cluster
root_cells <- colnames(obj)[obj$level3_identity == "HFSC1"]
cds <- order_cells(cds, root_cells = root_cells)

# --- 5. Visualize pseudotime trajectory ---
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 3) +
  ggtitle("Monocle3 Pseudotime: HFSC → ORS Lineage")

# --- 6. Compare pseudotime between Vertex and Occipital ---
pseudotime_df <- data.frame(
  pseudotime = pseudotime(cds),
  region = obj$Type,
  celltype = obj$level3_identity
)

ggplot(pseudotime_df, aes(x = pseudotime, fill = region)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~celltype, scales = "free_y") +
  scale_fill_manual(values = c("Occipital" = "#4DBBD5", "Vertex" = "#E64B35")) +
  theme_classic() +
  ggtitle("Pseudotime Distribution: Vertex vs Occipital")

# ============================================================================
# Part 2: CytoTRACE2 Developmental Potential
# ============================================================================

library(CytoTRACE2)  # v1.0.0

# --- 7. Run CytoTRACE2 on the subset ---
# CytoTRACE2 infers developmental potential from scRNA-seq data
# Higher score = higher stemness/developmental potential
cytotrace2_result <- cytotrace2(obj,
                                 species = "human",
                                 slot_type = "counts")

# --- 8. Add CytoTRACE2 scores to Seurat metadata ---
obj$CytoTRACE2_Score <- cytotrace2_result$CytoTRACE2_Score
obj$CytoTRACE2_Potency <- cytotrace2_result$CytoTRACE2_Potency

# --- 9. Visualize CytoTRACE2 scores on UMAP ---
FeaturePlot(obj, features = "CytoTRACE2_Score",
            reduction = "umap") +
  scale_color_viridis_c() +
  ggtitle("CytoTRACE2 Developmental Potential")

# --- 10. Compare stemness scores across cell types and regions ---
VlnPlot(obj, features = "CytoTRACE2_Score",
        group.by = "level3_identity",
        split.by = "Type",
        pt.size = 0) +
  ggtitle("CytoTRACE2 Score: Vertex vs Occipital")

# --- 11. Statistical comparison of CytoTRACE2 scores ---
# HFSC1 should show highest stemness; HFSC2 intermediate
cytotrace_stats <- obj@meta.data %>%
  group_by(level3_identity, Type) %>%
  summarise(
    mean_score = mean(CytoTRACE2_Score, na.rm = TRUE),
    sd_score   = sd(CytoTRACE2_Score, na.rm = TRUE),
    n          = n(),
    .groups = "drop"
  )
print(cytotrace_stats)

# sessionInfo()  # Uncomment to print session information
