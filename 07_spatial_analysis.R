# ============================================================================
# Spatial Transcriptomics: CosMx Data Loading, Seurat Label Transfer,
# and Giotto Spatial Visualization
# Project: AGA single-cell RNA-seq (Nature Communications)
# Input:  NanoString CosMx data, scRNA-seq reference (AGA Seurat object)
# Output: Spatially annotated cells, in situ plots
# ============================================================================

set.seed(42)

library(Seurat)
library(Giotto)
library(ggplot2)
library(dplyr)
library(data.table)

# ============================================================================
# Part 1: Load CosMx / NanoString Spatial Data into Seurat
# ============================================================================

# --- 1. Load NanoString CosMx data ---
# [USER: set path] to NanoString output directory
data_dir <- "data/cosmx"  # [USER: set path]
sce <- LoadNanostring(data.dir = data_dir, fov = "hair.follicle")

# --- 2. QC: compute mitochondrial percentage ---
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
VlnPlot(sce,
        features = c("nFeature_Nanostring", "nCount_Nanostring", "percent.mt"),
        ncol = 3, pt.size = 0)

# --- 3. Filter low-quality cells ---
sce <- subset(sce,
              subset = nCount_Nanostring > 30 &
                       nFeature_Nanostring > 20)

# --- 4. Normalize with SCTransform ---
sce <- SCTransform(sce, assay = "Nanostring",
                   clip.range = c(-10, 10), verbose = FALSE)

# --- 5. Dimensionality reduction and clustering ---
sce <- RunPCA(sce, assay = "SCT", npcs = 50, verbose = FALSE)
sce <- RunUMAP(sce, reduction = "pca", dims = 1:30,
               n.neighbors = 30, min.dist = 0.01)
sce <- FindNeighbors(sce, reduction = "pca", dims = 1:30)
sce <- FindClusters(sce, resolution = 0.5)

DimPlot(sce, reduction = "umap", label = TRUE)

# ============================================================================
# Part 2: Label Transfer from scRNA-seq Reference
# ============================================================================

# --- 6. Load scRNA-seq reference object ---
# sc_reference: annotated scRNA-seq Seurat object with cell_type metadata
# [USER: set path]
sc_reference <- readRDS("data/AGA_curated.rds")

# Verify reference annotations
DimPlot(sc_reference, group.by = "level3_identity", label = TRUE) + NoLegend()

# --- 7. Find transfer anchors ---
anchors <- FindTransferAnchors(
  reference = sc_reference,
  query = sce,
  normalization.method = "SCT",
  dims = 1:30
)

# --- 8. Transfer labels ---
predictions <- TransferData(
  anchorset = anchors,
  refdata = sc_reference$level3_identity,
  dims = 1:30
)

sce <- AddMetaData(sce, metadata = predictions)

# --- 9. Visualize transferred labels ---
DimPlot(sce, group.by = "predicted.id", label = TRUE, label.size = 3) +
  ggtitle("Transferred Cell Types (UMAP)")

SpatialDimPlot(sce, group.by = "predicted.id",
               label = TRUE, label.size = 3) +
  ggtitle("Transferred Cell Types (Spatial)")

# ============================================================================
# Part 3: Create Giotto Object from Seurat
# ============================================================================

# --- 10. Extract data for Giotto ---
raw_exprs <- GetAssayData(sce, assay = "Nanostring", slot = "counts")

# Spatial coordinates
spatial_locs <- GetTissueCoordinates(sce)
spatial_locs_dt <- data.table(
  cell_ID = rownames(spatial_locs),
  sdimx   = spatial_locs[, 1],
  sdimy   = spatial_locs[, 2]
)

# Cell metadata
cell_metadata <- data.table(
  cell_ID          = colnames(sce),
  level3_identity  = sce$predicted.id,
  level1_identity  = sce$predicted.id  # simplified; map as needed
)

# --- 11. Create Giotto object ---
g <- createGiottoObject(
  raw_exprs     = raw_exprs,
  spatial_locs  = spatial_locs_dt,
  cell_metadata = cell_metadata
)

# ============================================================================
# Part 4: Giotto Spatial Visualization
# ============================================================================

# --- 12. Spatial + UMAP side-by-side plot ---
spatDimPlot2D(g,
  plot_alignment    = "horizontal",
  dim_reduction_to_use = "umap",
  cell_color        = "level1_identity",
  show_image        = FALSE,
  dim_point_size    = 1.0,
  dim_point_shape   = "no_border",
  spat_point_size   = 0.3,
  spat_point_shape  = "no_border",
  dim_show_cluster_center = TRUE,
  dim_label_size    = 3,
  show_legend       = TRUE,
  title = "Cell Type Distribution in Molecular and Spatial Space"
)

# --- 13. In situ plot highlighting specific cell types ---
spatInSituPlotPoints(g,
  feats = list("rna" = c("AREG", "EGFR")),
  feats_color_code = c("AREG" = "red", "EGFR" = "blue"),
  point_size = 0.6,
  show_polygon = TRUE,
  polygon_feat_type = "cell",
  polygon_color = "black",
  polygon_fill = "level3_identity",
  polygon_fill_as_factor = TRUE,
  polygon_fill_code = c("Mast" = "#55efc4", "HFSC1" = "#0984e3",
                         "Others" = "grey"),
  polygon_alpha = 0.8,
  polygon_line_size = 0.05,
  background_color = "black",
  show_image = FALSE,
  title = "Mast-HFSC Spatial Proximity"
)

# --- 14. Highlight individual cell types in spatial context ---
all_level3_clusters <- unique(pDataDT(g)$level3_identity)

for (cluster_name in all_level3_clusters) {
  safe_name <- gsub("[^a-zA-Z0-9]", "_", cluster_name)

  spatPlot2D(g,
    cell_color           = "level3_identity",
    select_cell_groups   = cluster_name,
    show_other_cells     = TRUE,
    other_cell_color     = "grey",
    other_cell_alpha     = 0.2,
    spat_point_size      = 0.5,
    title = paste("Spatial Location of", cluster_name),
    save_plot = TRUE,
    save_param = list(
      save_name = paste0("spatPlot_", safe_name, ".png"),
      base_width = 12, base_height = 10, dpi = 300
    ),
    show_plot = FALSE
  )
}

# --- 15. Spatial gene expression ---
spatFeatPlot2D(g,
  feats = c("KRT15", "CD200", "WNT5A"),
  expression_values = "normalized",
  point_size = 0.8,
  save_plot = TRUE,
  save_param = list(save_name = "spatial_gene_expression.png")
)
