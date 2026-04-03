#!/usr/bin/env python3
# ============================================================================
# RNA Velocity and Pseudotime Analysis (scVelo + Scanpy)
# Project: AGA single-cell RNA-seq (Nature Communications)
# Input:  Merged loom files + Seurat-exported metadata/UMAP
# Output: Velocity stream plots, PAGA, pseudotime, driver gene ranking
# ============================================================================

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
import seaborn as sns
import os

np.random.seed(42)

# Publication-quality plot settings
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.family"] = "Arial"
scv.settings.verbosity = 3

# ============================================================================
# Part 1: Data Loading and Preprocessing
# ============================================================================

# --- 1. Read merged h5ad (loom data + Seurat annotations) ---
# [USER: set path] to the merged anndata file
adata = sc.read_h5ad("hfsc_merged_final.h5ad")

# Ensure categorical metadata
adata.obs["level2_identity"] = adata.obs["level2_identity"].astype("category")
adata.obs["Type"] = adata.obs["Type"].astype("category")

print(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

# --- 2. Filter and normalize ---
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# ============================================================================
# Part 2: Velocity Estimation (Dynamical Model)
# ============================================================================

# --- 3. Recover full transcriptional dynamics ---
scv.tl.recover_dynamics(adata, n_jobs=8)

# --- 4. Compute velocity ---
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)

# --- 5. Velocity stream plot on UMAP ---
scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="level2_identity",
    title="RNA Velocity Stream (Dynamical Model)",
    save="velocity_stream.pdf"
)

# ============================================================================
# Part 3: PAGA Trajectory Analysis
# ============================================================================

# --- 6. Compute PAGA with velocity information ---
scv.tl.paga(adata, groups="level2_identity", use_time_prior="velocity_pseudotime")

# --- 7. Plot PAGA directed graph ---
scv.pl.paga(
    adata,
    basis="umap",
    transitions="transitions_confidence",
    size=50,
    alpha=0.1,
    min_edge_width=2,
    node_size_scale=1.5,
    title="PAGA Directed Differentiation Trajectory",
    save="paga_trajectory.pdf"
)

# ============================================================================
# Part 4: Pseudotime (Diffusion Pseudotime)
# ============================================================================

# --- 8. Set root cell (HFSC1 = stem cell population) ---
root_cells = adata.obs[adata.obs["level2_identity"] == "HFSC1"].index
adata.uns["iroot"] = np.where(adata.obs_names == root_cells[0])[0][0]
print(f"Set HFSC1 as developmental root.")

# --- 9. Compute diffusion pseudotime ---
sc.tl.diffmap(adata)
sc.tl.dpt(adata)

# --- 10. Visualize pseudotime ---
sc.pl.umap(adata, color=["dpt_pseudotime", "level2_identity", "Type"],
           cmap="viridis", save="_pseudotime.pdf")

# ============================================================================
# Part 5: Velocity Driver Gene Ranking
# ============================================================================

# --- 11. Rank velocity driver genes by group ---
scv.tl.rank_velocity_genes(adata, groupby="Type", min_corr=0.3)

# Extract top drivers per group
names_df = pd.DataFrame(adata.uns["rank_velocity_genes"]["names"])
scores_df = pd.DataFrame(adata.uns["rank_velocity_genes"]["scores"])

print("Top 10 velocity driver genes (Vertex):")
print(names_df["Vertex"].head(10).tolist())

print("\nTop 10 velocity driver genes (Occipital):")
print(names_df["Occipital"].head(10).tolist())

# --- 12. Velocity confidence comparison ---
scv.tl.velocity_confidence(adata)
scv.tl.velocity_pseudotime(adata)
scv.tl.latent_time(adata)

# --- 13. Phase portrait for key genes ---
scv.pl.scatter(adata, basis=["FOS", "JUND"],
               color="level2_identity",
               frameon=False,
               save="phase_portrait_FOS_JUND.pdf")

# ============================================================================
# Part 6: Pseudotime Statistical Analysis
# ============================================================================

from scipy import stats

# --- 14. Compare pseudotime distributions between regions ---
vertex_time = adata.obs[adata.obs["Type"] == "Vertex"]["dpt_pseudotime"]
occipital_time = adata.obs[adata.obs["Type"] == "Occipital"]["dpt_pseudotime"]

stat, p_value = stats.ks_2samp(vertex_time, occipital_time)
print(f"\nKS Test (Vertex vs Occipital pseudotime):")
print(f"  Statistic: {stat:.4f}")
print(f"  P-value:   {p_value:.4e}")

# --- 15. Cell proportion comparison ---
counts = pd.crosstab(adata.obs["Type"], adata.obs["level2_identity"],
                     normalize="index") * 100
print("\nCell type proportions by region (%):")
print(counts)

# --- 16. Save processed object ---
adata.write("aga_velocity_processed.h5ad")
print("\nAnalysis complete. Saved to: aga_velocity_processed.h5ad")
