# ============================================================================
# Differential Expression Analysis and Pathway Enrichment
# Project: AGA single-cell RNA-seq (Nature Communications)
# Input:  Annotated Seurat object with level3_identity and Type (V/O)
# Output: DEG tables, volcano plots, GO/KEGG enrichment results
# ============================================================================

set.seed(42)

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(openxlsx)

# ============================================================================
# Part 1: Differential Expression (HFSC2 example: Vertex vs Occipital)
# ============================================================================

# --- 1. Subset HFSC2 cells ---
obj <- subset(AGA, subset = level3_identity == "HFSC2")
Idents(obj) <- obj$Type

# --- 2. Run FindMarkers (Vertex vs Occipital) ---
deg_results <- FindMarkers(obj,
                           ident.1 = "Vertex",
                           ident.2 = "Occipital",
                           min.pct = 0.1,
                           logfc.threshold = 0)

deg_results$gene <- rownames(deg_results)

# --- 3. Volcano plot ---
volcano_plot <- function(Diff_genes) {
  top5pos <- Diff_genes %>% arrange(desc(avg_log2FC)) %>% head(5)
  top5neg <- Diff_genes %>% arrange(avg_log2FC) %>% head(5)

  ggplot(Diff_genes, aes(x = pct.2 - pct.1, y = avg_log2FC)) +
    geom_point(color = "grey80") +
    geom_hline(yintercept = c(-0.25, 0.25), lty = "dashed",
               size = 1, color = "grey50") +
    geom_text_repel(data = top5pos, aes(label = gene),
                    color = "#CC5E50", size = 4) +
    geom_text_repel(data = top5neg, aes(label = gene),
                    color = "#394C7F", size = 4) +
    geom_point(data = top5pos, aes(color = "up"), size = 2) +
    geom_point(data = top5neg, aes(color = "down"), size = 2) +
    scale_color_manual(values = c("up" = "#CC5E50", "down" = "#394C7F")) +
    theme_bw(base_size = 14) +
    theme(panel.grid = element_blank()) +
    xlab(expression(Delta ~ "Percentage Difference")) +
    ylab("Log2-Fold Change")
}

p <- volcano_plot(deg_results)
print(p)

# ============================================================================
# Part 2: GO Biological Process Enrichment
# ============================================================================

# --- 4. Prepare gene lists ---
sig_genes <- deg_results %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) >= 0.25)

# Vertex-upregulated genes (putative disease-associated)
genes_up_symbol <- sig_genes %>% filter(avg_log2FC > 0) %>% pull(gene)
# Occipital-upregulated genes (putative protective)
genes_down_symbol <- sig_genes %>% filter(avg_log2FC < 0) %>% pull(gene)

# Convert gene symbols to Entrez IDs
symbol_to_entrez <- bitr(c(genes_up_symbol, genes_down_symbol),
                         fromType = "SYMBOL",
                         toType   = "ENTREZID",
                         OrgDb    = org.Hs.eg.db)

genes_up_entrez <- symbol_to_entrez %>%
  filter(SYMBOL %in% genes_up_symbol) %>% pull(ENTREZID)
genes_down_entrez <- symbol_to_entrez %>%
  filter(SYMBOL %in% genes_down_symbol) %>% pull(ENTREZID)

# --- 5. Define enrichment helper function ---
run_go_analysis <- function(entrez_ids, group_name) {
  res <- enrichGO(gene      = entrez_ids,
                  OrgDb     = org.Hs.eg.db,
                  ont       = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  readable  = TRUE)
  if (is.null(res) || nrow(res@result) == 0) return(NULL)

  final_df <- res@result %>%
    dplyr::select(ID, Description, p.adjust, Count, geneID) %>%
    mutate(Group = group_name)
  return(final_df)
}

# --- 6. Run enrichment for both directions ---
go_up   <- run_go_analysis(genes_up_entrez,   "Vertex_Upregulated")
go_down <- run_go_analysis(genes_down_entrez, "Occipital_Upregulated")

# ============================================================================
# Part 3: KEGG Pathway Enrichment
# ============================================================================

# --- 7. KEGG enrichment ---
run_kegg <- function(entrez_ids, group_name) {
  res <- enrichKEGG(gene     = entrez_ids,
                    organism = "hsa",
                    pvalueCutoff = 0.05)
  if (is.null(res) || nrow(res@result) == 0) return(NULL)

  final_df <- res@result %>%
    dplyr::select(ID, Description, p.adjust, Count, geneID) %>%
    mutate(Group = group_name)
  return(final_df)
}

kegg_up   <- run_kegg(genes_up_entrez,   "Vertex_Upregulated")
kegg_down <- run_kegg(genes_down_entrez, "Occipital_Upregulated")

# --- 8. Export results ---
output_list <- list(
  "GO_Vertex_Up"   = go_up,
  "GO_Occipital_Up" = go_down,
  "KEGG_Vertex_Up"  = kegg_up,
  "KEGG_Occipital_Up" = kegg_down
)
output_list <- output_list[!sapply(output_list, is.null)]

if (length(output_list) > 0) {
  write.xlsx(output_list, file = "HFSC2_GO_KEGG_Enrichment.xlsx")
  print("Results saved as: HFSC2_GO_KEGG_Enrichment.xlsx")
}
