# ============================================================================
# Cell-Cell Communication Analysis with CellChat
# Project: AGA single-cell RNA-seq (Nature Communications)
# Input:  Annotated Seurat object with level3_identity and Type (V/O)
# Output: CellChat objects, differential interaction analysis,
#         chord diagrams, heatmaps (Occipital vs Vertex)
# ============================================================================

set.seed(42)

library(CellChat)
library(patchwork)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# ============================================================================
# Part 1: Create CellChat Objects (one per region)
# ============================================================================

# --- 1. Split Seurat object by region ---
AGA_occ <- subset(AGA, subset = Type == "Occipital")
AGA_ver <- subset(AGA, subset = Type == "Vertex")

# --- 2. Create CellChat object for Occipital ---
cellchat_occ <- createCellChat(object = AGA_occ, group.by = "level3_identity")
CellChatDB.use <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
cellchat_occ@DB <- CellChatDB.use

# Preprocessing, communication inference
cellchat_occ <- subsetData(cellchat_occ)
cellchat_occ <- identifyOverExpressedGenes(cellchat_occ)
cellchat_occ <- identifyOverExpressedInteractions(cellchat_occ)
cellchat_occ <- computeCommunProb(cellchat_occ, type = "triMean")
cellchat_occ <- filterCommunication(cellchat_occ, min.cells = 5)
cellchat_occ <- computeCommunProbPathway(cellchat_occ)
cellchat_occ <- aggregateNet(cellchat_occ)
cellchat_occ <- netAnalysis_computeCentrality(cellchat_occ)

# --- 3. Repeat for Vertex ---
cellchat_ver <- createCellChat(object = AGA_ver, group.by = "level3_identity")
cellchat_ver@DB <- CellChatDB.use

cellchat_ver <- subsetData(cellchat_ver)
cellchat_ver <- identifyOverExpressedGenes(cellchat_ver)
cellchat_ver <- identifyOverExpressedInteractions(cellchat_ver)
cellchat_ver <- computeCommunProb(cellchat_ver, type = "triMean")
cellchat_ver <- filterCommunication(cellchat_ver, min.cells = 5)
cellchat_ver <- computeCommunProbPathway(cellchat_ver)
cellchat_ver <- aggregateNet(cellchat_ver)
cellchat_ver <- netAnalysis_computeCentrality(cellchat_ver)

# ============================================================================
# Part 2: Merge and Compare
# ============================================================================

# --- 4. Merge CellChat objects for comparison ---
object.list <- list(Occipital = cellchat_occ, Vertex = cellchat_ver)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# ============================================================================
# Part 3: Differential Interaction Visualization
# ============================================================================

# --- 5. Chord diagram: differential interactions ---
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = TRUE,
                          title.name = "Differential Interactions (Count)")
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight",
                          title.name = "Differential Interactions (Strength)")

# --- 6. Stacked bar: interaction count ranking ---
rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE)

# --- 7. Bubble plot: specific cell pair interactions ---
# Example: Mast -> HFSC signaling comparison
netVisual_bubble(cellchat,
                 sources.use = "Mast",
                 targets.use = c("HFSC1", "HFSC2"),
                 comparison = c(1, 2),
                 angle.x = 45,
                 remove.isolate = TRUE,
                 title.name = "Mast -> HFSC Signaling (Occ vs Ver)")

# --- 8. Chord diagram for specific pathway ---
pathways.show <- "SPP1"
for (i in seq_along(object.list)) {
  netVisual_aggregate(object.list[[i]],
                      signaling = pathways.show,
                      layout = "chord",
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

# --- 9. Validate ligand-receptor expression ---
plotGeneExpression(cellchat,
                  features = c("SPP1", "CD44"),
                  split.by = "datasets",
                  colors.ggplot = TRUE)

# --- 10. Incoming signal to Mast cells ---
par(mfrow = c(1, 1))
netVisual_aggregate(cellchat,
                    targets.use = "Mast",
                    layout = "chord",
                    remove.isolate = TRUE,
                    title.name = "Incoming Signals to Mast")

# ============================================================================
# Part 4: Signaling Role Heatmap
# ============================================================================

# --- 11. Signaling role scoring ---
pathway.union <- union(object.list[[1]]@netP$pathways,
                       object.list[[2]]@netP$pathways)

ht1 <- netAnalysis_signalingRole_heatmap(object.list[[1]],
                                          pattern = "outgoing",
                                          signaling = pathway.union,
                                          title = "Occipital: Outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(object.list[[2]],
                                          pattern = "outgoing",
                                          signaling = pathway.union,
                                          title = "Vertex: Outgoing")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

# --- 12. Save results ---
# pdf("CellChat_DiffInteraction.pdf", width = 10, height = 5)
# par(mfrow = c(1, 2), xpd = TRUE)
# netVisual_diffInteraction(cellchat, weight.scale = TRUE, ...)
# dev.off()
