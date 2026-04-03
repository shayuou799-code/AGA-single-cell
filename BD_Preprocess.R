library(readr)
library(tidyr)
library(plyr)
library(dplyr)
library(DoubletFinder)
library(Seurat)
library(celda)


matrix_formation = function(filename, SampleID, Unfiltered = FALSE){
  dat = read_tsv(filename, comment = "#")
  dat = dat[, c("Cell_Index", "Gene", "RSEC_Adjusted_Molecules")]
  
  if(Unfiltered){
    dat_N = dat %>% group_by(Cell_Index) %>% summarise(Gene_Num = n())
    retained_cells = dat_N$Cell_Index[dat_N$Gene_Num >= 500 & dat_N$Gene_Num <= 5000]
    dat = dat[which(dat$Cell_Index %in% retained_cells), ]
  }
  
  dat_matrix = dat %>% tidyr::pivot_wider(names_from = Cell_Index, values_from = RSEC_Adjusted_Molecules)
  genename = as.vector(t(dat_matrix[, 1]))
  
  dat_matrix = dat_matrix[, -1]
  dat_matrix = data.table::setnafill(dat_matrix, fill = 0)
  dat_matrix = as.matrix(dat_matrix)
  rownames(dat_matrix) = genename
  colnames(dat_matrix) = paste(SampleID, colnames(dat_matrix), sep ="_")
  return(dat_matrix)
}



scObject_create = function(matrix, SampleID, mtPercent = 20){
  #remove ambient RNAs using decontX from celda package. 
  matrix = SingleCellExperiment(list(counts = matrix))
  matrix = decontX(matrix)
  matrix = decontXcounts(matrix)
  #create the Seurat object. 
  scObject = CreateSeuratObject(matrix, project = SampleID, min.cells = 3, min.features = 200)
  scObject[["percent.mt"]] = PercentageFeatureSet(scObject, pattern = "^MT-|^mt-")
  scObject = subset(scObject, nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < mtPercent)
  #scObject = SetIdent(scObject, value = "SeuratObject")
  p = VlnPlot(scObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  cowplot::save_plot(p, filename = paste(SampleID, "_QC.jpeg", sep = ""), base_height = 6, base_width = 6, 
                     dpi = 300)
  scObject = scObject %>% NormalizeData() %>% ScaleData() %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    RunPCA() %>% FindNeighbors(reduction = "pca", dims = 1:50) %>%
  FindClusters(resolution = 0.6)
  
  return(scObject)
}

Doublet_removal = function(scObject, percent = 0.075){
  sweep.data = paramSweep(scObject, PCs = 1:10)
  sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
  bcmvn= find.pK(sweep.stats)
  homotypic.prop=modelHomotypic(scObject@meta.data$RNA_snn_res.0.6)         
  nExp_poi=round(percent*length(scObject$orig.ident))  
  nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
  scObject=doubletFinder(scObject, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                      nExp = nExp_poi.adj)
  doubletsID = scObject@meta.data[, grep("DF.classifications", colnames(scObject@meta.data))]
  scObject = subset(scObject, cells = colnames(scObject)[which(doubletsID == "Singlet")])
  scObject@meta.data[, grep("DF.classifications", colnames(scObject@meta.data))] = NULL
  return(scObject)
}

BD_Preprocessor = function(filename, SampleID, ...){
  BD_matrix = matrix_formation(filename, SampleID, ...)
  BD_object = scObject_create(BD_matrix, SampleID)
  BD_object = Doublet_removal(BD_object)
  write_rds(BD_object, path = paste(SampleID, ".rds", sep = ""))
}



