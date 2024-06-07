library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(SCINA)
library(future)
library(ggplot2)
library(scCustomize)
options(future.globals.maxSize = 256000 * 1024^2)
plan("multisession", workers = 10)
theme_set(theme_bw(base_size = 14))

#### ECMO-ARR.#####
dat_ECMOARR <- Read10X_h5("../data/TD006607-ECMO-ARR/outs/filtered_feature_bc_matrix.h5")
#### VSARR. #####
dat_VSARR <- Read10X_h5("../data/TD006607-VSARR/outs/filtered_feature_bc_matrix.h5")
#### VSSHAM. ####
dat_VSSHAM <- Read10X_h5("../data/TD006607-VSSHAM/outs/filtered_feature_bc_matrix.h5")

### Creating Seurat objects. ####
C_ECMOARR <- CreateSeuratObject(counts = dat_ECMOARR, project = "ECMOARR", min.cells = 3, min.features = 200)
C_ECMOARR <- PercentageFeatureSet(C_ECMOARR, "^MT-", col.name = "percent_mito")
C_ECMOARR <- PercentageFeatureSet(C_ECMOARR, "^RP[SL]", col.name = "percent_ribo")
VlnPlot(C_ECMOARR, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), pt.size = 0.1, ncol = 3) + NoLegend()
C_ECMOARR <- subset(C_ECMOARR, subset = nFeature_RNA < 5000 & nCount_RNA < 20000)

C_VSARR <- CreateSeuratObject(counts = dat_VSARR, project = "VSARR", min.cells = 3, min.features = 200)
C_VSARR <- PercentageFeatureSet(C_VSARR, "^MT-", col.name = "percent_mito")
C_VSARR <- PercentageFeatureSet(C_VSARR, "^RP[SL]", col.name = "percent_ribo")
VlnPlot(C_VSARR, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), pt.size = 0.1, ncol = 3) + NoLegend()
C_VSARR <- subset(C_VSARR, subset = nFeature_RNA < 5000 & nCount_RNA < 20000)

C_VSSHAM <- CreateSeuratObject(counts = dat_VSSHAM, project = "VSSHAM", min.cells = 3, min.features = 200)
C_VSSHAM <- PercentageFeatureSet(C_VSSHAM, "^MT-", col.name = "percent_mito")
C_VSSHAM <- PercentageFeatureSet(C_VSSHAM, "^RP[SL]", col.name = "percent_ribo")
VlnPlot(C_VSSHAM, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), pt.size = 0.1, ncol = 3) + NoLegend()
C_VSSHAM <- subset(C_VSSHAM, subset = nFeature_RNA < 5000 & nCount_RNA < 20000)

alldata <- merge(C_ECMOARR, c(C_VSARR, C_VSSHAM), add.cell.ids = c("ECMOARR", "VSARR", "VSSHAM"))
rm(dat_ECMOARR, dat_VSARR, dat_VSSHAM, C_ECMOARR, C_VSSHAM, C_VSARR)
class(alldata[["RNA"]])
alldata <- NormalizeData(alldata)
alldata <- FindVariableFeatures(alldata)
alldata <- ScaleData(alldata)
alldata <- RunPCA(alldata)
alldata <- FindNeighbors(alldata, dims = 1:30, reduction = "pca")
alldata <- FindClusters(alldata, resolution = 0.5, cluster.name = "unintegrated_clusters")
alldata <- RunUMAP(alldata, reduction = "pca", dims = 1:30, reduction.name = "umap.pca")
#p1 <- DimPlot(alldata, reduction = "umap.pca", group.by = "orig.ident")

alldata <- IntegrateLayers(
  object = alldata, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

alldata <- IntegrateLayers(
  object = alldata, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

###############################
############## cca #############
###############################
plan("multisession", workers = 1)
alldata <- FindNeighbors(alldata, reduction = "integrated.cca", dims = 1:30)
alldata <- FindClusters(alldata, resolution = 0.1, cluster.name = "cca_clusters")
alldata <- RunUMAP(alldata, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
p2 <- DimPlot(alldata, reduction = "umap.cca", group.by = "cca_clusters")
p2
###############################
############## rpca #############
###############################
alldata <- FindNeighbors(alldata, reduction = "integrated.rpca", dims = 1:30)
alldata <- FindClusters(alldata, resolution = 0.5, cluster.name = "rpca_clusters")
alldata <- RunUMAP(alldata, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
p3 <- DimPlot(alldata, reduction = "umap.rpca", group.by = "rpca_clusters")

saveRDS(alldata, "Brain_scRNA_seurat_v5.rds")

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
alldata <- JoinLayers(alldata)
sweep.res.list_kidney <- paramSweep(alldata, PCs = 1:10, sct = FALSE)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(alldata@meta.data$seurat_clusters)           ## ex: annotations <- alldata@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(alldata@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
alldata <- doubletFinder(alldata, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
alldata <- doubletFinder(alldata, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_8714", sct = FALSE)
data.filt = alldata[, alldata@meta.data$DF.classifications_0.25_0.09_8714 == "Singlet"]
saveRDS(data.filt, "Brain_scRNA_seurat_v5_filter_doublets.rds")

##############################################
#############Find All Markers#################
##############################################
data.filt<- readRDS("Brain_scRNA_seurat_v5.rds")
data.filt <- JoinLayers(data.filt)
Idents(data.filt) <- data.filt@meta.data$cca_clusters
markers <- FindAllMarkers(data.filt)
write.table(markers, file = "brain_markers.txt", sep = "\t", quote = F, row.names = F)
markers <- read.table("brain_markers.txt", header = T)

#Sys.setenv(OPENAI_API_KEY = 'sk-proj-0WWnKAwpapcwVUBbjUSGT3BlbkFJB1RaRW8nEK7E04NQwixg')
#library(GPTCelltype)
#library(openai)
#res <- gptcelltype(markers, 
#                   tissuename = 'pig brain', 
#                   model = 'gpt-4'
#)
#data.filt@meta.data$celltype <- as.factor(res[as.character(Idents(data.filt))])
marker_ref <- c("SLC1A3","SOX2", "GFAP", "HMGB2", "ASCL1", "REC4","MKI67", "SOX4", "NNAT","PROX1","NPY1R",
                "ERC2","CPLX2","RGS4","ELAVL2","GAD2","CSPG4","BCAS1","MOG","CIOA","RELN","FLT1",
                "FOXJ1","LYVET","DCN","AQP4", "C1QA")
marker_ref <- c("C1QA", "GFAP", "ASCL1", "NEUROD2", "GAD1", "PROX1", "MEIS2",
                "LHX6", "NR2F2", "OLIG2", "MBP", "PTPRC", "SPARC")
Stacked_VlnPlot(data.filt, features = marker_ref, group.by="cca_clusters")
VlnPlot(data.filt, features = "MEIS2", group.by = "cca_clusters")
### Cluster 5. ###
FeaturePlot(data.filt, features = "AQP4", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
FeaturePlot(data.filt, features = "GFAP", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)

### Cluster 4. ###
FeaturePlot(data.filt, features = "C1QA", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
FeaturePlot(data.filt, features = "CSF1R", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
FeaturePlot(data.filt, features = "CX3CR1", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
### Cluster 0. ##
FeaturePlot(data.filt, features = "SLC17A6", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
#FeaturePlot(data.filt, features = "GRIN1", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
#FeaturePlot(data.filt, features = "CAMK2A", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
### Cluster 2.##
FeaturePlot(data.filt, features = "MOG", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
FeaturePlot(data.filt, features = "MAG", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)

### Cluster 3.##
FeaturePlot(data.filt, features = "VCAN", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
FeaturePlot(data.filt, features = "CSPG4", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
FeaturePlot(data.filt, features = "PDGFRA", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)

### Cluster 6.##
FeaturePlot(data.filt, features = "LHX6", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
FeaturePlot(data.filt, features = "SOX6", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
FeaturePlot(data.filt, features = "GAD1", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
FeaturePlot(data.filt, features = "GAD2", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
FeaturePlot(data.filt, features = "SLC32A1", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
### Cluster 1.##
FeaturePlot(data.filt, features = "GABRB1", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
FeaturePlot(data.filt, features = "DAB1", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
FeaturePlot(data.filt, features = "GABRB3", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)

### Cluster 11.##
FeaturePlot(data.filt, features = "SOX4", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
FeaturePlot(data.filt, features = "HMGB2", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)

### Cluster 8.##
FeaturePlot(data.filt, features = "RUNX1", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
FeaturePlot(data.filt, features = "DOCK8", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
FeaturePlot(data.filt, features = "ARHGAP15", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
### Cluster 9.##
FeaturePlot(data.filt, features = "MOG", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
FeaturePlot(data.filt, features = "BCAS1", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
### Cluster 10.##
FeaturePlot(data.filt, features = "LHFPL3", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
FeaturePlot(data.filt, features = "PDGFRA", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
FeaturePlot(data.filt, features = "RUNX1", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
### Cluster 11.##
FeaturePlot(data.filt, features = "APOA1", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
FeaturePlot(data.filt, features = "CLDN5", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)
FeaturePlot(data.filt, features = "FLT1", reduction = "umap.cca", cols = c( "gray","white","#850101"), pt.size = 1)

cell_types <- c(`0` = "ExN",`1` = "InN", `2`= "OL", `3` = "OPC", `4` = "Microglia", 
                `5` = "Astrocytes", `6` = "MGE InN", `7` = "Astrocytes", `8` = "Blood cells", `9` = "NFOL", 
                `10` = "Progenitor", `11` = "Endothelial")
########## MGE InN: medial ganglionic eminence inhibitory neuron (LHX6)
########## CGE InN: caudal ganglionic eminence inhibitory neuron ()
########## LGE InN: lateral ganglionic erminence inhibitory neuron ()
########## ExN: excitatory neurons
########## OPC: Oligodendrocyte Progenitor Cells ()
########## OL: Oligodendrocytes ()
########## IPC: intermediate progenitor cells
########## NFOL: newly formed oligodendrocytes
data.filt <- RenameIdents(data.filt, cell_types)
seurat_clusters <- data.filt@meta.data$seurat_clusters
seurat_clusters <- as.numeric(as.character(seurat_clusters))
seurat_clusters[seurat_clusters == 0] <- "ExN"
seurat_clusters[seurat_clusters == 1] <- "InN"
seurat_clusters[seurat_clusters == 2] <- "OL"
seurat_clusters[seurat_clusters == 3] <- "OPC"
seurat_clusters[seurat_clusters == 4] <- "Microglia"
seurat_clusters[seurat_clusters == 5] <- "Astrocytes"
seurat_clusters[seurat_clusters == 6] <- "MGE InN"
seurat_clusters[seurat_clusters == 7] <- "Astrocytes"
seurat_clusters[seurat_clusters == 8] <- "Blood cells"
seurat_clusters[seurat_clusters == 9] <- "NFOL"
seurat_clusters[seurat_clusters == 10] <- "Progenitor"
seurat_clusters[seurat_clusters == 11] <- "Endothelial"
data.filt$cell_type <- seurat_clusters
saveRDS(data.filt, "Brain_scRNA_seurat_v5_filter_doublets_cell_type.rds")
DimPlot(data.filt, group.by = "cell_type", reduction = "umap.cca", label = T)
############## Differential analysis.#######
############## VSSHAM vs. VSARR #########
for (n in unique(data.filt@meta.data[, 'cell_type'])){
  print (n)
    cell_selection <- subset(data.filt, subset = (orig.ident == "VSSHAM" | orig.ident ==  "VSARR"), cells = colnames(data.filt)[data.filt@meta.data[, 'cell_type'] == n])
    if(length(table(cell_selection@meta.data$orig.ident)) == 2 && table(cell_selection@meta.data$orig.ident)[[1]] > 10 && table(cell_selection@meta.data$orig.ident)[[2]] > 10){
      cell_selection <- SetIdent(cell_selection, value = "orig.ident")
      DGE_cell_selection <- FindMarkers(cell_selection, log2FC.threshold = 0.1, ident.1 = "VSARR", ident.2 = "VSSHAM")
      n <- gsub(" ",".",n)
      x <- c(n, ".VSSHAM-vs-VSARR.gene.txt"); out <- paste(x, collapse=""); write.table(DGE_cell_selection, file = out, sep = "\t")
    }
}
############## VSSHAM vs. ECMOARR #########
for (n in unique(data.filt@meta.data[, 'cell_type'])){
  print (n)
  cell_selection <- subset(data.filt, subset = (orig.ident == "VSSHAM" | orig.ident ==  "ECMOARR"), cells = colnames(data.filt)[data.filt@meta.data[, 'cell_type'] == n])
  if(length(table(cell_selection@meta.data$orig.ident)) == 2 && table(cell_selection@meta.data$orig.ident)[[1]] > 10 && table(cell_selection@meta.data$orig.ident)[[2]] > 10){
    cell_selection <- SetIdent(cell_selection, value = "orig.ident")
    DGE_cell_selection <- FindMarkers(cell_selection, log2FC.threshold = 0.1, ident.1 = "ECMOARR", ident.2 = "VSSHAM")
    n <- gsub(" ",".",n)
    x <- c(n, ".VSSHAM-vs-ECMOARR.gene.txt"); out <- paste(x, collapse=""); write.table(DGE_cell_selection, file = out, sep = "\t")
  }
}
############## ECMOARR vs. VSARR #########
for (n in unique(data.filt@meta.data[, 'cell_type'])){
  print (n)
  cell_selection <- subset(data.filt, subset = (orig.ident == "ECMOARR" | orig.ident ==  "VSARR"), cells = colnames(data.filt)[data.filt@meta.data[, 'cell_type'] == n])
  if(length(table(cell_selection@meta.data$orig.ident)) == 2 && table(cell_selection@meta.data$orig.ident)[[1]] > 10 && table(cell_selection@meta.data$orig.ident)[[2]] > 10){
    cell_selection <- SetIdent(cell_selection, value = "orig.ident")
    DGE_cell_selection <- FindMarkers(cell_selection, log2FC.threshold = 0.1, ident.1 = "ECMOARR", ident.2 = "VSARR")
    n <- gsub(" ",".",n)
    x <- c(n, ".VSARR-vs-ECMOARR.gene.txt"); out <- paste(x, collapse=""); write.table(DGE_cell_selection, file = out, sep = "\t")
  }
}
