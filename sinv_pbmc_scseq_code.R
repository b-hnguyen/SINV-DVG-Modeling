# Single Cell RNA seq analysis of uninfected and 7dpi SINV TE-infected mouse PBMCs
## Ben Nguyen

library("Seurat")
library("dplyr")
library("patchwork")
library("cowplot")
library("ggplot2")
library("fields")
library("ROCR")
library("KernSmooth")
library("Matrix")
library("parallel")
library("clustree")
library("DoubletFinder")
library("EnhancedVolcano")
library("metap")
library('escape')
library('dittoSeq')
library('DESeq2')
library('SingleCellExperiment')
library('Matrix.utils')
library('magrittr')
library('purrr')
library('pheatmap')
library('topGO')
library('enrichR')
library('presto')
library('fgsea')
library('AnnotationDbi')
library('presto')

#Read in Data and QC
counts.data.ctrl = Read10X(data.dir="/dcl02/leased/giffin/SINV_scSeq/fastq/SCIBAR/run_cellranger/run_count_ETro_control/outs/filtered_feature_bc_matrix")
ctrl = CreateSeuratObject(counts = counts.data.ctrl, project = "Control", min.cells = 3)
DefaultAssay(ctrl) = "RNA"
ctrl[["percent.mt"]] = PercentageFeatureSet(ctrl, pattern = "^mt-") #calculate mitochondrial QC metrics and store into meta data
VlnPlot(ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) #violin plots of QC metrics
plot.percentmt.ctrl = FeatureScatter(ctrl, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot.nfeature.ctrl = FeatureScatter(ctrl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot.percentmt.ctrl + plot.nfeature.ctrl #plotting QC metrics
ctrl = subset(ctrl, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
ctrl = NormalizeData(ctrl, verbose = FALSE)
ctrl = FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)
ctrl <- ScaleData(ctrl, verbose = FALSE)
#Standard clustering workflow
ctrl <- RunPCA(ctrl, npcs = 30, verbose = FALSE)
ctrl <- RunUMAP(ctrl, reduction = "pca", dims = 1:30)
ctrl <- FindNeighbors(ctrl, reduction = "pca", dims = 1:30)
ctrl <- FindClusters(ctrl)
#Doublet removal
## pK Identification (no ground-truth)
sweep.res.list_ctrl <- paramSweep_v3(ctrl, PCs = 1:10, sct = FALSE)
sweep.stats_ctrl <- summarizeSweep(sweep.res.list_ctrl, GT = FALSE)
bcmvn_ctrl <- find.pK(sweep.stats_ctrl)
## Homotypic Doublet Proportion Estimate
annotations.ctrl <- ctrl@meta.data$seurat_clusters
homotypic.prop.ctrl <- modelHomotypic(annotations.ctrl)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi.ctrl <- round(0.075*nrow(ctrl@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj.ctrl <- round(nExp_poi.ctrl*(1-homotypic.prop.ctrl))
## Run DoubletFinder with varying classification stringencies
ctrl <- doubletFinder_v3(ctrl, PCs = 1:10, pN = 0.25, pK = 0.2, nExp = nExp_poi.ctrl, reuse.pANN = FALSE, sct = FALSE) #CHANGE pK HERE TO BE pK WITH HIGHEST CORRESPONDING BCmetric FROM bcmv TABLE -- proceeding with 0.2
ctrl <- doubletFinder_v3(ctrl, PCs = 1:10, pN = 0.25, pK = 0.2, nExp = nExp_poi.adj.ctrl, reuse.pANN = "pANN_0.25_0.2_694", sct = FALSE)
ctrl@meta.data[,"CellTypes_DF"] <- ctrl@meta.data$DF.classifications_0.25_0.2_645
ctrl@meta.data$CellTypes_DF[which(ctrl@meta.data$DF.classifications_0.25_0.2_645 == "Doublet")] <- "Doublet"
ctrl@meta.data$CellTypes_DF[which(ctrl@meta.data$DF.classifications_0.25_0.2_645 == "Singlet")] <- "Singlet"
ctrl.doublets = DimPlot(ctrl, group.by="CellTypes_DF", reduction="umap", pt.size=0.5, order=c("Coll.Duct.TC","Doublet"))
ctrl.with.doubs = ctrl #this will save the original Seurat object containing doublet information and both doublets and singlets
ctrl <- subset(ctrl, subset = CellTypes_DF == "Singlet") #this will subset out only singlets from the data - cluster from this step moving forward

counts.data.inf = Read10X(data.dir="/dcl02/leased/giffin/SINV_scSeq/fastq/SCIBAR/run_cellranger/run_count_ETro_infected/outs/filtered_feature_bc_matrix")
inf = CreateSeuratObject(counts = counts.data.inf, project = "Infected", min.cells = 3)
DefaultAssay(inf) = "RNA"
inf[["percent.mt"]] = PercentageFeatureSet(inf, pattern = "^mt-") #calculate mitochondrial QC metrics and store into meta data
VlnPlot(inf, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) #violin plots of QC metrics
plot.percentmt.inf = FeatureScatter(inf, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot.nfeature.inf = FeatureScatter(inf, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot.percentmt.inf + plot.nfeature.inf #plotting QC metrics
inf = subset(inf, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
inf = NormalizeData(inf, verbose = FALSE)
inf = FindVariableFeatures(inf, selection.method = "vst", nfeatures = 2000)
inf <- ScaleData(inf, verbose = FALSE)
#Standard clustering workflow
inf <- RunPCA(inf, npcs = 30, verbose = FALSE)
inf <- RunUMAP(inf, reduction = "pca", dims = 1:30)
inf <- FindNeighbors(inf, reduction = "pca", dims = 1:30)
inf <- FindClusters(inf)
#Doublet removal
## pK Identification (no ground-truth)
sweep.res.list_inf <- paramSweep_v3(inf, PCs = 1:10, sct = FALSE)
sweep.stats_inf <- summarizeSweep(sweep.res.list_inf, GT = FALSE)
bcmvn_inf <- find.pK(sweep.stats_inf)
## Homotypic Doublet Proportion Estimate
annotations.inf <- inf@meta.data$seurat_clusters
homotypic.prop.inf <- modelHomotypic(annotations.inf)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi.inf <- round(0.075*nrow(inf@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj.inf <- round(nExp_poi.inf*(1-homotypic.prop.inf))
## Run DoubletFinder with varying classification stringencies
inf <- doubletFinder_v3(inf, PCs = 1:10, pN = 0.25, pK = 0.24, nExp = nExp_poi.inf, reuse.pANN = FALSE, sct = FALSE) #CHANGE pK HERE TO BE pK WITH HIGHEST CORRESPONDING BCmetric FROM bcmv TABLE -- proceeding with 0.24
inf <- doubletFinder_v3(inf, PCs = 1:10, pN = 0.25, pK = 0.24, nExp = nExp_poi.adj.inf, reuse.pANN = "pANN_0.25_0.24_459", sct = FALSE)
inf@meta.data[,"CellTypes_DF"] <- inf@meta.data$DF.classifications_0.25_0.24_420
inf@meta.data$CellTypes_DF[which(inf@meta.data$DF.classifications_0.25_0.24_420 == "Doublet")] <- "Doublet"
inf@meta.data$CellTypes_DF[which(inf@meta.data$DF.classifications_0.25_0.24_420 == "Singlet")] <- "Singlet"
inf.doublets = DimPlot(inf, group.by="CellTypes_DF", reduction="umap", pt.size=0.5, order=c("Coll.Duct.TC","Doublet"))
inf.with.doubs = inf #this will save the original Seurat object containing doublet information and both doublets and singlets
inf <- subset(inf, subset = CellTypes_DF == "Singlet") #this will subset out only singlets from the data - cluster from this step moving forward

#Perform Integration
pbmc.list <- lapply(X = list(ctrl, inf), FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = pbmc.list)
immune.anchors = FindIntegrationAnchors(object.list = pbmc.list, anchor.features = features)
immune.combined = IntegrateData(anchorset = immune.anchors)
DefaultAssay(immune.combined) = "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.4)#change the resolution here from 0.1 to 0.8 for the Clustree analysis in the next step - can return to the original rds for clustering once appropriate resolution is chosen
## Proceeding with resolution of 0.4 based on clustree results from 20220222

# saveRDS(immune.combined, file = "/dcl02/leased/giffin/SINV_scSeq/ben/pbmc/20220223_sinv_pbmc.rds") #save rds file to avoid computationally intensive steps above

# Clustree analysis
head(immune.combined[[]])
clustree(immune.combined, prefix = "integrated_snn_res.")

# Visualization
plot.clus.treatment = DimPlot(immune.combined, reduction = "umap", group.by = "orig.ident")
plot.clus = DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
plot.clus.treatment + plot.clus

Idents(immune.combined) = "orig.ident" #switches UMAP labeling to samples instead of clusters
DimPlot(immune.combined, reduction = "umap", label = F, split.by="orig.ident") #plots UMAP with labeled samples
Idents(immune.combined) = "seurat_clusters" #switches back to labeling clusters
DimPlot(immune.combined, reduction="umap", split.by="orig.ident") #shows umaps side by side for each sample

# For performing differential expression after integration, we switch back to the original data
DefaultAssay(immune.combined) <- "RNA"
clus0.markers <- FindConservedMarkers(immune.combined, ident.1 = 0, grouping.var = "orig.ident", verbose = FALSE)
clus1.markers <- FindConservedMarkers(immune.combined, ident.1 = 1, grouping.var = "orig.ident", verbose = FALSE)
clus2.markers <- FindConservedMarkers(immune.combined, ident.1 = 2, grouping.var = "orig.ident", verbose = FALSE)
clus3.markers <- FindConservedMarkers(immune.combined, ident.1 = 3, grouping.var = "orig.ident", verbose = FALSE)
clus4.markers <- FindConservedMarkers(immune.combined, ident.1 = 4, grouping.var = "orig.ident", verbose = FALSE)
clus5.markers <- FindConservedMarkers(immune.combined, ident.1 = 5, grouping.var = "orig.ident", verbose = FALSE)
clus6.markers <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "orig.ident", verbose = FALSE)
clus7.markers <- FindConservedMarkers(immune.combined, ident.1 = 7, grouping.var = "orig.ident", verbose = FALSE)
clus8.markers <- FindConservedMarkers(immune.combined, ident.1 = 8, grouping.var = "orig.ident", verbose = FALSE)
clus9.markers <- FindConservedMarkers(immune.combined, ident.1 = 9, grouping.var = "orig.ident", verbose = FALSE)
clus10.markers <- FindConservedMarkers(immune.combined, ident.1 = 10, grouping.var = "orig.ident", verbose = FALSE)
clus11.markers <- FindConservedMarkers(immune.combined, ident.1 = 11, grouping.var = "orig.ident", verbose = FALSE)
clus12.markers <- FindConservedMarkers(immune.combined, ident.1 = 12, grouping.var = "orig.ident", verbose = FALSE)
clus13.markers <- FindConservedMarkers(immune.combined, ident.1 = 13, grouping.var = "orig.ident", verbose = FALSE)
clus14.markers <- FindConservedMarkers(immune.combined, ident.1 = 14, grouping.var = "orig.ident", verbose = FALSE)
clus15.markers <- FindConservedMarkers(immune.combined, ident.1 = 15, grouping.var = "orig.ident", verbose = FALSE)
clus16.markers <- FindConservedMarkers(immune.combined, ident.1 = 16, grouping.var = "orig.ident", verbose = FALSE)
clus17.markers <- FindConservedMarkers(immune.combined, ident.1 = 17, grouping.var = "orig.ident", verbose = FALSE)
clus18.markers <- FindConservedMarkers(immune.combined, ident.1 = 18, grouping.var = "orig.ident", verbose = FALSE)
clus19.markers <- FindConservedMarkers(immune.combined, ident.1 = 19, grouping.var = "orig.ident", verbose = FALSE)
clus20.markers <- FindConservedMarkers(immune.combined, ident.1 = 20, grouping.var = "orig.ident", verbose = FALSE)
## Data visulization of markers used for annotations of clusters -- CHANGE THESE ACCORDINGLY TO FIT YOUR DATA
FeaturePlot(immune.combined, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A", "CCL2", "PPBP"), min.cutoff = "q9")

immune.combined <- RenameIdents(immune.combined, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T",
                                `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "NK", `7` = "T activated", `8` = "DC", `9` = "B Activated",
                                `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets", `14` = "HSPC")
DimPlot(immune.combined, label = TRUE)

immune.combined <- RenameIdents(immune.combined, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T",
                                `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "NK", `7` = "T activated", `8` = "DC", `9` = "B Activated",
                                `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets", `14` = "HSPC")
DimPlot(immune.combined, label = TRUE)

# Identification of differentially expressed genes across conditions
Idents(immune.combined) = "seurat_clusters"
immune.combined$celltype.treatment <- paste(Idents(immune.combined), immune.combined$orig.ident, sep = "_")
Idents(immune.combined) <- "celltype.treatment"

clus0.deg <- FindMarkers(immune.combined, ident.1 = "0_Control", ident.2 = "0_Infected", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1.5))
clus1.deg <- FindMarkers(immune.combined, ident.1 = "1_Control", ident.2 = "1_Infected", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1.5))
clus2.deg <- FindMarkers(immune.combined, ident.1 = "2_Control", ident.2 = "2_Infected", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1.5))
clus3.deg <- FindMarkers(immune.combined, ident.1 = "3_Control", ident.2 = "3_Infected", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1.5))
clus4.deg <- FindMarkers(immune.combined, ident.1 = "4_Control", ident.2 = "4_Infected", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1.5))
clus5.deg <- FindMarkers(immune.combined, ident.1 = "5_Control", ident.2 = "5_Infected", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1.5))
clus6.deg <- FindMarkers(immune.combined, ident.1 = "6_Control", ident.2 = "6_Infected", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1.5))
clus7.deg <- FindMarkers(immune.combined, ident.1 = "7_Control", ident.2 = "7_Infected", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1.5))
clus8.deg <- FindMarkers(immune.combined, ident.1 = "8_Control", ident.2 = "8_Infected", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1.5))
clus9.deg <- FindMarkers(immune.combined, ident.1 = "9_Control", ident.2 = "9_Infected", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1.5))
clus10.deg <- FindMarkers(immune.combined, ident.1 = "10_Control", ident.2 = "10_Infected", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1.5))
clus11.deg <- FindMarkers(immune.combined, ident.1 = "11_Control", ident.2 = "11_Infected", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1.5))
clus12.deg <- FindMarkers(immune.combined, ident.1 = "12_Control", ident.2 = "12_Infected", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1.5))
clus13.deg <- FindMarkers(immune.combined, ident.1 = "13_Control", ident.2 = "13_Infected", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1.5))
clus14.deg <- FindMarkers(immune.combined, ident.1 = "14_Control", ident.2 = "14_Infected", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1.5))
clus15.deg <- FindMarkers(immune.combined, ident.1 = "15_Control", ident.2 = "15_Infected", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1.5))
clus16.deg <- FindMarkers(immune.combined, ident.1 = "16_Control", ident.2 = "16_Infected", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1.5))
clus17.deg <- FindMarkers(immune.combined, ident.1 = "17_Control", ident.2 = "17_Infected", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1.5))
clus18.deg <- FindMarkers(immune.combined, ident.1 = "18_Control", ident.2 = "18_Infected", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1.5))
clus19.deg <- FindMarkers(immune.combined, ident.1 = "19_Control", ident.2 = "19_Infected", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1.5))
clus20.deg <- FindMarkers(immune.combined, ident.1 = "20_Control", ident.2 = "20_Infected", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1.5))
## Visualization methods for DEGs -- CHANGE THESE ACCORDINGLY TO FIT YOUR DATA
### Plotting LFC between two conditions
theme_set(theme_cowplot())
Idents(immune.combined) = "seurat_clusters"
clus0cells <- subset(immune.combined, subset = seurat_clusters == c("0","2"))
Idents(clus0.cells) <- "orig.ident"
avg.clus0.cells <- as.data.frame(log1p(AverageExpression(clus0.cells, verbose = FALSE)$RNA))
avg.clus0.cells$gene <- rownames(avg.clus0.cells)
genes.to.label = c("Isg15", "Ly6e", "Ifi6", "Isg20", "Mx1", "Ifit2", "Ifit1", "Cxcl10", "Ccl8", "Irf7")
plot = ggplot(avg.clus0.cells, aes(IRF7, WT)) + geom_point() + ggtitle("Cluster 0")
LabelPoints(plot = plot, points = genes.to.label, repel = TRUE)
### Expression on cluster UMAP and violin plot
FeaturePlot(immune.combined, features = c("Irf7"), split.by = "orig.ident", max.cutoff = 3, cols = c("grey", "red"))

### Comparing DGE of one cluster against all others
Idents(immune.combined) = "seurat_clusters"
Idents(immune.combined, WhichCells(object = immune.combined, idents = 2)) <- 'clus2.pos'
Idents(immune.combined, WhichCells(object = immune.combined, idents = c(0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20))) <- 'clus2.neg'
genes <- FindMarkers(immune.combined, ident.1 = 'clus2.pos', ident.2 = 'clus2.neg')
EnhancedVolcano(genes, lab = rownames(genes), x = 'avg_log2FC', y = 'p_val_adj')
# Can also do by expression
# Idents(pbmc_small, WhichCells(object = pbmc_small, expression = CD8A > 1, slot = 'data')) <- 'cd8.pos'
# Idents(pbmc_small, WhichCells(object = pbmc_small, expression = CD8A <= 1, slot = 'data')) <- 'cd8.neg'


### Highlighting groups
Idents(immune.combined) = "seurat_clusters"
c0 = WhichCells(immune.combined, idents="0")
c1 = WhichCells(immune.combined, idents="1")
c2 = WhichCells(immune.combined, idents="2")
c3 = WhichCells(immune.combined, idents="3")
c4 = WhichCells(immune.combined, idents="4")
c5 = WhichCells(immune.combined, idents="5")
c6 = WhichCells(immune.combined, idents="6")
c7 = WhichCells(immune.combined, idents="7")
c8 = WhichCells(immune.combined, idents="8")
c9 = WhichCells(immune.combined, idents="9")
c10 = WhichCells(immune.combined, idents="10")
c11 = WhichCells(immune.combined, idents="11")
c12 = WhichCells(immune.combined, idents="12")
c13 = WhichCells(immune.combined, idents="13")
c14 = WhichCells(immune.combined, idents="14")
c15 = WhichCells(immune.combined, idents="15")
c16 = WhichCells(immune.combined, idents="16")
c17 = WhichCells(immune.combined, idents="17")
c18 = WhichCells(immune.combined, idents="18")
c19 = WhichCells(immune.combined, idents="19")
c20 = WhichCells(immune.combined, idents="20")
DimPlot(immune.combined, cells.highlight = list(c1, c11, c18, c8, c6, c4, c16, c12, c15, c3, c0, c9, c14, c19, c7, c17, c5), 
        cols.highlight = c("green2", "green", "darkgreen", "maroon", "darkred", "red4", "red3", "red2", "orange3", "orange2", "orange",
                                   "blue4", "blue3", "blue2", "blue", "green3", "red"), 
                                   cols = "grey")

### Plotting average expression across clusters
Idents(immune.combined) = "seurat_clusters"
Tcells = subset(immune.combined, subset = seurat_clusters == c("1","11","18","8","6","4"))
NKcells = subset(immune.combined, subset = seurat_clusters == c("16","12","15","3"))
Bcells = subset(immune.combined, subset = seurat_clusters == c("0","9","14","19","2"))
myeloid = subset(immune.combined, subset = seurat_clusters == c("7","17","5"))
avgTcells = AverageExpression(Tcells, assay = "RNA", return.seurat=T)
avgNKcells = AverageExpression(NKcells, assay = "RNA", return.seurat=T)
avgBcells = AverageExpression(Bcells, assay = "RNA", return.seurat=T)
avgMyecells = AverageExpression(myeloid, assay = "RNA", return.seurat=T)
clus0.sort = clus0.markers[order(clus1.markers$Infected_avg_log2FC, decreasing = TRUE),]
clus1.sort = clus1.markers[order(clus1.markers$Control_avg_log2FC, decreasing = TRUE),]
clus2.sort = clus2.markers[order(clus1.markers$Infected_avg_log2FC, decreasing = TRUE),]
clus3.sort = clus3.markers[order(clus1.markers$Infected_avg_log2FC, decreasing = TRUE),]
clus4.sort = clus4.markers[order(clus1.markers$Infected_avg_log2FC, decreasing = TRUE),]
clus5.sort = clus5.markers[order(clus1.markers$Infected_avg_log2FC, decreasing = TRUE),]
clus6.sort = clus6.markers[order(clus1.markers$Infected_avg_log2FC, decreasing = TRUE),]
clus7.sort = clus7.markers[order(clus1.markers$Infected_avg_log2FC, decreasing = TRUE),]
clus8.sort = clus8.markers[order(clus1.markers$Infected_avg_log2FC, decreasing = TRUE),]
clus9.sort = clus9.markers[order(clus1.markers$Infected_avg_log2FC, decreasing = TRUE),]
clus11.sort = clus11.markers[order(clus1.markers$Infected_avg_log2FC, decreasing = TRUE),]
clus12.sort = clus12.markers[order(clus1.markers$Infected_avg_log2FC, decreasing = TRUE),]
clus14.sort = clus14.markers[order(clus1.markers$Infected_avg_log2FC, decreasing = TRUE),]
clus15.sort = clus15.markers[order(clus1.markers$Infected_avg_log2FC, decreasing = TRUE),]
clus16.sort = clus16.markers[order(clus1.markers$Infected_avg_log2FC, decreasing = TRUE),]
clus17.sort = clus17.markers[order(clus1.markers$Infected_avg_log2FC, decreasing = TRUE),]
clus18.sort = clus18.markers[order(clus1.markers$Infected_avg_log2FC, decreasing = TRUE),]
clus19.sort = clus19.markers[order(clus1.markers$Infected_avg_log2FC, decreasing = TRUE),]
DoHeatmap(avgTcells, features = c(head(rownames(clus1.sort, 5)), head(rownames(clus4.sort, 5)), head(rownames(clus6.sort, 5)), 
                                  head(rownames(clus8.sort, 5)), head(rownames(clus11.sort, 5)), head(rownames(clus18.sort, 5))), 
          slot = "data", draw.lines=F, group.bar.height=0) + scale_fill_gradientn(colors = c("gray", "red3", "darkred"))
DoHeatmap(avgNKcells, features = c(head(rownames(clus3.sort, 5)), head(rownames(clus12.sort, 5)), head(rownames(clus15.sort, 5)), 
                                   head(rownames(clus16.sort, 5))), 
          slot = "data", draw.lines=F, group.bar.height=0) + scale_fill_gradientn(colors = c("gray", "red3", "darkred"))
DoHeatmap(avgBcells, features = c(head(rownames(clus0.sort, 5)), head(rownames(clus2.sort, 5)), head(rownames(clus9.sort, 5)), 
                                  head(rownames(clus14.sort, 5)), head(rownames(clus19.sort, 5))), 
          slot = "data", draw.lines=F, group.bar.height=0) + scale_fill_gradientn(colors = c("gray", "red3", "darkred"))
DoHeatmap(avgMyecells, features = c(head(rownames(clus5.sort, 5)), head(rownames(clus7.sort, 5)), head(rownames(clus17.sort, 5))), 
          slot = "data", draw.lines=F, group.bar.height=0) + scale_fill_gradientn(colors = c("gray", "red3", "darkred"))


plot <- VlnPlot(immune.combined, features = c("LYZ", "ISG15", "CXCL10"), split.by = "orig.ident", group.by = "seurat_clusters", pt.size = 0.05, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

### Volcano plots
EnhancedVolcano(clus4.deg, lab = rownames(clus4.deg), x = 'avg_log2FC', y = 'p_val_adj') #log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group
EnhancedVolcano(clus4.deg, lab = rownames(clus4.deg), x = 'avg_log2FC', y = 'p_val_adj', selectLab = gene.list.hallmark$HALLMARK_INTERFERON_ALPHA_RESPONSE, title = 'IFNa Response', FCcutoff = 0.5, pointSize = 1.0, labSize = 6.0, drawConnectors = T)

### Built-in enrichR analysis
Idents(immune.combined) <- "orig.ident"
DEenrichRPlot(immune.combined, ident.1 = "Infected", ident.2 = "Control", max.genes = 100, enrich.database = "GO_Biological_Process_2021", num.pathway = 20)
DEenrichRPlot(immune.combined, ident.1 = "Infected", ident.2 = "Control", max.genes = 100, enrich.database = "KEGG_2019_Mouse", num.pathway = 20)
Idents(immune.combined) <- "celltype.treatment"
DEenrichRPlot(immune.combined, ident.1 = "11_Infected", ident.2 = "11_Control", max.genes = 50, enrich.database = "GO_Biological_Process_2021", num.pathway = 10)
DEenrichRPlot(immune.combined, ident.1 = "11_Infected", ident.2 = "11_Control", max.genes = 50, enrich.database = "KEGG_2019_Mouse", num.pathway = 10)


### DotPlot of annotation markers
Idents(immune.combined) = "seurat_clusters"
annotation.markers = as.list(c("Cd3e","Lat","Rag1","Cd8a","Mki67","Isg20","Ccr7","Ccr9","Ccl5","Nkg7","Gzmb","Cd79a","Ms4a1","Ighm","Ighd",'Sell',"Ly6a","Il7r","Cd74","Mzb1", "Igkc","Itgam","Ccr2","Ccl6","C1qa","Cd52","Siglech","Itgax","Irf7","Itga2b","Ppbp","Hbb-bt","Hba-a1"))
annotation.markers = as.character(annotation.markers)
levels(immune.combined) <- c("1","11","18","8","6","4","16","12","15","3","0","9","14","19","7","17","5","20","13","10","2")
#levels(immune.combined) <- c('T cell (immature)','T cell CD8 (Ki67+)','T cell CD8 (activated)','T cell CD8 (naïve)','T cell CD4 (naïve)','T cell CD4 (CCR9+)','NKT cell 1','NKT cell 2','NK cell 1','NK cell 2','B cell (naïve)','B cell (early activated)','B cell (Ly6a+)','B cell (immature)','Monocyte','Macrophage (M1)','Macrophage (M2)','pDC','Platelet','Erythrocyte','Fibroblast')
DotPlot(immune.combined, features = annotation.markers, cols = c("blue", "red"), dot.scale = 6) + RotatedAxis()


#Tabulating the number of cells
## Number of infected cells & control cells
table(paste0(immune.combined@meta.data$orig.ident))
## Total number of cells, i.e. infected + control cells, in each cluster
table(Idents(immune.combined))
## Number of control versus infected cells in each cluster
table(Idents(immune.combined), immune.combined$orig.ident)
## Number of doublets and singlets
table(paste0(immune.combined.with.doubs@meta.data$CellTypes_DF))

# Making plots of cluster markers
# Example below is for Iglc1; change gene name to gene of interest
Idents(immune.combined) = "seurat_clusters"

VlnPlot(immune.combined, features = c("Iglc1"), pt.size = 0, combine = FALSE) 
wrap_plots(plots = plots, ncol = 1)

FeaturePlot(immune.combined, features = c("Iglc1"), split.by = "orig.ident", max.cutoff = 3, cols = c("grey", "red"))
FeaturePlot(immune.combined, features = c("Iglc1"), max.cutoff = 3, cols = c("grey", "red"))

VlnPlot(immune.combined, features = c("Irf7"), pt.size = .005, split.by="orig.ident", combine = FALSE, assay="RNA") #split by genotype
VlnPlot(immune.combined, features = c("FCER1A"), pt.size = 0.005, combine = FALSE, sort=T, assay="RNA") #all clusters and sorted
FeaturePlot(immune.combined, features = c("Irf7"), split.by = "orig.ident", min.cutoff = 0.02, cols = c("grey", "red"))

# FeaturePlot(immune.combined, features = c("Cd3e"), min.cutoff = 0.02, cols = c("grey", "red"), reduction = "umap") + VlnPlot(immune.combined, features = c("Cd3e"), pt.size = 0.005, combine = FALSE, sort=T, assay="RNA")


#-------------------------------------------------------------------------------------------------------------------------------
#GSEA USING ESCAPE

#Setting up the single cell experiment object and acquiring gene sets
msigdbr::msigdbr_collections() #view what gene set collections are available for download
msigdbr::msigdbr_species() #view which species are available
gs.immune <- getGeneSets(library="H", species="Mus musculus")
gene.list.hallmark <- GSEABase::geneIds(gs.immune)

#Enrichment on RNA count data
es.seurat <- enrichIt(obj = sce, gene.sets = gs.immune, groups = 400, cores = 2)
immune.combined.gsea <- AddMetaData(immune.combined, es.seurat)

#Visualization
dittoHeatmap(immune.combined.gsea, genes = NULL, metas = names(es.seurat),
             annot.by = c("seurat_clusters", "orig.ident"),
             cluster_cols = FALSE, order.by = c("orig.ident", "seurat_clusters"),
             fontsize = 7)
dittoHeatmap(immune.combined.gsea, genes = NULL, 
             metas = c("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
                       "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_COMPLEMENT", "HALLMARK_APOPTOSIS"),
             annot.by = c("seurat_clusters", "orig.ident"),
             cluster_cols = FALSE, order.by = c("orig.ident", "seurat_clusters"),
             fontsize = 7)
multi_dittoPlot(immune.combined.gsea, vars = c("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
                                               "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_COMPLEMENT", "HALLMARK_APOPTOSIS"), 
                group.by = "orig.ident", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Enrichment Scores", 
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))
enrichmentPlot(immune.combined.gsea, gene.set = "HALLMARK_TNFA_SIGNALING_VIA_NFKB", gene.sets = "gs.immune", group = "orig.ident")

##Look at global differences by cluster in individual genes in a gene set
Idents(immune.combined.gsea) = "seurat_clusters"
DotPlot(immune.combined.gsea, features = gene.list.hallmark$HALLMARK_INTERFERON_ALPHA_RESPONSE, cols = c("blue", "red"), dot.scale = 6, split.by = "orig.ident", group.by="seurat_clusters") + RotatedAxis()
##Look at global differences in individual genes in a gene set
Idents(immune.combined.gsea) = "orig.ident"
DotPlot(immune.combined.gsea, features = gene.list.hallmark$HALLMARK_INTERFERON_ALPHA_RESPONSE, cols = c("blue", "red"), dot.scale = 6) + RotatedAxis() #can then pull outshort list of interesting genes based on global differences
short.list.ifna <- as.character(c("B2m","Cd47","Cd74","Cnp","Elf1","Epsti1","H2-D1","H2-O7","H2-T23","Irf2","Ly6e","Psma3","Psmb8","Psmb9","Psme1","Psme2","Txnip"))
##Look at cluster/sample differences in indiviudal genes in a gene set
VlnPlot(immune.combined.gsea, features = "Cd47", split.by = "orig.ident", group.by = "seurat_clusters", pt.size = 0.05, combine = FALSE)
Idents(immune.combined.gsea) = "seurat_clusters"
DotPlot(immune.combined.gsea, features = short.list.ifna, cols = c("blue", "red"), dot.scale = 6, split.by = "orig.ident", group.by="seurat_clusters") + RotatedAxis()
## If need to remove genes not found then use the following
# test <- gene.list.hallmark$HALLMARK_INTERFERON_ALPHA_RESPONSE
# test <- test[! test %in% c("C1s2", "Cxcl11", "Ifi27l2b", "Ifi44l", "Lamp3", "Mx2")]

ridgeEnrichment(immune.combined.gsea[[]], gene.set = "HALLMARK_INFLAMMATORY_RESPONSE", group = "seurat_clusters", facet = "orig.ident", add.rug = TRUE)

#Plotting differential expression of genes in a gene set
cytokines = read.csv("cytokines.csv")
DefaultAssay(immune.combined) = "RNA"
Idents(immune.combined) = "orig.ident"
all.markers <- FindMarkers(immune.combined, ident.1 = "Control", ident.2 = "Infected", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1.5))
cytokine.markers = subset(all.markers, rownames(all.markers) %in% cytokines[,1])
cytokine.markers <- cytokine.markers %>% mutate(pos = avg_log2FC>0)
ggplot(cytokine.markers,aes(fill=pos)) + geom_col(aes(x=rownames(cytokine.markers),y=avg_log2FC),position="identity") + coord_flip() + theme_bw() + 
  labs(x=NULL, y="Avg Difference L2FC") + 
  theme(text=element_text(family="Helvetica", size=20), panel.grid.minor = element_blank(), panel.border = element_blank()) + 
  scale_fill_manual(values=c("red3","green3")) + theme(legend.position="none")
test = subset(all.markers, rownames(all.markers) == "Ifit1")
head(test)

DoHeatmap(immune.combined, features = gene.list.hallmark$HALLMARK_INTERFERON_ALPHA_RESPONSE, group.by = "orig.ident", slot = "data")

PCA <- performPCA(enriched = immune.combined.gsea[[]], gene.sets = names(gs.immune), groups = c("orig.ident", "seurat_clusters"))
pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = TRUE)
pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = FALSE, facet = "cluster") 
masterPCAPlot(immune.combined.gsea[[]], gene.sets = names(gs.immune), PCx = "PC1", PCy = "PC2", top.contribution = 10)

#---------------------

#GSEA using fgsea
### Filter out unnecessary genes to make the process run faster
genes = wilcoxauc(immune.combined, 'seurat_clusters')
dplyr::count(genes, group) #gives number of genes in each cluster

### Download gsea lists
msigdbr::msigdbr_collections() #view what gene set collections are available for download
msigdbr::msigdbr_species() #view which species are available
gs.immune <- getGeneSets(library="H", species="Mus musculus")
gene.list.hallmark <- GSEABase::geneIds(gs.immune)

### Total GSEA
total.gsea = genes %>% arrange(desc(auc)) %>% dplyr::select(feature, auc)
total.ranks = tibble::deframe(total.gsea)

fgsea = fgsea(gene.list.hallmark, stats = total.ranks, nperm=1000)
fgseatidy = fgsea %>% as_tibble() %>% arrange(desc(NES))

ggplot(fgseatidy %>% filter(padj < 0.008) %>% head(n= 50), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark Pathways NES from GSEA", size=10) + 
  theme_minimal() +
  scale_fill_manual(values = c("green3", "red3"))

### Split cluster of interest out into gene-level statistics
clus11.gsea = genes %>% dplyr::filter(group == "11") %>% arrange(desc(auc)) %>% dplyr::select(feature, auc)
ranks = tibble::deframe(clus11.gsea)

### Calculate GSEA
fgsea11 = fgsea(gene.list.hallmark, stats = ranks, nperm=1000)
fgsea11tidy = fgsea11 %>% as_tibble() %>% arrange(desc(NES))
fgsea11tidy %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj) %>% head()

### Make plots
ggplot(fgsea11tidy %>% filter(padj < 0.008) %>% head(n= 50), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal() +
  scale_fill_manual(values = c("green", "red"))

plotEnrichment(gene.list.hallmark[["HALLMARK_E2F_TARGETS"]],
               ranks) + labs(title="HALLMARK_E2F_TARGETS")


