# Adapted from https://github.com/Amandahsr/BulkTissueDeconvolution/tree/master/Guided%20Clustering%20(Seurat)
library(Seurat)
library(Matrix)

#Read normalised 10X scRNA-seq dataset. The scRNA-seq data contains 287,269 cells from the 4 chambers of normal human heart samples (LV, LA, RV, RA) and 33,694 genes.
sc.counts <- Read10X('C:\\Users\\Admin\\Documents\\Github\\ABC\\data')
#Sample selection step: For the purpose of the project, retain only cells from the LV.
LV.cells <- sc.counts[,substr(colnames(sc.counts), 1, 2) == "LV"]

#Pre-process genes to ensure we only cluster genes with meaningful data. Retain genes that meet two criteria: 
#1) Has counts in at least 10 cells and 2) Belongs to the top 30% highly expressed genes.
LV.cells <- LV.cells[rowSums(LV.cells != 0) >= 10, ,drop = FALSE]

#Due to RAM limitations, only retain the top 30% highly-expressed genes.  
LV.cells <- LV.cells[rev(order(rowSums(LV.cells))),]
LV.cells <- LV.cells[1:7494,]

#gene names extracted from SNPs left ventricle eQTL file and to retain only the genes related to SNPs.
gene.names.snps<-read.csv("gene.names.csv")
LV.cells<-LV.cells[rownames(LV.cells)%in%gene.names.snps[[1]],]

#Initialise UMI counts as Seurat object. min.cells and min.features are set to 0 because pre-processing was already performed.
seurat.data <- CreateSeuratObject(LV.cells, project = 'Seurat Guided Clustering', min.cells = 0, min.features = 0)

#Select highly variable genes with average expression set between 0.03 and 5, and z-score >= 5. This highlights meaningful biological signals during clustering.
variable.genes <- FindVariableFeatures(seurat.data, selection.method = 'mvp', mean.cutoff = c(0.03, 5), dispersion.cutoff = c(0.5, Inf))
gene.names <- rownames(variable.genes)


scaled.data <- ScaleData(variable.genes, features = gene.names)

#Linear dimensional reduction step: PCA is performed to filter out noise by extracting the first 50 PCs for guided clustering.
pca <- RunPCA(scaled.data, features = VariableFeatures(object = variable.genes))

#Cluster cells by creating a KNN graph and optimizing it.
clusters <- FindNeighbors(pca)
clusters <- FindClusters(clusters, resolution = 1)

umap <- RunUMAP(clusters, dims = 1:50, n.neighbors = 25, min.dist = 1.5, spread = 3.5)
DimPlot(umap, reduction = "umap", label = TRUE)

#DEG analysis: Find genes that are preferentially expressed in each cluster with p-value <= 0.05.
clusters.markers <- FindAllMarkers(umap, only.pos=TRUE, return.thresh = 0.05)

#Feature plots: Visualise cell-type markers to identify cell types. A subset of markers used to identify cell-types are shown below.
Cardiomyocyte.markers <- FeaturePlot(umap, features = c('MYH7', 'TTN', 'TNNT2', 'RYR2'))
Macrophage.markers <- FeaturePlot(umap, features = c('CD163', 'MRC1', 'COLEC12', 'MARCH1', 'SLC11A1', 'RBPJ', 'F13A1'))
SmoothMuscleCell.markers <- FeaturePlot(umap, features = c('MYH11', 'LMOD1'))
Fibroblast.markers <- FeaturePlot(umap, features = c('POSTN', 'IGF1', 'ADAMTS4', 'DCN', 'ELN', 'NOX4', 'VCAN'))
Endothelial.markers <- FeaturePlot(umap, features = c('PECAM1', 'PROX1', 'VWF'))
Pericyte.markers <- FeaturePlot(umap, features = c('PDGFRB', 'STEAP4'))
TCell.markers <- FeaturePlot(umap, features = c('PTPRC', 'SKAP1'))
BCell.markers <- FeaturePlot(umap, features = c('BANK1'))
NeuronalCell.markers <- FeaturePlot(umap, features = c('NRXN1', 'NRXN3', 'NCAM2'))


