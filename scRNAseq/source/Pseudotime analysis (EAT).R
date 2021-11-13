##### Code for
##### Author:  Eun Seo Park (evergreen@dgist.ac.kr)
###Load Data
library(monocle)
library(scater)
library(scran)
library(Seurat)
seurat <- readRDS("I:/서울대_infla/EAT_palantir_scenic_result_seurat2.rds")

#Convert single cell experiment object
sce <- as.SingleCellExperiment(seurat)
set.seed(123)
clusters <- quickCluster(sce, method="igraph")
table(clusters)
sce <- computeSumFactors(sce, cluster=clusters)
sizeFactors(sce)
seurat$sizefactor <- sizeFactors(sce)

#Make cell dataset
hvgGenes <- VariableFeatures(seurat)
rawcounts <- seurat@assays$RNA@counts
newclusters <- Idents(seurat)

pd.df = seurat@meta.data
pd =new('AnnotatedDataFrame', data =pd.df)
fd.df = data.frame(gene_short_name=rownames(rawcounts), row.names = rownames(rawcounts))
fd =new('AnnotatedDataFrame', data =fd.df)
monosetFiltered<- newCellDataSet(as.matrix(rawcounts),
                                 phenoData = pd,
                                 featureData = fd,
                                 expressionFamily=negbinomial.size())

#Infer trajectory using HVG
monosetFiltered <- setOrderingFilter(monosetFiltered, hvgGenes)

monosetFiltered@phenoData@data$Size_Factor = sizeFactors(sce)
monosetFiltered <- estimateDispersions(monosetFiltered)
monosetFiltered <- reduceDimension(monosetFiltered, residualModelFormulaStr = ~nCount_RNA+nFeature_RNA)
monosetFiltered <- orderCells(monosetFiltered,root_state = 1)
pData(monosetFiltered)$cluster = newclusters

#Plot results
plot_cell_trajectory(monosetFiltered, color_by = "seurat_clusters", cell_size = 1,
                     show_state_number = TRUE)
plot_cell_trajectory(monosetFiltered, color_by = "diet", cell_size = 1)
plot_cell_trajectory(monosetFiltered, color_by = "state", cell_size = 1)
monocle_dim <- t(as.matrix(monosetFiltered@reducedDimS))
colnames(monocle_dim) <- c("monocle1","monocle2")

#Add monocle2 results to Seurat object

seurat[['monocle']] <- CreateDimReducObject(embeddings = monocle_dim, key = "monocle", assay = DefaultAssay(seurat))
DimPlot_eunseo(seurat, group.by = "seurat_clusters", reduction = "monocle",
               label = F, path = "EAT_monocle2_batch_corrected_seurat_cluster.png")
DimPlot_eunseo(seurat, group.by = "seurat_clusters",label = T)

