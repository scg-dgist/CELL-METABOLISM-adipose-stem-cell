##### Code for
##### Author:  Eun Seo Park (evergreen@dgist.ac.kr)
library("scran")
library("scater")

#load data
H9 <- readRDS("/home/espark/project/mouse_iAT_eAT/data/dropletfiltered_mitox/dropletfilterH9set.rds")
H10 <- readRDS("/home/espark/project/mouse_iAT_eAT/data/dropletfiltered_mitox/dropletfilterH10set.rds")
H11 <- readRDS("/home/espark/project/mouse_iAT_eAT/data/dropletfiltered_mitox/dropletfilterH11set.rds")
H12 <- readRDS("/home/espark/project/mouse_iAT_eAT/data/dropletfiltered_mitox/dropletfilterH12set.rds")
A1 <- readRDS("/home/espark/project/mouse_iAT_eAT/data/dropletfiltered_mitox/dropletfilterA1set.rds")
A2 <- readRDS("/home/espark/project/mouse_iAT_eAT/data/dropletfiltered_mitox/dropletfilterA2set.rds")

#label sample name
colnames(H9) <- paste0(substring(colnames(H9),1,16),"-H9")
colnames(H10) <- paste0(substring(colnames(H10),1,16),"-H10")
colnames(H11) <- paste0(substring(colnames(H11),1,16),"-H11")
colnames(H12) <- paste0(substring(colnames(H12),1,16),"-H12")
colnames(A1) <- paste0(substring(colnames(A1),1,16),"-A1")
colnames(A2) <- paste0(substring(colnames(A2),1,16),"-A2")

rowData(H9) <- rowData(H9)[,-3]
rowData(H10) <- rowData(H10)[,-3]
rowData(H11) <- rowData(H11)[,-3]
rowData(H12) <- rowData(H12)[,-3]
rowData(A1) <- rowData(A1)[,-3]
rowData(A2) <- rowData(A2)[,-3]

mcols(H9) <- NULL
mcols(H10) <- NULL
mcols(H11) <- NULL
mcols(H12) <- NULL
mcols(A1) <- NULL
mcols(A2) <- NULL

#merge sample
sce_merge <- cbind(H9,H10,H11,H12,A1,A2)
### filtering genes
library(Matrix)
# keep_feature <- rowMeans(sce_merge@assays$data$logcounts)!=0
keep_feature <- rowSums(counts(sce_merge) != 0) > 3
sce_merge <- sce_merge[keep_feature, ]

#Normalization
library(scran)
set.seed(12345)
clusters <- quickCluster(sce_merge, method="igraph")
table(clusters)
sce_merge <- computeSumFactors(sce_merge, cluster=clusters)

library(scater)
sce_merge <- logNormCounts(sce_merge)

#Feature Selection
dec <- modelGeneVar(sce_merge)
hvg <- getTopHVGs(dec, fdr.threshold = 0.05)


library(Seurat)
first_seurat <- as.Seurat(sce_merge)
first_seurat@assays$RNA@var.features <- hvg

#convert ensembl gene name to Gene symbol
load("/home/espark/project/mouse_iAT_eAT/data/qcfiltered/ensemblGenes2019-07-11.RData")
rownames(first_seurat@assays$RNA@counts) <- uniquifyFeatureNames(rownames(first_seurat@assays$RNA@counts),
                                                           ensemblGenes[rownames(first_seurat@assays$RNA@counts),"external_gene_name"])
rownames(first_seurat@assays$RNA@data) <- uniquifyFeatureNames(rownames(first_seurat@assays$RNA@data),
                                                         ensemblGenes[rownames(first_seurat@assays$RNA@data),"external_gene_name"])
first_seurat@assays$RNA@var.features <- uniquifyFeatureNames(first_seurat@assays$RNA@var.features,
                                                       ensemblGenes[first_seurat@assays$RNA@var.features,"external_gene_name"])

#Make Seurat object
first_seurat <- ScaleData(first_seurat)
first_seurat@reductions$PCA_coldata <- NULL
PCA = 15
first_seurat<- RunPCA(first_seurat, npcs = 50)
first_seurat<- FindNeighbors(first_seurat, reduction="pca", dims= 1:PCA)
first_seurat<- FindClusters(first_seurat, reduction="pca", dims= 1:PCA)
first_seurat<- RunUMAP(first_seurat, dims = 1:PCA, seed.use = 42)


###Remove immune and mesothelial cells
rm_immune_subset <- subset(first_seurat, ident=c(14,16),invert=TRUE)
rm_immune_sce <- as.SingleCellExperiment(rm_immune_subset)

### filtering genes
library(Matrix)
# keep_feature <- rowMeans(sce_merge@assays$data$logcounts)!=0
library(scran)
set.seed(123)
clusters <- quickCluster(rm_immune_sce, method="igraph")
table(clusters)
rm_immune_sce <- computeSumFactors(rm_immune_sce, cluster=clusters)
library(scater)
rm_immune_sce <- logNormCounts(rm_immune_sce)

dec <- modelGeneVar(rm_immune_sce)
hvg <- getTopHVGs(dec, fdr.threshold = 0.05)
library(Seurat)

rm_immune_seurat <- as.Seurat(rm_immune_sce)
rm_immune_seurat@assays$RNA@var.features <- hvg
rm_immune_seurat <- ScaleData(rm_immune_seurat)
rm_immune_seurat@reductions$PCA_coldata <- NULL
rm_immune_seurat@reductions$pca <- NULL
rm_immune_seurat@reductions$umap <- NULL
rm_immune_seurat@reductions$PCA <- NULL
rm_immune_seurat@reductions$UMAP <- NULL

PCA <- 20
rm_immune_seurat<- RunPCA(rm_immune_seurat, npcs = PCA)
rm_immune_seurat<- FindNeighbors(rm_immune_seurat, reduction="pca", dims= 1:PCA)
rm_immune_seurat<- FindClusters(rm_immune_seurat, reduction="pca", dims= 1:PCA)
rm_immune_seurat<- RunUMAP(rm_immune_seurat, dims = 1:PCA, seed.use = 42)

#Add sample information to metadata
rm_immune_seurat$sample <- substring(colnames(rm_immune_seurat),18)
rm_immune_seurat$region <- "IAT"
rm_immune_seurat$region[rm_immune_seurat$sample %in% c("A2","H10","H12")] <- "EAT"
rm_immune_seurat$diet <- "NCD"
rm_immune_seurat$diet[rm_immune_seurat$sample %in% c("H11","H12")] <- "stHFD"
rm_immune_seurat$diet[rm_immune_seurat$sample %in% c("A1","A2")] <- "ltHFD"

saveRDS(rm_immune_seurat,"all_merge_rm_immune_mseo_seurat.rds")

#Make EAT seurat

EAT_subset <- subset(rm_immune_seurat, cells=colnames(rm_immune_seurat)[rm_immune_seurat$region=="EAT"])
EAT_sce <- as.SingleCellExperiment(EAT_subset)


### filtering genes
library(Matrix)
# keep_feature <- rowMeans(sce_merge@assays$data$logcounts)!=0

#Normalization
library(scran)
set.seed(123)
clusters <- quickCluster(EAT_sce, method="igraph")
table(clusters)
EAT_sce <- computeSumFactors(EAT_sce, cluster=clusters)
library(scater)
EAT_sce <- logNormCounts(EAT_sce)

#Feature Selection
dec <- modelGeneVar(EAT_sce)
hvg <- getTopHVGs(dec,fdr.threshold = 0.05)

library(Seurat)
EAT_seurat <- as.Seurat(EAT_sce)
EAT_seurat@assays$RNA@var.features <- hvg
EAT_seurat <- ScaleData(EAT_seurat)
EAT_seurat@reductions$PCA_coldata <- NULL
EAT_seurat@reductions$pca <- NULL
EAT_seurat@reductions$umap <- NULL
EAT_seurat@reductions$PCA <- NULL
EAT_seurat@reductions$UMAP <- NULL
PCA = 10
EAT_seurat<- RunPCA(EAT_seurat, npcs = PCA)
EAT_seurat<- FindNeighbors(EAT_seurat, reduction="pca", dims= 1:PCA)
EAT_seurat<- FindClusters(EAT_seurat, reduction="pca", dims= 1:PCA)
EAT_seurat<- RunUMAP(EAT_seurat, dims = 1:PCA, seed.use = 42)

#Remove testis cells in EAT samples

seurat_subset_recluster <- function(seurat_obj,cells=NULL, idents=NULL,invert=FALSE,
                                    PCA=15,n = NULL, prop = NULL,
                                    var.threshold = 0,  fdr.threshold = NULL){
  library(Seurat)
  seurat_subset <- subset(seurat_obj, cells=cells, idents=idents,invert=invert)
  sce <- as.SingleCellExperiment(seurat_subset)
  library(scran)
  set.seed(123)
  clusters <- quickCluster(sce, method="igraph")
  table(clusters)
  sce <- computeSumFactors(sce, cluster=clusters)
  
  library(scater)
  
  sce <- logNormCounts(sce)
  dec <- modelGeneVar(sce)
  hvg <- getTopHVGs(dec, n=n,
                    var.threshold = var.threshold, fdr.field = "FDR",
                    fdr.threshold = fdr.threshold)
  
  print(paste0("The number of hvg is ",length(hvg)))
  
  library(Seurat)
  seurat <- as.Seurat(sce)
  seurat@assays$RNA@var.features <- hvg
  seurat <- ScaleData(seurat)
  
  seurat@reductions$PCA_coldata <- NULL
  seurat@reductions$pca <- NULL
  seurat@reductions$umap <- NULL
  seurat@reductions$PCA <- NULL
  seurat@reductions$UMAP <- NULL
  
  PCA = PCA
  
  seurat<- RunPCA(seurat, npcs = PCA)
  seurat<- FindNeighbors(seurat, reduction="pca", dims= 1:PCA)
  seurat<- FindClusters(seurat, reduction="pca", dims= 1:PCA)
  seurat<- RunUMAP(seurat, dims = 1:PCA, seed.use = 42)
  
  return(seurat)
  
}

rm_testis_EAT_seurat <- seurat_subset_recluster(seurat, idents=15,invert=TRUE, fdr.threshold=0.05)


NCD_seurat <- seurat_subset_recluster(rm_immune_seurat, cells=colnames(rm_immune_seurat)[rm_immune_seurat$diet %in% c("NCD")], n=1000,
                                      PCA=25)


