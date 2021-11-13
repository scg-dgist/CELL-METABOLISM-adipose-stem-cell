##### Code for
##### Author:  Eun Seo Park (evergreen@dgist.ac.kr)
###Load Data
library(devtools)
library(harmony)
library(dplyr)
library(Seurat)
sapply(list.files("D:/OneDrive - dgist.ac.kr/Function/",full.names = T),source)
NCD_seurat <- readRDS("D:/OneDrive - dgist.ac.kr/Adipocyte_Precursor/paper_final/NCD_Seurat.rds")

#Check Batch Errct
library(cowplot)
DimPlot(NCD_seurat, reduction = 'umap', label = F, label.size = 7.5, group.by = 'region')
fig1=DimPlot(NCD_seurat, reduction = 'umap', label = T, label.size = 7.5)
fig2=DimPlot(NCD_seurat, reduction = 'umap', label = T, label.size = 7.5, split.by = 'region')
plot_grid(fig1,fig2)
 
#Batch correction using Harmony package
set.seed(123456)
NCD_seurat_H <- NCD_seurat %>% 
  RunHarmony("region", plot_convergence = F, max.iter.harmony = 100)

NCD_seurat_H <- NCD_seurat_H %>% 
  RunUMAP(reduction = "harmony", dims = 1:15) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:15) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()


fig1=DimPlot(NCD_seurat_H, reduction = 'umap', label = T, label.size = 7.5)
fig2=DimPlot(NCD_seurat_H, reduction = 'umap', label = T, label.size = 7.5, split.by = 'region')



