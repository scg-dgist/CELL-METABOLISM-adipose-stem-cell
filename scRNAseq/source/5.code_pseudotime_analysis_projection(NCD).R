##### Code for
##### Author:  Eun Seo Park (evergreen@dgist.ac.kr)
###Load EAT NCD Data
library(Seurat)
NCD_Seurat <- readRDS("D:/OneDrive - dgist.ac.kr/Adipocyte_Precursor/paper_final/NCD_Seurat.rds")
EAT_NCD_seurat <- readRDS("D:/OneDrive - dgist.ac.kr/Adipocyte_Precursor/paper_final/EAT_NCD_seurat.rds")
EAT_NCD_seurat$annot <- NCD_seurat$annot[colnames(EAT_NCD_seurat)]
#ADD Palantir TSNE reduction coordinate ON Seurat
EAT_NCD_palantir1 <- read.csv("D:/OneDrive - dgist.ac.kr/Adipocyte_Precursor/paper_final/palantir/EAT_NCD/EAT_NCD_tsne_200_300.csv",row.names = 1)
colnames(EAT_NCD_palantir1) <- c("Palantir1","Palantir2")

EAT_NCD_palantir1 <- as.matrix(EAT_NCD_palantir1)
EAT_NCD_seurat[['palantir1']] <- CreateDimReducObject(embeddings = EAT_NCD_palantir1, key = "Palantir_I", assay = DefaultAssay(EAT_NCD_seurat))
DimPlot(EAT_NCD_seurat, reduction = "palantir1", group.by = "annot",
        label = T)

#ADD Pseudotime to metadata
pseudotime <- read.csv("D:/OneDrive - dgist.ac.kr/Adipocyte_Precursor/paper_final/palantir/EAT_NCD/EAT_NCD_pseudotime_200.csv",row.names = 1, header = F)
EAT_NCD_seurat$pseudotime <- pseudotime[colnames(EAT_NCD_seurat),"V2"]

###Load IAT NCD Data
library(Seurat)
IAT_NCD_seurat <- readRDS("D:/OneDrive - dgist.ac.kr/Adipocyte_Precursor/paper_final/IAT_NCD_seurat.rds")
IAT_NCD_seurat$annot <- NCD_seurat$annot[colnames(IAT_NCD_seurat)]
IAT_NCD_seurat$stage <- NCD_seurat$stage[colnames(IAT_NCD_seurat)]

#ADD Palantir TSNE reduction coordinate ON Seurat
IAT_NCD_palantir1 <- read.csv("D:/OneDrive - dgist.ac.kr/Adipocyte_Precursor/paper_final/palantir/IAT_NCD/IAT_NCD_tsne_100_200.csv",row.names = 1)
colnames(IAT_NCD_palantir1) <- c("Palantir1","Palantir2")
IAT_NCD_palantir1 <- as.matrix(IAT_NCD_palantir1)
IAT_NCD_seurat[['palantir1']] <- CreateDimReducObject(embeddings = IAT_NCD_palantir1, key = "Palantir_I", assay = DefaultAssay(IAT_NCD_seurat))
DimPlot(IAT_NCD_seurat, reduction = "palantir1", group.by = "stage",
        label = T)
DimPlot(IAT_NCD_seurat, 
        label = T)
seurat_visualization(IAT_NCD_seurat,reduction = "palantir1",
                     genelist = c("Fabp4","Bmp7","Il6","Cebpd"))
#ADD Pseudotime to metadata
pseudotime <- read.csv("D:/OneDrive - dgist.ac.kr/Adipocyte_Precursor/paper_final/palantir/IAT_NCD/IAT_NCD_pseudotime_100.csv",row.names = 1, header = F)
IAT_NCD_seurat$pseudotime <- pseudotime[colnames(IAT_NCD_seurat),"V2"]

#ADD path probability to metadata
Branch_Probs <- read.csv("D:/OneDrive - dgist.ac.kr/Adipocyte_Precursor/paper_final/palantir/IAT_NCD/IAT_NCD_branch_probs_100.csv",row.names = 1)
IAT_NCD_seurat$Path1 <- Branch_Probs[colnames(IAT_NCD_seurat),"CATCCCAAGGAGTCTG.H9"]
IAT_NCD_seurat$Path2 <- Branch_Probs[colnames(IAT_NCD_seurat),"GGTTCTCAGTATGGAT.H9"]

library(Seurat)
#load Seurat
EAT_NCD_seurat <- readRDS("D:/OneDrive - dgist.ac.kr/Adipocyte_Precursor/paper_final/EAT_NCD_seurat.rds")
IAT_NCD_seurat <- readRDS("D:/OneDrive - dgist.ac.kr/Adipocyte_Precursor/paper_final/IAT_NCD_seurat.rds")
IAT_NCD_seurat <- ScaleData(IAT_NCD_seurat, features = rownames(IAT_NCD_seurat))
#Palantir load 
IAT_NCD_palantir <- as.matrix(IAT_NCD_seurat@reductions$palantir1@cell.embeddings)


library(SingleCellExperiment)
library(KernelKnn)

### Select Highly variable genes (feature selection)
hvg <- VariableFeatures(IAT_NCD_seurat)
hvg_mat <- as.matrix(IAT_NCD_seurat@assays$RNA@data)[hvg,]
hvg_test_mat <- as.matrix(EAT_NCD_seurat@assays$RNA@data[hvg,])

#Find nearest cells of EAT cells from IAT cells

indexN <- KernelKnn::knn.index.dist(t(hvg_mat),t(hvg_test_mat),k=10,threads = 4,
                                    method = "pearson_correlation")
iN2 <- indexN$test_knn_idx
rownames(iN2) <- colnames(hvg_test_mat)

#Average palantir tsne coordinate of cells
iN3 <- apply(iN2, 2, function(x) colnames(hvg_mat)[x])
idpalantirx <- apply(iN3, 2, function(x) IAT_NCD_palantir[x,1])
idpalantiry <- apply(iN3, 2, function(x) IAT_NCD_palantir[x,2])
prjpalantirx <- rowMeans(idpalantirx)
prjpalantiry <- rowMeans(idpalantiry)

names(prjpalantirx) <- rownames(iN2)
names(prjpalantiry) <- rownames(iN2)
En_EAT_cellname <- colnames(EAT_NCD_seurat)[EAT_NCD_seurat$annot %in% c("E1","E2","E3","E4")]


df_test=data.frame(x=rbind(as.matrix(IAT_NCD_palantir[,1]),as.matrix(prjpalantirx[En_EAT_cellname])), 
                   y=rbind(as.matrix(IAT_NCD_palantir[,2]),as.matrix(prjpalantiry[En_EAT_cellname])), 
                   expression= c(rep("IAT_NCD",length(colnames(IAT_NCD_seurat))),EAT_NCD_seurat$annot[En_EAT_cellname]))

# Load the "scales" package
require(scales)
# Create vector of default ggplot2 colors
EAT_annot_col <- gg_color_hue(11)[c(1,2,3,4)]
library(ggplot2)

ggplot(df_test,
       aes(x=x, y=y, color = expression)) +
  geom_point(size=0.8) +
  scale_color_manual(values = c(EAT_annot_col,"grey79"))+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())





