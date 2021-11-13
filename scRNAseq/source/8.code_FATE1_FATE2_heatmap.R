##### Code for
##### Author:  Eun Seo Park (evergreen@dgist.ac.kr)
###Load Data
library(Seurat)
EAT_seurat <- readRDS("D:/OneDrive - dgist.ac.kr/Adipocyte_Precursor/paper_final/EAT_palantir_scenic_seurat.rds")
branch_probs <- read.csv("D:/OneDrive - dgist.ac.kr/Adipocyte_Precursor/paper_final/EAT_batch_corrected_palantir/EAT_all_batch_corrected_branch_probs_300.csv",row.names = 1)
colnames(branch_probs) <- c("A", "B")

#Normalize path probability of cells
normWeights <- sweep(branch_probs, 1, FUN = "/",
                     STATS = apply(branch_probs, 1, sum))
head(normWeights)
set.seed(10)

#Multinomial sampling based on path probability
wSamp <- apply(normWeights, 1, 
               function(prob){rmultinom(n = 1, prob = prob, size = 1)})
head(wSamp)

if(is.null(dim(wSamp))){
  wSamp <- matrix(wSamp, ncol = 1)
} else{
  wSamp <- t(wSamp)
}

colnames(wSamp) <- colnames(normWeights)
meta <- cbind(EAT_seurat@meta.data, wSamp)

#Add path information to metadata
EAT_seurat@meta.data <- meta
EAT_seurat$path <- "No"
EAT_seurat$path[EAT_seurat$A==1] <- "A"
EAT_seurat$path[EAT_seurat$B==1] <- "B"


#Simplify pseudotime information
EAT_seurat$pseudotime_simple <- "Early"
EAT_seurat$pseudotime_simple[EAT_seurat$pseudotime>=0.35&EAT_seurat$pseudotime<0.75] <- "Intermediate"
EAT_seurat$pseudotime_simple[EAT_seurat$pseudotime>=0.75&EAT_seurat$pseudotime<=1.0] <- "Late"

test_seurat_bc <- EAT_seurat
Idents(test_seurat_bc) <- paste0(test_seurat_bc$pseudotime_simple,"_",test_seurat_bc$path)

#Find DEG between Two path
Early_gene_bc <- FindMarkers(test_seurat_bc, ident.1 = "Early_A",ident.2 = "Early_B",logfc.threshold = 0)
Inter_gene_bc <-FindMarkers(test_seurat_bc, ident.1 = "Intermediate_A",ident.2 = "Intermediate_B",logfc.threshold = 0)
Late_gene_bc <- FindMarkers(test_seurat_bc, ident.1 = "Late_A",ident.2 = "Late_B",logfc.threshold = 0)

sig_Inter_gene_bc <- subset(Inter_gene_bc, abs(Inter_gene_bc$avg_logFC)>0.25&Inter_gene_bc$p_val_adj<0.05)
sig_Late_gene_bc <- subset(Late_gene_bc,abs(Late_gene_bc$avg_logFC)>0.25&Late_gene_bc$p_val_adj<0.05)


check_seurat <- EAT_seurat
check_seurat$pseudotime[check_seurat$pseudotime>0.8] <-0.8
Idents(check_seurat) <- paste0(check_seurat$path,"_",round(check_seurat$pseudotime,1))

#Average gene expression per path
check_seurat <- ScaleData(check_seurat, features = rownames(check_seurat))
avg_gene <- AverageExpression(check_seurat, features = unique(c(rownames(sig_Inter_gene_bc),
                                                         rownames(sig_Late_gene_bc),c("Tgfb1")))
                              ,slot = "scale.data")

#Visualize
mat <- avg_gene$RNA
library(pheatmap)
library(RColorBrewer)
palette_length = 100
my_color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(palette_length)

my_breaks <- c(seq(min(mat), 0,
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(mat)/palette_length,
                   max(mat),
                   length.out=floor(palette_length/2)))

colnames(mat) <- gsub("A","Path1",colnames(mat))
colnames(mat) <- gsub("B","Path2",colnames(mat))
g <-pheatmap(mat[,c(rev(paste0(rep("Path1",9),seq(0,0.8,0.1))),paste0(rep("Path2",9),seq(0,0.8,0.1)))]
         ,cluster_cols = F, color = my_color,breaks = my_breaks)
g <- pheatmap2edit(g)
