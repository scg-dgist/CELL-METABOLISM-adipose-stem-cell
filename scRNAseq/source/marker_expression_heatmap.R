##### Code for
##### Author:  Eun Seo Park (evergreen@dgist.ac.kr)
###Load Data
library(Seurat)
library(pheatmap)
library(dplyr)
sapply(list.files("D:/OneDrive - dgist.ac.kr/Function/",full.names  = T),source)
###Load Data
library(dplyr)
library(Seurat)
NCD_seurat <- readRDS("D:/OneDrive - dgist.ac.kr/Adipocyte_Precursor/paper_final/NCD_Seurat.rds")
Idents(NCD_seurat) <- NCD_seurat$new_annot

#Find cluster marker
NCD_annot_marker <- FindAllMarkers(NCD_seurat)
NCD_annot_marker_sig <- subset(NCD_annot_marker, NCD_annot_marker$p_val_adj<0.05 & abs(NCD_annot_marker$avg_logFC)>0.5)

#Scailing gene expression profile
NCD_seurat <- ScaleData(NCD_seurat,scale.max = 3,features = rownames(NCD_seurat))
NCD_annot_marker_sig$cluster <- as.character(NCD_annot_marker_sig$cluster)
NCD_annot_marker_sig <- NCD_annot_marker_sig[order(NCD_annot_marker_sig$cluster),]

NCD_cluster_levels <- unique(NCD_annot_marker_sig$cluster)

#Average expression per cluster
cluster_means <- AverageExpression(NCD_seurat,features = unique(NCD_annot_marker_sig$gene),
                                   slot = "data")
scaled_tmp <- cluster_means$RNA

library(RColorBrewer)
library(pheatmap)
palette_length = 100
my_color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(palette_length)

my_breaks <- c(seq(min(scaled_tmp), 0,
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(scaled_tmp)/palette_length,
                   max(scaled_tmp),
                   length.out=floor(palette_length/2)))
g <- pheatmap(scaled_tmp[,NCD_cluster_levels], color = my_color, breaks = my_breaks,cluster_rows = F,cluster_cols = F)
g <- pheatmap2edit(g)


NCD_seurat <- readRDS("D:/OneDrive - dgist.ac.kr/Adipocyte_Precursor/paper_final/NCD_Seurat.rds")
annot_seurat <- NCD_seurat
Idents(annot_seurat) <- factor(annot_seurat$annot, levels = c("E1","E2","E3","E4",
                                                              "I1","I2","I3","I4",
                                                              "I5","I6","I7"))
annot_seurat <- ScaleData(annot_seurat, features = rownames(annot_seurat))
public_marker <- read.csv("I:/SNU_adipo_infla/public_scRNA_marker.csv")

names_PD <- as.vector(unique(gsub("_.*","",names(public_marker))))
for(k in 1:length(names_PD)){
  marker_df <- public_marker[,grep(names_PD[k],names(public_marker))]
  marker_gene <- levels(unlist(public_marker[,grep(names_PD[k],names(public_marker))]))
  marker_gene <- marker_gene[marker_gene %in% rownames(NCD_seurat)]
  
  mean_scaled_exp <-AverageExpression(annot_seurat,features = marker_gene,
                                      slot = "scale.data") 
  gene_df <- data.frame(row.names =marker_gene)
  for( i in 1:nrow(gene_df)){
    for( j in 1:ncol(marker_df)){
      if(rownames(gene_df)[i] %in% marker_df[,colnames(marker_df)[j]]){
        gene_df[i,"label"] <- colnames(marker_df)[j]
      }
    }
  }
  
  cols <- gg_color_hue(ncol(public_marker))
  names(cols) <- colnames(public_marker)
  
  
  
  
  scaled_tmp <-mean_scaled_exp$RNA
  if(nrow(scaled_tmp)!=55){
    mat_tmp <- matrix(data = 0,nrow = 55-nrow(scaled_tmp),ncol = nlevels(Idents(annot_seurat)))
    rownames(mat_tmp) <- 1:nrow(mat_tmp)
    colnames(mat_tmp) <- levels(Idents(annot_seurat))
    final_tmp <- rbind(scaled_tmp,mat_tmp)
  }else{
    final_tmp <-scaled_tmp
  }
  
  library(RColorBrewer)
  palette_length = 100
  my_color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(palette_length)
  my_breaks <- c(seq(min(scaled_tmp), 0,
                     length.out=ceiling(palette_length/2) + 1),
                 seq(max(scaled_tmp)/palette_length,
                     max(scaled_tmp),
                     length.out=floor(palette_length/2)))
  
  g<- pheatmap::pheatmap(final_tmp,color = my_color,breaks = my_breaks,cluster_cols = F,
                         cluster_rows = F,height = 40,
                         width = 40,fontsize = 5,
                         annotation_row = gene_df,annotation_colors = list(label=cols[colnames(marker_df)]),
                         border_color = NA,na_col = "grey")
  g <- pheatmap2edit(g)
  g <- create_pptx(g,"supple_public_heatmap.pptx")
}

for(k in 1:length(names_PD)){
  marker_df <- public_marker[,grep(names_PD[k],names(public_marker))]
  marker_gene <- levels(unlist(public_marker[,grep(names_PD[k],names(public_marker))]))
  marker_gene <- marker_gene[marker_gene %in% rownames(NCD_seurat)]
  
  mean_scaled_exp <-AverageExpression(annot_seurat,features = marker_gene,
                                      slot = "data") 
  gene_df <- data.frame(row.names =marker_gene)
  for( i in 1:nrow(gene_df)){
    for( j in 1:ncol(marker_df)){
      if(rownames(gene_df)[i] %in% marker_df[,colnames(marker_df)[j]]){
        gene_df[i,"label"] <- colnames(marker_df)[j]
      }
    }
  }
  
  cols <- gg_color_hue(ncol(public_marker))
  names(cols) <- colnames(public_marker)
  
  
  
  
  scaled_tmp <-mean_scaled_exp$RNA
  if(nrow(scaled_tmp)!=55){
    mat_tmp <- matrix(data = 0,nrow = 55-nrow(scaled_tmp),ncol = nlevels(Idents(annot_seurat)))
    rownames(mat_tmp) <- 1:nrow(mat_tmp)
    colnames(mat_tmp) <- levels(Idents(annot_seurat))
    final_tmp <- rbind(scaled_tmp,mat_tmp)
  }else{
    final_tmp <-scaled_tmp
  }
  
  library(RColorBrewer)
  palette_length = 100
  my_color = colorRampPalette(brewer.pal(9, "Reds"))(palette_length)

  g<- pheatmap::pheatmap(final_tmp,color = my_color,cluster_cols = F,
                         cluster_rows = F,height = 40,
                         width = 40,fontsize = 5,
                         annotation_row = gene_df,annotation_colors = list(label=cols[colnames(marker_df)]),
                         border_color = NA,na_col = "grey")
  g <- pheatmap2edit(g)
  g <- create_pptx(g,"supple_public_heatmap.pptx")
}



