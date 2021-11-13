##### Code for
##### Author:  Eun Seo Park (evergreen@dgist.ac.kr)
###Load Data
library(Seurat)
EAT_seurat <- readRDS("EAT_seurat.rds")
Idents(EAT_seurat) <- EAT_seurat$seurat_clusters
E10_marker <- FindMarkers(EAT_seurat, ident.1 = 3)
up_sig_E10_marker <- subset(E10_marker,E10_marker$avg_logFC>0&E10_marker$p_val_adj<0.05)
down_sig_E10_marker <- subset(E10_marker,E10_marker$avg_logFC<0&E10_marker$p_val_adj<0.05)
sig_E10_marker <- subset(E10_marker, E10_marker$p_val_adj<0.05)

library(topGO)
library(plyr)
topGO_function <- function(targetGenes, backgroundGenes, onts=c("BP"), topnodes=1000){
  
  allGene = factor(as.integer(backgroundGenes %in% targetGenes))
  names(allGene) = backgroundGenes
  
  tab = as.list(onts)
  names(tab) = onts
  for(i in 1:length(onts)){
    tgd = new("topGOdata", ontology=onts[i], allGenes=allGene, nodeSize=3,
              annot=annFUN.org, mapping="org.Mm.eg.db", ID="Symbol")
    resultTopGO.elim = runTest(tgd, algorithm="elim", statistic="Fisher")
    tab[[i]] = GenTable(tgd, Fisher.elim=resultTopGO.elim,
                        orderBy="Fisher.classic", topNodes = topnodes,
                        numChar=1000L)
    tab[[i]]$onts = onts[i]
  }
  topGOResults = rbind.fill(tab)
  topGOResults$Fisher.elim = as.numeric(topGOResults$Fisher.elim)
  topGOResults$Fisher.elim[is.na(topGOResults$Fisher.elim)] = 0
  #topGOResults = subset(topGOResults, Fisher.elim < padjval)
  topGOResults$log10P = -log10(topGOResults$Fisher.elim)
  return(topGOResults)
}
UP_GO_E10 <- topGO_function(rownames(up_sig_E10_marker), rownames(EAT_seurat),
                           onts = c("MF","BP","CC"))
Down_GO_E10 <- topGO_function(rownames(down_sig_E10_marker), rownames(EAT_seurat),
                           onts = c("MF","BP","CC"))

