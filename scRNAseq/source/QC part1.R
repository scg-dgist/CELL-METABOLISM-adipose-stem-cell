##### Code for
##### Author:  Eun Seo Park (evergreen@dgist.ac.kr)
library(DropletUtils)
setwd("D:/espark/Adipocyte_precursor/A1")
dir.A1 <- "D:/espark/Adipocyte_precursor/A1/raw_feature_bc_matrix"
list.files(dir.A1)
set.seed(0)
A1 <- read10xCounts(dir.A1)
class(counts(A1))

br.out <- barcodeRanks(counts(A1))
sort(A1$total_counts,decreasing = TRUE)
# # Making a plot.
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")
abline(h=br.out$knee, col="dodgerblue", lty=2)
abline(h=br.out$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
       legend=c("knee", "inflection"))

set.seed(100)
e.out.A1 <- emptyDrops(counts(A1))
is.cell.A1 <- e.out.A1$FDR <= 0.01
sum(is.cell.A1, na.rm=TRUE)
colnames(A1) = colData(A1)$Barcode

cd.A1 = counts(A1)[,which(e.out.A1$FDR <= 0.01)]
# write.csv(as.matrix(cd.A1), file="counts.csv")

hist(A1$total_counts, 
     breaks=100, col="grey80",
     xlab="Log-total UMI count")

hist(H9$total_counts, 
     breaks=100, col="grey80",
     xlab="Log-total UMI count")

A1 <- A1[,which(e.out.A1$FDR < 0.01)]

library(scater)
library(EnsDb.Mmusculus.v79)
location <- mapIds(EnsDb.Mmusculus.v79, keys=rowData(A1)$ID, 
                   column="SEQNAME", keytype="GENEID")

rowData(A1)$CHR <- location
summary(location=="MT")
A1 <- calculateQCMetrics(A1, feature_controls=list(Mito=which(location=="MT")))

# Dimensionality reduction plots
A1 <- runPCA(A1, use_coldata=TRUE)
# reducedDimNames(A1)
plotReducedDim(A1, use_dimred = "PCA_coldata")


# par(mfrow=c(1,3))

dev.off()

hist(A1$log10_total_counts, 
     breaks=100, col="grey80",
     xlab="Log-total UMI count")

hist(A1$log10_total_features_by_counts, 
     breaks=100, col="grey80",
     xlab="Log-total number of expressed features")

hist(A1$pct_counts_Mito, 
     breaks=100, col="grey80",
     xlab="Proportion of reads in mitochondrial genes")

filtering <- numeric()
filtering[which(A1$log10_total_counts > 3.0 & A1$pct_counts_Mito <= 3 & A1$log10_total_features_by_counts > 3)] <- 1
filtering[which(A1$log10_total_counts <= 3.0 | A1$pct_counts_Mito > 3 | A1$log10_total_features_by_counts <= 3)] <- 0
A1$filtering <- filtering
table(A1$filtering)

ggplot(as.data.frame(A1@reducedDims$PCA_coldata), 
       aes(x=PC1, y=PC2, color = as.factor(A1$filtering))) +
  geom_point() + 
  theme_bw() + 
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

ggplot(as.data.frame(A1@reducedDims$PCA_coldata),
       aes(x=PC1, y=PC2, color = A1$log10_total_counts)) +
  geom_point() + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())




####A2
library(DropletUtils)
setwd("D:/espark/Adipocyte_precursor/A2")
dir.A2 <- "D:/espark/Adipocyte_precursor/A2/raw_feature_bc_matrix"
list.files(dir.A2)
set.seed(0)
A2 <- read10xCounts(dir.A2)
class(counts(A2))


br.out <- barcodeRanks(counts(A2))

# # Making a plot.
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")
abline(h=br.out$knee, col="dodgerblue", lty=2)
abline(h=br.out$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
       legend=c("knee", "inflection"))

set.seed(100)
e.out.A2 <- emptyDrops(counts(A2))
is.cell.A2 <- e.out.A2$FDR <= 0.01
sum(is.cell.A2, na.rm=TRUE)
colnames(A2) = colData(A2)$Barcode
cd.A2 = counts(A2)[,which(e.out.A2$FDR <= 0.01)]

A2 <- A2[,which(e.out.A2$FDR < 0.01)]
library(scater)
library(EnsDb.Mmusculus.v79)

location <- mapIds(EnsDb.Mmusculus.v79, keys=rowData(A2)$ID, 
                   column="SEQNAME", keytype="GENEID")

rowData(A2)$CHR <- location
summary(location=="MT")
A2 <- calculateQCMetrics(A2, feature_controls=list(Mito=which(location=="MT")))

# Dimensionality reduction plots
A2 <- runPCA(A2, use_coldata=TRUE)
# reducedDimNames(A2)
plotReducedDim(A2, use_dimred = "PCA_coldata")


# par(mfrow=c(1,3))
dev.off()

hist(A2$log10_total_counts, 
     breaks=100, col="grey80",
     xlab="Log-total UMI count")

hist(A2$log10_total_features_by_counts, 
     breaks=100, col="grey80",
     xlab="Log-total number of expressed features")

hist(A2$pct_counts_Mito, 
     breaks=100, col="grey80",
     xlab="Proportion of reads in mitochondrial genes")

filtering <- numeric()
filtering[which(A2$log10_total_counts > 3& A2$pct_counts_Mito <= 5 & A2$log10_total_features_by_counts > 3)] <- 1
filtering[which(A2$log10_total_counts <= 3 | A2$pct_counts_Mito > 5 | A2$log10_total_features_by_counts <= 3)] <- 0
A2$filtering <- filtering
table(A2$filtering)

ggplot(as.data.frame(A2@reducedDims$PCA_coldata), 
       aes(x=PC1, y=PC2, color = as.factor(A2$filtering))) +
  geom_point() + 
  theme_bw() + 
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())


ggplot(as.data.frame(A2@reducedDims$PCA_coldata),
       aes(x=PC1, y=PC2, color = A2$log10_total_features_by_counts)) +
  geom_point() + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())


##H9

library(DropletUtils)
setwd("D:/espark/Adipocyte_precursor/H9")
dir.H9 <- "D:/espark/Adipocyte_precursor/H9/raw_feature_bc_matrix"
list.files(dir.H9)
set.seed(0)
H9 <- read10xCounts(dir.H9)
class(counts(H9))


br.out <- barcodeRanks(counts(H9))

# # Making a plot.
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")
abline(h=br.out$knee, col="dodgerblue", lty=2)
abline(h=br.out$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
       legend=c("knee", "inflection"))

set.seed(100)
e.out.H9 <- emptyDrops(counts(H9))
is.cell.H9 <- e.out.H9$FDR <= 0.01
sum(is.cell.H9, na.rm=TRUE)
colnames(H9) = colData(H9)$Barcode
cd.H9 = counts(H9)[,which(e.out.H9$FDR <= 0.01)]

H9 <- H9[,which(e.out.H9$FDR < 0.01)]
library(scater)
library(EnsDb.Mmusculus.v79)

location <- mapIds(EnsDb.Mmusculus.v79, keys=rowData(H9)$ID, 
                   column="SEQNAME", keytype="GENEID")

rowData(H9)$CHR <- location
summary(location=="MT")
H9 <- calculateQCMetrics(H9, feature_controls=list(Mito=which(location=="MT")))

# Dimensionality reduction plots
H9 <- runPCA(H9, use_coldata=TRUE)
# reducedDimNames(H9)
plotReducedDim(H9, use_dimred = "PCA_coldata")


# par(mfrow=c(1,3))
dev.off()

hist(H9$log10_total_counts, 
     breaks=100, col="grey80",
     xlab="Log-total UMI count")

hist(H9$log10_total_features_by_counts, 
     breaks=100, col="grey80",
     xlab="Log-total number of expressed features")

hist(H9$pct_counts_Mito, 
     breaks=100, col="grey80",
     xlab="Proportion of reads in mitochondrial genes")

filtering <- numeric()
filtering[which(H9$log10_total_counts > 3 & H9$pct_counts_Mito <= 5 & H9$log10_total_features_by_counts > 3)] <- 1
filtering[which(H9$log10_total_counts <= 3 | H9$pct_counts_Mito > 5 | H9$log10_total_features_by_counts <= 3)] <- 0
H9$filtering <- filtering
table(H9$filtering)

ggplot(as.data.frame(H9@reducedDims$PCA_coldata), 
       aes(x=PC1, y=PC2, color = as.factor(H9$filtering))) +
  geom_point() + 
  theme_bw() + 
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())




ggplot(as.data.frame(H9@reducedDims$PCA_coldata),
       aes(x=PC1, y=PC2, color = H9$pct_counts_Mito)) +
  geom_point() + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

library(DropletUtils)
setwd("D:/espark/Adipocyte_precursor/H10")
dir.H10 <- "D:/espark/Adipocyte_precursor/H10/raw_feature_bc_matrix"
list.files(dir.H10)
set.seed(0)
H10 <- read10xCounts(dir.H10)
class(counts(H10))

br.out <- barcodeRanks(counts(H10))
sort(H10$total_counts,decreasing = TRUE)
# # Making a plot.
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")
abline(h=br.out$knee, col="dodgerblue", lty=2)
abline(h=br.out$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
       legend=c("knee", "inflection"))

set.seed(100)
e.out.H10 <- emptyDrops(counts(H10))
is.cell.H10 <- e.out.H10$FDR <= 0.01
sum(is.cell.H10, na.rm=TRUE)
colnames(H10) = colData(H10)$Barcode

cd.H10 = counts(H10)[,which(e.out.H10$FDR <= 0.01)]
# write.csv(as.matrix(cd.H10), file="counts.csv")

hist(H10$total_counts, 
     breaks=100, col="grey80",
     xlab="Log-total UMI count")

hist(H9$total_counts, 
     breaks=100, col="grey80",
     xlab="Log-total UMI count")

H10 <- H10[,which(e.out.H10$FDR < 0.01)]

library(scater)
library(EnsDb.Mmusculus.v79)
location <- mapIds(EnsDb.Mmusculus.v79, keys=rowData(H10)$ID, 
                   column="SEQNAME", keytype="GENEID")

rowData(H10)$CHR <- location
summary(location=="MT")
H10 <- calculateQCMetrics(H10, feature_controls=list(Mito=which(location=="MT")))

# Dimensionality reduction plots
H10 <- runPCA(H10, use_coldata=TRUE)
# reducedDimNames(H10)
plotReducedDim(H10, use_dimred = "PCA_coldata")


# par(mfrow=c(1,3))

dev.off()

hist(H10$log10_total_counts, 
     breaks=100, col="grey80",
     xlab="Log-total UMI count")

hist(H10$log10_total_features_by_counts, 
     breaks=100, col="grey80",
     xlab="Log-total number of expressed features")

hist(H10$pct_counts_Mito, 
     breaks=100, col="grey80",
     xlab="Proportion of reads in mitochondrial genes")

filtering <- numeric()
filtering[which(H10$log10_total_counts > 3.0 & H10$pct_counts_Mito <= 5 & H10$log10_total_features_by_counts > 3)] <- 1
filtering[which(H10$log10_total_counts <= 3.0 | H10$pct_counts_Mito > 5 | H10$log10_total_features_by_counts <= 3)] <- 0
H10$filtering <- filtering
table(H10$filtering)

ggplot(as.data.frame(H10@reducedDims$PCA_coldata), 
       aes(x=PC1, y=PC2, color = as.factor(H10$filtering))) +
  geom_point() + 
  theme_bw() + 
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

ggplot(as.data.frame(H10@reducedDims$PCA_coldata),
       aes(x=PC1, y=PC2, color = H10$log10_total_counts)) +
  geom_point() + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

###########H11
######
library(DropletUtils)
setwd("D:/espark/Adipocyte_precursor/H11")
dir.H11 <- "D:/espark/Adipocyte_precursor/H11/raw_feature_bc_matrix"
list.files(dir.H11)
set.seed(0)
H11 <- read10xCounts(dir.H11)
class(counts(H11))

br.out <- barcodeRanks(counts(H11))
sort(H11$total_counts,decreasing = TRUE)
# # Making a plot.
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")
abline(h=br.out$knee, col="dodgerblue", lty=2)
abline(h=br.out$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
       legend=c("knee", "inflection"))

set.seed(100)
e.out.H11 <- emptyDrops(counts(H11))
is.cell.H11 <- e.out.H11$FDR <= 0.01
sum(is.cell.H11, na.rm=TRUE)
colnames(H11) = colData(H11)$Barcode

cd.H11 = counts(H11)[,which(e.out.H11$FDR <= 0.01)]
# write.csv(as.matrix(cd.H11), file="counts.csv")

hist(H11$total_counts, 
     breaks=100, col="grey80",
     xlab="Log-total UMI count")

hist(H9$total_counts, 
     breaks=100, col="grey80",
     xlab="Log-total UMI count")

H11 <- H11[,which(e.out.H11$FDR < 0.01)]

library(scater)
library(EnsDb.Mmusculus.v79)
location <- mapIds(EnsDb.Mmusculus.v79, keys=rowData(H11)$ID, 
                   column="SEQNAME", keytype="GENEID")

rowData(H11)$CHR <- location
summary(location=="MT")
H11 <- calculateQCMetrics(H11, feature_controls=list(Mito=which(location=="MT")))

# Dimensionality reduction plots
H11 <- runPCA(H11, use_coldata=TRUE)
# reducedDimNames(H11)
plotReducedDim(H11, use_dimred = "PCA_coldata")


# par(mfrow=c(1,3))

dev.off()

hist(H11$log10_total_counts, 
     breaks=100, col="grey80",
     xlab="Log-total UMI count")

hist(H11$log10_total_features_by_counts, 
     breaks=100, col="grey80",
     xlab="Log-total number of expressed features")

hist(H11$pct_counts_Mito, 
     breaks=100, col="grey80",
     xlab="Proportion of reads in mitochondrial genes")

filtering <- numeric()
filtering[which(H11$log10_total_counts > 3.0 & H11$pct_counts_Mito <= 3 & H11$log10_total_features_by_counts > 3)] <- 1
filtering[which(H11$log10_total_counts <= 3.0 | H11$pct_counts_Mito > 3 | H11$log10_total_features_by_counts <= 3)] <- 0
H11$filtering <- filtering
table(H11$filtering)

ggplot(as.data.frame(H11@reducedDims$PCA_coldata), 
       aes(x=PC1, y=PC2, color = as.factor(H11$filtering))) +
  geom_point() + 
  theme_bw() + 
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

ggplot(as.data.frame(H11@reducedDims$PCA_coldata),
       aes(x=PC1, y=PC2, color = H11$log10_total_counts)) +
  geom_point() + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

######H12
#######

library(DropletUtils)
setwd("D:/espark/Adipocyte_precursor/H12")
dir.H12 <- "D:/espark/Adipocyte_precursor/H12/raw_feature_bc_matrix"
list.files(dir.H12)
set.seed(0)
H12 <- read10xCounts(dir.H12)
class(counts(H12))

br.out <- barcodeRanks(counts(H12))
sort(H12$total_counts,decreasing = TRUE)
# # Making a plot.
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")
abline(h=br.out$knee, col="dodgerblue", lty=2)
abline(h=br.out$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
       legend=c("knee", "inflection"))

set.seed(100)
e.out.H12 <- emptyDrops(counts(H12))
is.cell.H12 <- e.out.H12$FDR <= 0.01
sum(is.cell.H12, na.rm=TRUE)
colnames(H12) = colData(H12)$Barcode

cd.H12 = counts(H12)[,which(e.out.H12$FDR <= 0.01)]
# write.csv(as.matrix(cd.H12), file="counts.csv")

hist(H12$total_counts, 
     breaks=100, col="grey80",
     xlab="Log-total UMI count")

hist(H9$total_counts, 
     breaks=100, col="grey80",
     xlab="Log-total UMI count")

H12 <- H12[,which(e.out.H12$FDR < 0.01)]

library(scater)
library(EnsDb.Mmusculus.v79)
location <- mapIds(EnsDb.Mmusculus.v79, keys=rowData(H12)$ID, 
                   column="SEQNAME", keytype="GENEID")

rowData(H12)$CHR <- location
summary(location=="MT")
H12 <- calculateQCMetrics(H12, feature_controls=list(Mito=which(location=="MT")))

# Dimensionality reduction plots
H12 <- runPCA(H12, use_coldata=TRUE)
# reducedDimNames(H12)
plotReducedDim(H12, use_dimred = "PCA_coldata")


# par(mfrow=c(1,3))

dev.off()

hist(H12$log10_total_counts, 
     breaks=100, col="grey80",
     xlab="Log-total UMI count")

hist(H12$log10_total_features_by_counts, 
     breaks=100, col="grey80",
     xlab="Log-total number of expressed features")

hist(H12$pct_counts_Mito, 
     breaks=100, col="grey80",
     xlab="Proportion of reads in mitochondrial genes")

filtering <- numeric()
filtering[which(H12$log10_total_counts > 3.0 & H12$pct_counts_Mito <= 5 & H12$log10_total_features_by_counts > 3)] <- 1
filtering[which(H12$log10_total_counts <= 3.0 | H12$pct_counts_Mito > 5 | H12$log10_total_features_by_counts <= 3)] <- 0
H12$filtering <- filtering
table(H12$filtering)

ggplot(as.data.frame(H12@reducedDims$PCA_coldata), 
       aes(x=PC1, y=PC2, color = as.factor(H12$filtering))) +
  geom_point() + 
  theme_bw() + 
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

ggplot(as.data.frame(H12@reducedDims$PCA_coldata),
       aes(x=PC1, y=PC2, color = H12$log10_total_counts)) +
  geom_point() + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

saveRDS(A1, "D:/espark/Adipocyte_precursor/dropletfiltered_mitox/dropletfilterA1set.rds")
saveRDS(A2, "D:/espark/Adipocyte_precursor/dropletfiltered_mitox/dropletfilterA2set.rds")
saveRDS(H9, "D:/espark/Adipocyte_precursor/dropletfiltered_mitox/dropletfilterH9set.rds")
saveRDS(H10, "D:/espark/Adipocyte_precursor/dropletfiltered_mitox/dropletfilterH10set.rds")
saveRDS(H11, "D:/espark/Adipocyte_precursor/dropletfiltered_mitox/dropletfilterH11set.rds")
saveRDS(H12, "D:/espark/Adipocyte_precursor/dropletfiltered_mitox/dropletfilterH12set.rds")

