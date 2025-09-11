#############################
# edgeR: exploring the data #
#############################
#Install package on R3.6.1
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("edgeR")

#Load package on R3.6.1
library("edgeR")
library("dplyr")
library("tidyr")


#############################
####  EXPLORE ABUNDANCE  ####
#############################
#Import data
x1=read.table(snakemake@input[['matrix1']],row.names="Name.2", stringsAsFactors=FALSE, header=T)
x2=read.table(snakemake@input[['matrix2']],row.names="Name.2", stringsAsFactors=FALSE, header=T)
design = readTargets(snakemake@input[['design']], sep="\t")

y1 = DGEList(counts=x1,group=design$stage)
y2 = DGEList(counts=x2,group=design$stage)

#Normalizing
y1 = calcNormFactors(y1, method = "TMM")
y2 = calcNormFactors(y2, method = "TMM")

#Estimating the dispersion
y1 = estimateCommonDisp(y1, verbose=TRUE)
y1 = estimateTagwiseDisp(y1)

y2 = estimateCommonDisp(y2, verbose=TRUE)
y2 = estimateTagwiseDisp(y2)

#MDS plot of the TOP 15000 genes
pdf(snakemake@output[["sirna"]])
plotMDS(y1, top = 1000, dim=c(1,2), col=c(rep("black",9), rep("pink",9), rep("blue",9), rep("purple",9), rep("red",9), rep("green",9)), pch=c(1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3), cex= 1.5, cex.axis= 1.0, main= list("MDS Plot: sRNA in DCL5 durum wheat genotypes", cex = 1.25), xlim=c(-5,12), ylim=c(-5,12),  xlab="logFC dim1", ylab="logFC dim2")
legend("topleft", legend=c("(AABB) Pre-meiotic", "(AABB) Mid meiotic", "(AABB) Post-meiotic", "(aabb) Pre-meiotic", "(aabb) Mid meiotic", "(aabb) Post-meiotic", "(aAbb) Pre-meiotic", "(aAbb) Mid meiotic", "(aAbb) Post-meiotic", "(aabB) Pre-meiotic", "(aabB) Mid meiotic", "(aabB) Post-meiotic", "(aabb-18C) Pre-meiotic", "(aabb-18C) Mid meiotic", "(aabb-18C) Post-meiotic", "(aabb-22C) Pre-meiotic", "(aabb-22C) Mid meiotic", "(aabb-22C) Post-meiotic"), col=c("black", "black", "black", "pink", "pink", "pink", "blue", "blue", "blue", "purple", "purple", "purple", "red", "red", "red", "green", "green", "green"), pch=c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3),box.lty=0, cex= 0.75)
dev.off()

pdf(snakemake@output[["ambiguous"]])
plotMDS(y2, top = 1000, dim=c(1,2), col=c(rep("black",9), rep("pink",9), rep("blue",9), rep("purple",9), rep("red",9), rep("green",9)), pch=c(1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3), cex= 1.5, cex.axis= 1.0, main= list("MDS Plot: sRNA in DCL5 durum wheat genotypes", cex = 1.25), xlim=c(-5,12), ylim=c(-5,12),  xlab="logFC dim1", ylab="logFC dim2")
legend("topleft", legend=c("(AABB) Pre-meiotic", "(AABB) Mid meiotic", "(AABB) Post-meiotic", "(aabb) Pre-meiotic", "(aabb) Mid meiotic", "(aabb) Post-meiotic", "(aAbb) Pre-meiotic", "(aAbb) Mid meiotic", "(aAbb) Post-meiotic", "(aabB) Pre-meiotic", "(aabB) Mid meiotic", "(aabB) Post-meiotic", "(aabb-18C) Pre-meiotic", "(aabb-18C) Mid meiotic", "(aabb-18C) Post-meiotic", "(aabb-22C) Pre-meiotic", "(aabb-22C) Mid meiotic", "(aabb-22C) Post-meiotic"), col=c("black", "black", "black", "pink", "pink", "pink", "blue", "blue", "blue", "purple", "purple", "purple", "red", "red", "red", "green", "green", "green"), pch=c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3),box.lty=0, cex= 0.75)
dev.off()
