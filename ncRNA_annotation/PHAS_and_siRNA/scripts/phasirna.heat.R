###########################################
# pheatmap: visualizing expression change #
###########################################
#Load R package
library("pheatmap")
library("RColorBrewer")


#Load data
x1 = read.delim(snakemake@input[['phas21']],row.names=1)
x2 = read.delim(snakemake@input[['phas24']],row.names=1)
x3 = read.delim(snakemake@input[['ambiguous21']],row.names=1)
x4 = read.delim(snakemake@input[['ambiguous24']],row.names=1)
x5 = read.delim(snakemake@input[['sirna']],row.names=1)

head(x1, n=5)
head(x2, n=5)
head(x3, n=5)
head(x4, n=5)
head(x5, n=5)

depth.NULL <- function(x, ...) 1

#Define metrics for clustering
drows1 <- "euclidean"
dcols1 <- "euclidean"

#Calculate and draw heatmaps
pdf(snakemake@output[["phas21"]],width=20,height=10)
phas21 <- list(
    x1,
    color                 = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(50),
    cellwidth             = 20,
    cellheight            = NA,
    scale                 = "row",
    treeheight_row        = 50,
    #breaks               = mat_breaks,
    #legend_breaks        = -1:10,
    kmeans_k              = NA,
    show_rownames         = F,
    show_colnames         = T,
    main                  = "21-nt phasiRNAs accumulating in anthers",
    clustering_method     = "ward.D2",
    cluster_rows          = FALSE,
    cluster_cols          = FALSE,
    clustering_distance_rows    = drows1,
    clustering_distance_cols    = dcols1)
do.call("pheatmap", phas21)
dev.off()


pdf(snakemake@output[["phas24"]],width=20,height=10)
phas24 <- list(
    x2,
    color                 = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(50),
    cellwidth             = 20,
    cellheight            = NA,
    scale                 = "row",
    treeheight_row        = 50,
    #breaks               = mat_breaks,
    #legend_breaks        = -1:10,
    kmeans_k              = NA,
    show_rownames         = F,
    show_colnames         = T,
    main                  = "24-nt phasiRNAs accumulating in anthers",
    clustering_method     = "ward.D2",
    cluster_rows          = FALSE,
    cluster_cols          = FALSE,
    clustering_distance_rows    = drows1,
    clustering_distance_cols    = dcols1)
do.call("pheatmap", phas24)
dev.off()


pdf(snakemake@output[["ambiguous21"]],width=20,height=10)
ambiguous21 <- list(
    x3,
    color                 = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(50),
    cellwidth             = 20,
    cellheight            = NA,
    scale                 = "row",
    treeheight_row        = 50,
    #breaks               = mat_breaks,
    #legend_breaks        = -1:10,
    kmeans_k              = NA,
    show_rownames         = F,
    show_colnames         = T,
    main                  = "Ambiguous 21-nt phasiRNAs accumulating in anthers",
    clustering_method     = "ward.D2",
    cluster_rows          = FALSE,
    cluster_cols          = FALSE,
    clustering_distance_rows    = drows1,
    clustering_distance_cols    = dcols1)
do.call("pheatmap", ambiguous21)
dev.off()


pdf(snakemake@output[["ambiguous24"]],width=20,height=10)
ambiguous24 <- list(
    x4,
    color                 = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(50),
    cellwidth             = 20,
    cellheight            = NA,
    scale                 = "row",
    treeheight_row        = 50,
    #breaks               = mat_breaks,
    #legend_breaks        = -1:10,
    kmeans_k              = NA,
    show_rownames         = F,
    show_colnames         = T,
    main                  = "Ambiguous 24-nt phasiRNAs accumulating in anthers",
    clustering_method     = "ward.D2",
    cluster_rows          = FALSE,
    cluster_cols          = FALSE,
    clustering_distance_rows    = drows1,
    clustering_distance_cols    = dcols1)
do.call("pheatmap", ambiguous24)
dev.off()


pdf(snakemake@output[["sirna"]],width=20,height=10)
sirna <- list(
    x5,
    color                 = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(50),
    cellwidth             = 20,
    cellheight            = NA,
    scale                 = "row",
    treeheight_row        = 50,
    #breaks               = mat_breaks,
    #legend_breaks        = -1:10,
    kmeans_k              = NA,
    show_rownames         = F,
    show_colnames         = T,
    main                  = "24-nt hc-siRNAs accumulating in anthers",
    clustering_method     = "ward.D2",
    cluster_rows          = FALSE,
    cluster_cols          = FALSE,
    clustering_distance_rows    = drows1,
    clustering_distance_cols    = dcols1)
do.call("pheatmap", sirna)
dev.off()
