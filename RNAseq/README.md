# Transcriptome analyses





## Software version

```
fastp v0.24.0
salmon v1.10.3
r v4.3.2
edgeR v4.0.16
```

---
### 1. Quality Control

**📥 Inputs**  
• `*.fastq`: raw RNA-seq data from NCBI

**📥 Outputs**  
• `*.filtered.fastq`: filtered, trimmed reads

⚙️ **Trim with fastp**  
```
#for paired-end 
for fq1 in *_1.fastq; do
  fq2=$(echo $fq1  | sed 's/_1.fastq/_2.fastq/g')
  out1=$(echo $fq1 | sed 's/fastq/filtered.fastq/g')
  out2=$(echo $fq2 | sed 's/fastq/filtered.fastq/g')
  fastp --in1 $fq1 --out1 $out1 --in2 $fq2 --out2 $out2 -q 20 --length_required 50 --detect_adapter_for_pe -w 16
done

#for single end
for fq1 in *.fastq; do
  out1=$(echo $fq1 | sed 's/fastq/filtered.fastq/g')
  fastp --in1 $fq1 --out1 $out1 -q 20 --length_required 50 -w 16
done
```
---
### 2. Quantification
**📥 Inputs**  
• `*.filtered.fastq`: filtered, trimmed reads  
• `Kronos.collapsed.chromosomes.masked.v1.1.fa`: reference genome  
• `Kronos.v2.1.gff3`: v2.1 annotations  

**📥 Outputs**  
• `quant.genes.sf`: Gene-level quantifications    

⚙️ **Indexing**  
```
#prepare database
grep ">" Kronos.collapsed.chromosomes.masked.v1.1.fa | cut -d " " -f 1 | cut -d ">" -f 2 > decoys.txt
gffread -w Kronos.v2.1.transcripts.fa -g Kronos.collapsed.chromosomes.masked.v1.1.fa Kronos.v2.1.gff3
cat Kronos.v2.1.transcripts.fa Kronos.collapsed.chromosomes.masked.v1.1.fa > Kronos.gentrome.fa

#index gentrome
salmon index -t Kronos.gentrome.fa -d decoys.txt -p 30 -i salmon_index
```
⚙️ **Quantification**  
```
#for paired-end data:
indir=$1
for fq1 in ${indir}/*_1.filtered.fastq; do
  prefix=$(basename "$fq1" | cut -d "_" -f 1) 
  fq2=$(echo $fq | sed 's/_1.filtered/_2.filtered/g')

  salmon quant -l A -1 "$fq1" -2 "$fq2" -p 40 \
    -g Kronos.v2.1.gtf \
    -i salmon_index/ \
    -o "${prefix}" --validateMappings
done
```
```
#for single_end
indir=$1

for fq1 in ${indir}/*.filtered.fastq; do
  prefix=$(basename "$fq1" | cut -d "." -f 1) 

  salmon quant -l A -r "$fq1" -p 40 \
    -g Kronos.v2.1.gtf \
    -i salmon_index/ \
    -o "${prefix}" --validateMappings
done
```
----

### 3. Counts to percentiles
**📥 Inputs**  
• `quant.genes.sf`: Gene-level quantifications  

**📥 Outputs**  
• `NLR_GlobalrankPercentile`: Percentiles for NLR genes


```
#!/usr/bin/env Rscript
library(tximport)
library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(stringr)
library(edgeR)
library(purrr)
library(ggplot2)
library(matrixStats)
library(ggrepel)


# ----------- set input names -------------
folder_name = "PRJNA673987"  # <------------- a folder that contains all SRR* within a single bioproject
nlr_file    = "nlr_gene_ids.txt"
sample_tbl  = "samples.tsv" # includes SRR accessions and their conditions

# ----------- set output names -------------
raw_counts  = paste0(folder_name, "_salmon.rawcounts.txt")

# ----------- set opath -------------
setwd(paste0('./', folder_name))

# ----------- Load sample metadata -----------
samples <- read.table("samples.tsv", header= T, sep = "\t") 

##########################################
# ----------------  PCA  -----------------
##########################################
# ----------- Read quant.genes.sf counts ----------
read_salmon_counts <- function(fpath) {
  df <- read_tsv(file.path(fpath, "quant.genes.sf"), show_col_types = FALSE)
  df <- df %>% select(gene_id = Name, NumReads)
  return(df)
}

cts_list <- lapply(samples$path, read_salmon_counts)
names(cts_list) <- samples$sample_id

# Merge into count matrix
counts_df <- Reduce(function(x, y) full_join(x, y, by="gene_id"), cts_list)
colnames(counts_df)[-1] <- samples$sample_id
counts_df[is.na(counts_df)] <- 0
cts_mat <- as.matrix(counts_df[,-1])
rownames(cts_mat) <- counts_df$gene_id
storage.mode(cts_mat) <- "numeric"

## ---------- export raw count table ---------------------
write_tsv(counts_df, file = raw_counts)

# ----------- edgeR: DGEList + TMM normalization ----------
group <- factor(samples$condition) 
dge <- DGEList(counts = cts_mat, group = group)
dge <- calcNormFactors(dge, method = "TMM")
logCPM <- cpm(dge, log=TRUE, prior.count=1)

# ----------- PCA on logCPM -----------
logCPM_filtered <- logCPM[apply(logCPM, 1, function(x) var(x) > 0), ]
pca <- prcomp(t(logCPM_filtered), scale. = TRUE)
pca_df <- as.data.frame(pca$x)
pca_df$sample_id <- rownames(pca_df)
pca_df <- left_join(pca_df, samples, by = "sample_id")

# ----------- Plot PCA -----------
# Compute variance explained
explained <- pca$sdev^2 / sum(pca$sdev^2)

# Plot PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = as.factor(condition))) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = sample_id), size = 3, max.overlaps = 50) +
  theme_bw() +
  labs(
    x = sprintf("PC1: %.1f%%", 100 * explained[1]),
    y = sprintf("PC2: %.1f%%", 100 * explained[2])
  ) +
  theme(plot.title = element_text(hjust = 0.5))

# -----------  PCA without outliers -----------
outlier_ids <- c("") # <------------- remove outlier samples
keep_ids <- setdiff(colnames(logCPM), outlier_ids)   # your chosen list
scores <- pca_df[pca_df$sample_id %in% keep_ids, ]
ggplot(scores, aes(PC1, PC2, color = as.factor(condition))) +
  geom_point(size=3) + ggrepel::geom_text_repel(aes(label=sample_id), size=3) +
  theme_bw() +
  labs(title="PCA (frozen basis, outliers hidden)",
       x=sprintf("PC1: %.1f%%", 100* (pca_df$sdev[1]^2/sum(pca_df$sdev^2))),
       y=sprintf("PC2: %.1f%%", 100* (pca_df$sdev[2]^2/sum(pca_df$sdev^2))))


###################################################################
############ export normalized NLR percentiles#####################
###################################################################
# ----------- set input  -------------
outlier_ids <- c("") # <------------- remove outlier samples
samples_clean <- samples[!samples$sample_id %in% outlier_ids, , drop=FALSE]
stopifnot(all(file.exists(file.path(samples_clean$path, "quant.genes.sf"))))

files <- setNames(file.path(samples_clean$path, "quant.genes.sf"),
                  samples_clean$sample_id)

# NLR list
nlr_ids <- unique(read.csv(nlr_file, header = FALSE)[,1])

# Group labels
samples_clean$grp <- paste(samples_clean$project, samples_clean$condition, sep="|")

# ----------- Normalization  -------------
txi <- tximport(files,
                type = "salmon",
                txOut = TRUE,                        # rows = entries from quant.genes.sf (genes)
                countsFromAbundance = "lengthScaledTPM",
                ignoreTxVersion = TRUE)

dge <- DGEList(counts = txi$counts)

# Keep lowly-expressed genes per edgeR + ALWAYS keep NLRs
keep_base <- filterByExpr(dge, group = samples_clean$condition)
keep_nlr  <- rownames(dge) %in% nlr_ids
keep      <- keep_base | keep_nlr

dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

# Normalized counts
cpm_mat    <- cpm(dge, log = FALSE)                     # TMM-normalized CPM
logCPM_mat <- cpm(dge, log = TRUE, prior.count = 1)

# ---------------- Subset to NLRs--------------------------------------------------
nlr_present <- intersect(nlr_ids, rownames(cpm_mat))
nlr_missing <- setdiff(nlr_ids, rownames(cpm_mat))

nlr_cpm <- rbind(
  cpm_mat[nlr_present, , drop=FALSE],
  matrix(0, nrow = length(nlr_missing), ncol = ncol(cpm_mat),
         dimnames = list(nlr_missing, colnames(cpm_mat)))
)[nlr_ids, , drop=FALSE]

nlr_logCPM <- rbind(
  logCPM_mat[nlr_present, , drop=FALSE],
  matrix(0, nrow = length(nlr_missing), ncol = ncol(logCPM_mat),
         dimnames = list(nlr_missing, colnames(logCPM_mat)))
)[nlr_ids, , drop=FALSE]

#--------------- NLR ranking --------------------------------------
rank_mat <- apply(nlr_cpm, 2, function(x) rank(-x, ties.method = "average"))
nlr_rank_pct <- 1 - (rank_mat - 1) / (nrow(rank_mat) - 1)

# Optionally aggregate replicates: median by (project × condition)
split_idx <- split(seq_len(ncol(nlr_cpm)), samples_clean$grp)

rowMed_by_group <- function(mat, idx) {
  out <- do.call(cbind, lapply(idx, function(ix) matrixStats::rowMedians(mat[, ix, drop = FALSE])))
  colnames(out) <- names(idx); rownames(out) <- rownames(mat); out
}

nlr_rank_pct_grp <- rowMed_by_group(nlr_rank_pct, split_idx)     # for clustering / global comparisons
nlr_CPM_grp   <- sapply(split_idx, function(ix) rowMeans(nlr_cpm[, ix, drop=FALSE]))
nlr_logCPM_grp   <- sapply(split_idx, function(ix) rowMeans(nlr_logCPM[, ix, drop=FALSE]))

#--------------- global ranking --------------------------------------
rank_all <- apply(cpm_mat, 2, function(x) rank(-x, ties.method = "average"))
rank_pct_all <- 1 - (rank_all - 1) / (nrow(rank_all) - 1)   # 0..1, 1 = highest

## NLRs in *global* rank context
nlr_rank_pct_global <- rank_pct_all[intersect(nlr_ids, rownames(rank_pct_all)), , drop = FALSE]

## For consistency with your within-NLR object order, fill & re-order:
nlr_missing_global <- setdiff(nlr_ids, rownames(nlr_rank_pct_global))
if (length(nlr_missing_global)) {
  nlr_rank_pct_global <- rbind(
    nlr_rank_pct_global,
    matrix(NA_real_, nrow = length(nlr_missing_global), ncol = ncol(rank_pct_all),
           dimnames = list(nlr_missing_global, colnames(rank_pct_all)))
  )
}
nlr_rank_pct_global <- nlr_rank_pct_global[nlr_ids, , drop = FALSE]
nlr_rank_pct_global_grp <- rowMed_by_group(nlr_rank_pct_global, split_idx)


# ----------- Write output  -------------
groups_keep <- c("|Base", "|Middle", "|Top"   ) # <------------- Conditions to export
cols_keep <- samples_clean$grp %in% groups_keep
sample_ids_keep <- samples_clean$sample_id[cols_keep]

# Group-averaged
write_tsv(
  as_tibble(tibble::rownames_to_column(as.data.frame(nlr_rank_pct_global_grp[, groups_keep, drop=FALSE]),
                                       var = "gene_id")),
  "NLR_GlobalrankPercentile.avg_by_group.SELECTED.tsv"
)

write_tsv(
  as_tibble(tibble::rownames_to_column(as.data.frame(nlr_CPM_grp[nlr_ids, groups_keep, drop=FALSE]),
                                       var = "gene_id")),
  "NLR_TMM_CPM.avg_by_group.SELECTED.tsv"
)



###################################################################
#########     percentile vs CPM comparison    #####################
###################################################################
# thresholds
cutoffs <- c(0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90)

# function to get CPM values closest to each cutoff
get_cpm_at_percentile <- function(cpm_vec, rank_vec, cutoffs) {
  sapply(cutoffs, function(ct) {
    # pick the gene(s) closest to the cutoff
    idx <- which.min(abs(rank_vec - ct))
    cpm_vec[idx]
  })
}

# apply across all samples
cpm_cutoffs <- sapply(colnames(cpm_mat), function(sample) {
  get_cpm_at_percentile(cpm_mat[, sample], rank_pct_all[, sample], cutoffs)
})

# tidy into data.frame
cpm_cutoffs_df <- as.data.frame(t(cpm_cutoffs))
colnames(cpm_cutoffs_df) <- paste0("Percentile_", cutoffs)
cpm_cutoffs_df$Sample <- rownames(cpm_cutoffs_df)

cpm_cutoffs_df
```


----

### 4. Analyze NLR datasets

**📥 Outputs**    
• `NLR_GlobalrankPercentile`: median NLR percentile ranks per condition  
• `NLR_CPM`: mean TMM-normalized CPM counts for NLRs per condition  

```
#!/usr/bin/env Rscript

library(matrixStats)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggtern)

# ----------- set path -------------
setwd('C:/Users/skyun/Desktop/Research/25.Kronos_analysis/06.RNAseq/percentile_ranking_analysis')

# ----------- Process input -------------
data <- read.csv(file = "global_percentiles.txt", sep = "\t", header= TRUE, row.names = 1)
data[is.na(data)] <- 0
data <- as.matrix(data)

# Parse tissue as leading alpha characters before first digit
parse_tissue <- function(nm) {
  # e.g., "leaf1" -> "leaf", "endosperm10" -> "endosperm", "meristem_3" -> "meristem"
  m <- regmatches(nm, regexpr("^[A-Za-z_]+", nm))
  tiss <- ifelse(nchar(m) > 0, gsub("_+$","", m), nm)
  tolower(tiss)
}
samples <- colnames(data)
tissue  <- vapply(samples, parse_tissue, character(1))
meta <- data.frame(sample = samples, tissue = factor(tissue), stringsAsFactors = FALSE)

# ------------------- Gene filtering -------------------
# Remove zero-variance genes; keep top variable genes to stabilize PCA
vars <- matrixStats::rowVars(data)
keep <- vars > 0
data  <- data[keep, , drop=FALSE]
vars <- vars[keep]

# ------------------- Helpers -------------------
plot_pca <- function(pc_df, ve, title){
  ggplot(pc_df, aes(PC1, PC2, color=tissue, label=sample)) +
    geom_point(size=3, alpha=0.95) +
    ggrepel::geom_text_repel(show.legend=FALSE, max.overlaps=12) +
    stat_ellipse(level=0.68, linetype=2, linewidth=0.4) +
    labs(
      title = title,
      x = paste0("PC1 (", scales::percent(ve[1]), ")"),
      y = paste0("PC2 (", scales::percent(ve[2]), ")")
    ) +
    theme_minimal()
}

plot_scree <- function(ve, title){
  ggplot(data.frame(PC=seq_along(ve), VarExp=ve), aes(PC, VarExp)) +
    geom_col() + geom_point() +
    scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
    labs(title=title, x="Principal Component", y="% Variance Explained") +
    theme_minimal()
}

top_loadings <- function(rot, pc=1, n=50){
  v <- rot[, pc]
  n <- min(n, length(v))
  pos <- sort(v, decreasing=TRUE)[seq_len(n)]
  neg <- sort(v, decreasing=FALSE)[seq_len(n)]
  rbind(
    data.frame(gene=names(pos), loading=as.numeric(pos), PC=paste0("PC",pc), direction="positive"),
    data.frame(gene=names(neg), loading=as.numeric(neg), PC=paste0("PC",pc), direction="negative")
  )
}


###################################################################
####################    PCA                   #####################
###################################################################
# No row scaling; default PCA centers features implicitly via prcomp(..., center=TRUE)
pca_cov <- prcomp(t(data), center=TRUE, scale.=FALSE)
ve_cov <- (pca_cov$sdev^2) / sum(pca_cov$sdev^2)

pc_cov <- as.data.frame(pca_cov$x)
pc_cov$sample <- rownames(pc_cov)
pc_cov <- merge(pc_cov, meta, by="sample")

p_cov_scree <- plot_scree(ve_cov, "Scree (Covariance-PCA, raw percentiles)")
p_cov_pc12  <- plot_pca(pc_cov, ve_cov, "Covariance-PCA (raw percentiles)")

p_cov_scree
p_cov_pc12


###################################################################
####################   Pearson Correlations   #####################
###################################################################
cor <- cor(data, method = "pearson")

# Heatmap
pheatmap(cor, 
         main = "Pearson correlations",
         clustering_method = "ward.D2", 
         display_numbers = TRUE)


###################################################################
####################   Tissue-level medians   #####################
###################################################################
mat <- read.delim("global_percentiles.txt", row.names = 1, check.names = FALSE)
mat[is.na(mat)] <- 0
mat <- as.matrix(mat)

# --- Tissue parsing ---
parse_tissue <- function(nm) tolower(sub("_+$","", regmatches(nm, regexpr("^[A-Za-z_]+", nm))))
tissue_vec <- vapply(colnames(mat), parse_tissue, character(1))

leaf_cols <- which(tissue_vec == "leaf")
seed_cols <- which(tissue_vec == "endosperm")
stopifnot(length(leaf_cols) > 0, length(seed_cols) > 0)

# --- Meristem: keep ONLY 9–14  ---
# meristem development in wild type (VanGessel et al.): 5 samples
# + one sample from Xu et al.
meristem_all <- which(tissue_vec == "meristem")
m_nums <- suppressWarnings(as.integer(sub(".*?(\\d+)$","\\1", colnames(mat)[meristem_all])))
meristem_keep <- meristem_all[!is.na(m_nums) & m_nums >= 9 & m_nums <= 14]
if (length(meristem_keep) == 0) stop("No meristem9–14 columns found.")


# --- Per-tissue median percentiles (global percentiles) ---
leaf_med     <- rowMedians(mat[, leaf_cols,     drop = FALSE])
seed_med     <- rowMedians(mat[, seed_cols,     drop = FALSE])
meristem_med <- rowMedians(mat[, meristem_keep, drop = FALSE])

# --- Threshold and sets ---
thr_global <- 0.70   # <--------------------- use 0.40 or 0.70
leaf_expressed      <- rownames(mat)[leaf_med     >= thr_global]
seed_expressed      <- rownames(mat)[seed_med     >= thr_global]
meristem_expressed  <- rownames(mat)[meristem_med >= thr_global]

# --- Min within threshold per tissue ---
min_within_threshold <- function(med_vec, thr, rn) {
  sel <- which(is.finite(med_vec) & med_vec >= thr)
  if (!length(sel)) return(list(gene=NA_character_, value=NA_real_))
  i <- sel[ which.min(med_vec[sel]) ]
  list(gene = rn[i], value = med_vec[i])
}
leaf_min     <- min_within_threshold(leaf_med,     thr_global, rownames(mat))
seed_min     <- min_within_threshold(seed_med,     thr_global, rownames(mat))
meristem_min <- min_within_threshold(meristem_med, thr_global, rownames(mat))

############ Check normalized CPM counts to see if the cutoff makes sense
cat(sprintf("Lowest *within* threshold (thr=%.2f):\n", thr_global))
cat(sprintf("Check normalized counts"))
cat(sprintf("  leaf:     %s (median=%.3f)\n", leaf_min$gene,     leaf_min$value))
cat(sprintf("  seed:     %s (median=%.3f)\n", seed_min$gene,     seed_min$value))
cat(sprintf("  meristem: %s (median=%.3f)\n", meristem_min$gene, meristem_min$value))

# Sanity: minima must be >= thr
stopifnot(is.na(leaf_min$value)     || leaf_min$value     >= thr_global)
stopifnot(is.na(seed_min$value)     || seed_min$value     >= thr_global)
stopifnot(is.na(meristem_min$value) || meristem_min$value >= thr_global)

# --- Venn counts for 3 tissues ---
L <- unique(leaf_expressed)
S <- unique(seed_expressed)
M <- unique(meristem_expressed)

L_only <- setdiff(L, union(S, M))
S_only <- setdiff(S, union(L, M))
M_only <- setdiff(M, union(L, S))

LS     <- setdiff(intersect(L, S), M)  # Leaf ∩ Seed only
LM     <- setdiff(intersect(L, M), S)  # Leaf ∩ Meristem only
SM     <- setdiff(intersect(S, M), L)  # Seed ∩ Meristem only

LSM    <- Reduce(intersect, list(L, S, M))  # all three

venn_segments <- list(
  leaf_only      = L_only,
  seed_only      = S_only,
  meristem_only  = M_only,
  leaf_seed      = LS,
  leaf_meristem  = LM,
  seed_meristem  = SM,
  all_three      = LSM
)

# Print counts for all regions
for (nm in names(venn_segments)) {
  cat(sprintf("%-16s = %d\n", nm, length(venn_segments[[nm]])))
}

# Optional totals
cat("leaf total     =", length(L), "\n")
cat("seed total     =", length(S), "\n")
cat("meristem total =", length(M), "\n")
cat("union total    =", length(unique(c(L, S, M))), "\n")

# Union of all genes
union_genes <- sort(unique(c(L, S, M)))

# Build presence/absence table
presence_df <- data.frame(
  gene = union_genes,
  L = as.integer(union_genes %in% L),
  S = as.integer(union_genes %in% S),
  M = as.integer(union_genes %in% M)
)

# Write to file
write.table(presence_df,
            file = paste0("NLR_presence_matrix.", thr_global, ".txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)


###################################################################
########   Optional: NLR percentile distributions   ###############
###################################################################
med_df <- tibble(
  gene = rownames(mat),
  Leaf = leaf_med,
  Seed = seed_med,
  Meristem = meristem_med
) |>
  pivot_longer(-gene, names_to="tissue", values_to="median_pct")


# Summary table
all_summ <- med_df |>
  group_by(tissue) |>
  summarize(
    n = n(),
    q25 = quantile(median_pct, .25, na.rm=TRUE),
    median = median(median_pct, na.rm=TRUE),
    q75 = quantile(median_pct, .75, na.rm=TRUE),
    q90 = quantile(median_pct, .90, na.rm=TRUE),
    mean = mean(median_pct, na.rm=TRUE),
    sd = sd(median_pct, na.rm=TRUE),
    .groups="drop"
  )
print(all_summ)

# Violin/box (all NLRs)
ggplot(med_df, aes(tissue, median_pct, fill=tissue)) +
  geom_violin(scale="width", trim=FALSE) +
  geom_boxplot(width=0.12, alpha=0.7, outlier.size=0.4) +
  coord_cartesian(ylim=c(0,1)) +
  labs(title="NLR percentile distributions by tissue (all genes)",
       x=NULL, y="Median global percentile (per tissue)") +
  theme_minimal() + theme(legend.position="none")


# --------  Distributions for selected thresholds --------
thr_expr <- 0.70  # <--------------------- use 0.40 or 0.70
expr_df <- med_df |>
  group_by(tissue) |>
  filter(median_pct >= thr_expr) |>
  ungroup()

expr_summ <- expr_df |>
  group_by(tissue) |>
  summarize(
    n = n(),
    q25 = quantile(median_pct, .25, na.rm=TRUE),
    median = median(median_pct, na.rm=TRUE),
    q75 = quantile(median_pct, .75, na.rm=TRUE),
    q90 = quantile(median_pct, .90, na.rm=TRUE),
    mean = mean(median_pct, na.rm=TRUE),
    sd = sd(median_pct, na.rm=TRUE),
    .groups="drop"
  )
print(expr_summ)

ggplot(expr_df, aes(tissue, median_pct, fill=tissue)) +
  geom_violin(scale="width", trim=FALSE) +
  geom_boxplot(width=0.12, alpha=0.7, outlier.size=0.4) +
  coord_cartesian(ylim=c(thr_expr,1)) +
  labs(title=paste0("NLR percentile distributions by tissue (expressed ≥ ", thr_expr, ")"),
       x=NULL, y="Median global percentile (per tissue)") +
  theme_minimal() + theme(legend.position="none")

# Pairwise KS (shape/location) + BH FDR and Cliff’s delta
pairs <- combn(unique(expr_df$tissue), 2, simplify=FALSE)
ks_tab <- do.call(rbind, lapply(pairs, function(ps){
  x <- expr_df$median_pct[expr_df$tissue==ps[1]]
  y <- expr_df$median_pct[expr_df$tissue==ps[2]]
  out <- ks.test(x,y, exact=FALSE)
  data.frame(comp=paste(ps,collapse=" vs "),
             ks_stat=unname(out$statistic),
             p=out$p.value)
}))
ks_tab$p_adj <- p.adjust(ks_tab$p, method="BH")
print(ks_tab)


###################################################################
########################   Ternary plot   #########################
###################################################################
# 1) Assemble scores
scores <- cbind(
  Leaf      = as.numeric(leaf_med),
  Endosperm = as.numeric(seed_med),
  Meristem  = as.numeric(meristem_med)
)
rownames(scores) <- rownames(mat)

# 2) Universe = union of expressed genes
univ <- sort(unique(c(leaf_expressed, seed_expressed, meristem_expressed)))
scores_u <- scores[univ, , drop = FALSE]

# 3) Winner tissue with margin (ties -> "Shared")
best_idx   <- max.col(scores_u, ties.method = "first")
best_val   <- scores_u[cbind(seq_len(nrow(scores_u)), best_idx)]
second_val <- apply(scores_u, 1, function(x) sort(x, decreasing = TRUE)[2])
margin     <- best_val - second_val
winner     <- colnames(scores_u)[best_idx]
winner[margin < 0.05] <- "Shared"

# 4) Normalize to barycentric coords (sum = 1)
xyz <- pmax(scores_u, 0)
xyz <- sweep(xyz, 1, rowSums(xyz), "/")
plot_df <- data.frame(gene = rownames(scores_u), xyz, winner, check.names = FALSE)

# 5) Plot (consistent colors with other panels)
plot_df$winner <- factor(plot_df$winner, levels = c("Leaf","Endosperm","Meristem","Shared"))
cols <- c(Leaf="#8FBC5A", Endosperm="#6AAAE2", Meristem="#C97ABB", Shared="grey60")  # match your palette

ggtern(plot_df, aes(Leaf, Endosperm, Meristem, color = winner)) +
  geom_point(size = 1.6, alpha = 0.65) +
  scale_color_manual(values = cols, name = NULL) +
  labs(title = "Expressed NLRs (median ≥ 0.40): tissue composition",
       T = "Leaf", L = "Endosperm", R = "Meristem") +
  theme_bw() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 10, face = "bold")
  )

ggtern(plot_df, aes(x = Leaf, y = Endosperm, z = Meristem, color = winner)) +
  geom_point(alpha=0.65, size=1.6) +
  scale_color_manual(values = cols, name = NULL) +
  labs(title = "Expressed NLRs (median ≥ 0.40): tissue composition",
       x = "Leaf", y = "Endosperm", z = "Meristem") +
  theme_bw()


###################################################################
################   ITOL visualization   #########################
###################################################################
add_dot1 <- function(ids) sub("(\\.[0-9]+)?$", ".1", ids, perl = TRUE)
agg_fun <- function(x) max(x, na.rm = TRUE) 
prep_vals <- function(scores, ids) {
  v <- setNames(as.numeric(scores), ids)
  v <- v[is.finite(v) & v >= 0.4]
  tapply(v, add_dot1(names(v)), agg_fun)
}

# 0,4: "#FED976", 0.7: "#FC4E2A", 1.0: "#BD0026"
write_itol_gradient3 <- function(filename, label, values,
                                 col_min = "#FED976",
                                 col_mid = "#FC4E2A",
                                 col_max = "#BD0026",
                                 v_min = 0.40, v_mid = 0.70, v_max = 1.00,
                                 sep = " ") {
  con <- file(filename, "w"); on.exit(close(con))
  writeLines(c(
    "DATASET_GRADIENT",
    paste("SEPARATOR", "SPACE"),
    paste("DATASET_LABEL", label),
    "COLOR #000000",
    paste("COLOR_MIN", col_min),
    paste("COLOR_MID", col_mid),
    paste("COLOR_MAX", col_max),
    paste("USER_MIN_VALUE", format(v_min, trim=TRUE)),
    paste("USER_MID_VALUE", format(v_mid, trim=TRUE)),
    paste("USER_MAX_VALUE", format(v_max, trim=TRUE)),
    paste("STRIP_WIDTH", "90"),
    "DATA"
  ), con)
  df <- data.frame(id = names(values), val = as.numeric(values), check.names = FALSE)
  utils::write.table(df, con, sep = sep, row.names = FALSE, col.names = FALSE, quote = FALSE)
}


lv <- prep_vals(leaf_med,     rownames(mat))
write_itol_gradient3("iTOL_leaf_median.high_expression.txt",      "Leaf median percentile",      lv)
sv <- prep_vals(seed_med,     rownames(mat)) 
write_itol_gradient3("iTOL_endosperm_median.high_expression.txt", "Endosperm median percentile", sv)
mv <- prep_vals(meristem_med, rownames(mat))
write_itol_gradient3("iTOL_meristem_median.high_expression.txt",  "Meristem median percentile",  mv)



###################################################################
######   Differential expressions across developmental stages #####
###################################################################
# ------------------------- inputs -------------------------
# ---- CPM (endosperm) ----
cpm <- read.csv("cpm_counts.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

# ---- Pick endosperm/seed columns in the percentile matrix 'mat' and call "expressed" ----
# (Assumes 'mat' is already in memory: genes x samples, values in [0,1])
pc <- 1
seed_cols_pct <- grep("^(seed|endosperm)\\d+$", tolower(colnames(mat)))
thr_expr <- 0.40
seed_selected <- names(seed_med)[seed_med >= thr_expr]   # expressed in endosperm

# Build consecutive pairs: 2->3, 3->4, ...
pairs <- Map(c, 2:8, 3:9)
pair_lab <- function(a,b) paste0("Seed", b, "/Seed", a)


# Collect log2FC for each consecutive pair
out_list <- list()
for (k in seq_along(pairs)) {
  a <- pairs[[k]][1]  # previous stage
  b <- pairs[[k]][2]  # next stage
  colA <- paste0("Seed", a)
  colB <- paste0("Seed", b)
  x <- cpm[seed_selected, colA, drop=FALSE]  # stage a
  y <- cpm[seed_selected, colB, drop=FALSE]  # stage b
  
  log2fc <- log2( (y[ , 1] + pc) / (x[, 1] + pc) )
  out_list[[pair_lab(a,b)]] <- data.frame(
    gene   = rownames(y),
    pair   = pair_lab(a,b),
    log2FC = as.numeric(log2fc),
    stringsAsFactors = FALSE
  )
}

df_long <- dplyr::bind_rows(out_list)
df_long$pair <- factor(df_long$pair,
                       levels = sapply(pairs, function(v) pair_lab(v[1], v[2])))

# ----------------- Plots -----------------
# assign colors in a sequential palette
n_pairs <- length(levels(df_long$pair))
pal <- scales::brewer_pal(palette = "YlGnBu")(n_pairs)   # or try viridis::viridis(n_pairs)
names(pal) <- levels(df_long$pair)

p <- ggplot(df_long, aes(x = log2FC, color = pair)) +
  geom_density(size = 1.1, adjust = 1) +
  xlim(-2, 2) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "grey60") +
  labs(
    title = "Consecutive endosperm stages (Seed b / Seed a)",
    x = "log2 fold-change", y = "Density", color = "Stage pair"
  ) +
  scale_color_manual(values = pal) +
  theme_minimal()
print(p)

summ <- df_long %>%
  group_by(pair) %>%
  summarize(
    n_genes = n(),
    median_log2FC = median(log2FC, na.rm=TRUE),
    IQR_log2FC    = IQR(log2FC, na.rm=TRUE),
    pct_above_neg0p5 = mean(log2FC > -0.5) * 100,  # main metric
    pct_within_2x    = mean(abs(log2FC) < 1) * 100,
    .groups = "drop"
  )
print(summ)


#--------------- meristem -------------------------------
meristem_cols_pct <- grep("^(meristem)\\d+$", tolower(colnames(mat)))
thr_expr <- 0.40
meristem_selected <- names(meristem_med)[meristem_med >= thr_expr]   # expressed in endosperm

# Build consecutive pairs: 2->3, 3->4, ...
pairs <- Map(c, 9:12, 10:13)
pair_lab <- function(a,b) paste0("Meristem", b, "/Meristem", a)

# Collect log2FC for each consecutive pair
out_list <- list()
for (k in seq_along(pairs)) {
  a <- pairs[[k]][1]  # previous stage
  b <- pairs[[k]][2]  # next stage
  colA <- paste0("Meristem", a)
  colB <- paste0("Meristem", b)
  x <- cpm[meristem_selected, colA, drop=FALSE]  # stage a
  y <- cpm[meristem_selected, colB, drop=FALSE]  # stage b
  
  log2fc <- log2( (y[ , 1] + pc) / (x[, 1] + pc) )
  out_list[[pair_lab(a,b)]] <- data.frame(
    gene   = rownames(y),
    pair   = pair_lab(a,b),
    log2FC = as.numeric(log2fc),
    stringsAsFactors = FALSE
  )
}

df_long <- dplyr::bind_rows(out_list)
df_long$pair <- factor(df_long$pair,
                       levels = sapply(pairs, function(v) pair_lab(v[1], v[2])))

# ----------------- Plots -----------------
# assign colors in a sequential palette
n_pairs <- length(levels(df_long$pair))
pal <- scales::brewer_pal(palette = "YlGnBu")(n_pairs)   # or try viridis::viridis(n_pairs)
names(pal) <- levels(df_long$pair)

ggplot(df_long, aes(x = log2FC, color = pair)) +
  geom_density(size = 1.1, adjust = 1) +
  xlim(-2, 2) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "grey60") +
  labs(
    title = "Consecutive endosperm stages (Seed b / Seed a)",
    x = "log2 fold-change", y = "Density", color = "Stage pair"
  ) +
  scale_color_manual(values = pal) +
  theme_minimal()


summ <- df_long %>%
  group_by(pair) %>%
  summarize(
    n_genes = n(),
    median_log2FC = median(log2FC, na.rm=TRUE),
    IQR_log2FC    = IQR(log2FC, na.rm=TRUE),
    pct_above_neg0p5 = mean(log2FC > -0.5) * 100,  # main metric
    pct_within_2x    = mean(abs(log2FC) < 1) * 100,
    .groups = "drop"
  )
print(summ)
```

