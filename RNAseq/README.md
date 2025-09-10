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

#--------------- Group expression --------------------------------------
split_idx <- split(seq_len(ncol(nlr_cpm)), samples_clean$grp)

rowMed_by_group <- function(mat, idx) {
  out <- do.call(cbind, lapply(idx, function(ix) matrixStats::rowMedians(mat[, ix, drop = FALSE])))
  colnames(out) <- names(idx); rownames(out) <- rownames(mat); out
}

CPM_grp   <- sapply(split_idx, function(ix) rowMeans(cpm_mat[, ix, drop=FALSE]))

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
  as_tibble(tibble::rownames_to_column(as.data.frame(CPM_grp[nlr_ids, groups_keep, drop=FALSE]),
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

