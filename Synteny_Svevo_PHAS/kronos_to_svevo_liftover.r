## kronos_to_svevo_liftover.R
## Refactored pipeline: Kronos (query) -> Svevo (target)
cmd_makeblastdb <- "makeblastdb -in Ttu_Svevov1.fa -dbtype nucl -parse_seqids -out Ttu_Svevov1_genome_db"
system(cmd_makeblastdb)

## ========== 0) User-configurable paths ==========
base_prefix   <- "Ttu_Kronos_PHAS_to_Svevo"
out_dir       <- "results_Kronos_PHAS_to_Svevo"

## Query (Kronos)
query_label   <- "Kronos"
query_bed     <- "Ttu_Kronos_PHAS.bed"              # PHAS loci in Kronos BED (4 or 6 cols)
query_fa      <- "Ttu_Kronos_genome_clean.fa"       # Kronos genome FASTA (with .fai alongside)
query_fai     <- paste0(query_fa, ".fai")
query_PHAS_fa <- "Ttu_Kronos_PHAS.fa"               # will be created by bedtools getfasta

## Target (Svevo)
target_label  <- "Svevo"
target_fa     <- "Ttu_Svevov1.fa"                   # Svevo genome FASTA (with .fai alongside)
target_fai    <- paste0(target_fa, ".fai")
target_db     <- "Ttu_Svevov1_genome_db"            # makeblastdb-prepared DB name for Svevo

## ========== 1) Prep output dir ==========
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ========== 2) Extract Kronos PHAS sequences ==========
## (If you already have query_PHAS_fa, this will overwrite it)
cmd_getfasta <- sprintf(
  "bedtools getfasta -fi %s -fo %s -bed %s -name",
  shQuote(query_fa), shQuote(query_PHAS_fa), shQuote(query_bed)
)
cat("[CMD]", cmd_getfasta, "\n")
system(cmd_getfasta)

## ========== 3) BLAST Kronos→Svevo ==========
blast_tab <- file.path(out_dir, paste0(base_prefix, "_blast.tab"))
cmd_blast <- sprintf(
  "blastn -query %s -db %s -outfmt 6 -out %s -num_threads 12",
  shQuote(query_PHAS_fa), shQuote(target_db), shQuote(blast_tab)
)
cat("[CMD]", cmd_blast, "\n")
system(cmd_blast)

## ========== 4) Read BLAST + query BED and pick top hits (simple, fast) ==========
cn <- c("query","subject","identity","alignment_length","mismatches","gap_opens",
        "q_start","q_end","s_start","s_end","evalue","bit_score")

df <- read.table(blast_tab, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                 col.names = cn, quote = "", comment.char = "", check.names = FALSE)
if (nrow(df) == 0) stop("No BLAST hits found.")

read_bed4_or_6 <- function(path) {
  x <- read.table(path, header = FALSE, sep = "\t", quote = "", comment.char = "",
                  stringsAsFactors = FALSE, fill = TRUE, check.names = FALSE)
  if (!ncol(x) %in% c(4, 6)) stop(sprintf("Expected 4 or 6 columns in %s, found %d.", path, ncol(x)))
  if (ncol(x) == 4) {
    colnames(x) <- c("chrom","start","end","name")
  } else {
    colnames(x) <- c("chrom","start","end","name","score","strand")
    x <- x[, c("chrom","start","end","name")]
  }
  x$chrom <- trimws(x$chrom)
  x$name  <- trimws(x$name)
  x$start <- as.integer(x$start)
  x$end   <- as.integer(x$end)
  x
}
bed_q <- read_bed4_or_6(query_bed)

## --- match BLAST query to BED name (strip bedtools suffix like "::1A:123-456") ---
df$query_base    <- sub("::.*$", "", df$query)
bed_q$name_base  <- sub("::.*$", "", bed_q$name)

idx <- match(df$query_base, bed_q$name_base)   # vectorized, memory-light
df$chrom_query <- bed_q$chrom[idx]
df$chrom_query <- trimws(df$chrom_query)

## standardize labels on both sides for comparison ("1A", "Chr1A", "Ttu_1A" -> "1A")
std_chr <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y <- sub("^CHR",  "", y)
  y <- sub("^TTU_", "", y)
  y <- gsub("[^0-9A-Z]+", "", y)
  y <- ifelse(y %in% c("UN","UNPLACED"), "UN", y)
  y
}
df$subject_std     <- std_chr(df$subject)
df$chrom_query_std <- std_chr(df$chrom_query)

## label same-chromosome hits
df$Chr_same <- ifelse(!is.na(df$chrom_query_std) & df$chrom_query_std == df$subject_std, "Same", "Different")

## rank within each query: Same first, then bit_score, then identity
same_flag <- as.integer(df$Chr_same == "Same")
bs <- suppressWarnings(as.numeric(df$bit_score)); bs[is.na(bs)] <- -Inf
id <- suppressWarnings(as.numeric(df$identity));  id[is.na(id)] <- -Inf

ord <- order(df$query, -same_flag, -bs, -id)
df1 <- df[ord, , drop = FALSE]                         # keep name df1 if you use it elsewhere
top_hits <- df1[!duplicated(df1$query), , drop = FALSE]

## diagnostics
cat("Queries total:", length(unique(df$query)), "\n")
cat("chrom_query known for:", sum(!is.na(df$chrom_query)), "rows\n")
cat("Top hits chosen:", nrow(top_hits), "\n")
cat("  — Same-chromosome:", sum(top_hits$Chr_same == "Same", na.rm = TRUE), "\n")

## ========== 5) Write BED on the TARGET (Svevo) side ==========
start0 <- pmin(as.numeric(top_hits$s_start), as.numeric(top_hits$s_end)) - 1L
end1   <- pmax(as.numeric(top_hits$s_start), as.numeric(top_hits$s_end))
strand <- ifelse(as.numeric(top_hits$s_start) <= as.numeric(top_hits$s_end), "+", "-")

bed_t <- data.frame(
  chrom  = top_hits$subject,   # raw subject; Step 6 will normalize to FASTA headers
  start  = start0,
  end    = end1,
  name   = top_hits$query,
  score  = ".",
  strand = strand,
  stringsAsFactors = FALSE
)

final_bed <- sub("\\.tab$", ".bed", blast_tab)
write.table(bed_t, file = final_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("BED written:", final_bed, "\n")

## ========== 6) Extract mapped sequences from Svevo (no remap needed) ==========
## Read the BED we just wrote (chrom = BLAST subject)
bed_out <- read.table(final_bed, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                      quote = "", comment.char = "", check.names = FALSE)

## sanity check: all chromosomes in BED must exist in the target FASTA headers
fa_hdrs <- system(sprintf("grep '^>' %s | sed 's/^>//; s/\\s.*$//'", shQuote(target_fa)),
                  intern = TRUE)
missing <- setdiff(unique(bed_out[[1]]), fa_hdrs)
if (length(missing)) {
  warning(sprintf("These chromosomes are not present in FASTA: %s",
                  paste(missing, collapse = ", ")))
}

## write temp and run bedtools (strand-aware)
tmp_bed <- tempfile(fileext = ".bed")
write.table(bed_out, file = tmp_bed, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

extracted_fa <- sub("\\.bed$", ".fa", final_bed)
cmd_get_t <- sprintf("bedtools getfasta -fi %s -bed %s -fo %s -s",
                     shQuote(target_fa), shQuote(tmp_bed), shQuote(extracted_fa))
cat("[CMD]", cmd_get_t, "\n")
system(cmd_get_t)

upper_fa <- sub("\\.bed$", "_Up.fa", final_bed)
cmd_upper <- sprintf("awk '{if ($0 ~ /^>/) print $0; else print toupper($0)}' %s > %s",
                     shQuote(extracted_fa), shQuote(upper_fa))
cat("[CMD]", cmd_upper, "\n")
system(cmd_upper)

## ========== 7) Build join table for plotting (Kronos vs Svevo) ==========
to_bed4 <- function(path) {
  x <- read_bed4_or_6(path)
  x[, c("chrom","start","end","name")]
}
bed_query  <- to_bed4(query_bed)      # Kronos
bed_target <- to_bed4(final_bed)      # Svevo (from BLAST)

## normalize IDs by stripping "::..." so BLAST queries match BED names
bed_query$name_base  <- sub("::.*$", "", bed_query$name)
bed_target$name_base <- sub("::.*$", "", bed_target$name)

bed_query$genome     <- query_label
bed_target$genome    <- target_label
bed_query$midpoint   <- (bed_query$start  + bed_query$end)  / 2
bed_target$midpoint  <- (bed_target$start + bed_target$end) / 2

merged <- merge(bed_query, bed_target,
                by = "name_base", suffixes = c("_query","_target"),
                all = FALSE, sort = FALSE)

cat("Rows in Kronos (query) BED:", nrow(bed_query),  "\n")
cat("Rows in Svevo (target) BED:", nrow(bed_target), "\n")
cat("Inner-join rows:", nrow(merged), "\n")
if (nrow(merged) == 0) stop("No rows after inner-join; cannot plot.")

## ========== 8) Normalize chromosomes for plotting and filter UN ==========
## With matching headers, just trim; optionally collapse UN/Un/ChrUn to "UN"
norm_chr <- function(x) {
  y <- trimws(as.character(x))
  y <- ifelse(toupper(y) %in% c("UN","CHRUN","UNPLACED"), "UN", y)
  y
}
merged$chr_query_std  <- norm_chr(merged$chrom_query)
merged$chr_target_std <- norm_chr(merged$chrom_target)

## read .fai and keep only used (non-UN) chromosomes
fai_q <- read.table(query_fai, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                    col.names = c("chr","length","offset","linebases","linewidth"))
fai_t <- read.table(target_fai, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                    col.names = c("chr","length","offset","linebases","linewidth"))
fai_q$chr_std <- norm_chr(fai_q$chr)
fai_t$chr_std <- norm_chr(fai_t$chr)

merged <- subset(merged, chr_query_std != "UN" & chr_target_std != "UN")
fai_q  <- subset(fai_q,  chr_std != "UN")
fai_t  <- subset(fai_t,  chr_std != "UN")

## consistent ordering: extract the first number anywhere, then arm letters (e.g., "1A", "12D")
order_chr <- function(x) {
  x   <- as.character(x)
  num <- suppressWarnings(as.integer(sub(".*?([0-9]+).*", "\\1", x)))
  arm <- sub(".*?[0-9]+([A-Za-z]+).*", "\\1", x)
  x[order(num, arm, na.last = NA)]
}

used_chrs <- order_chr(sort(unique(c(merged$chr_query_std, merged$chr_target_std))))

fai_q <- fai_q[fai_q$chr_std %in% used_chrs, , drop = FALSE]
fai_t <- fai_t[fai_t$chr_std %in% used_chrs, , drop = FALSE]
fai_q$chr <- fai_q$chr_std; fai_q$chr_std <- NULL
fai_t$chr <- fai_t$chr_std; fai_t$chr_std <- NULL

cat("Used chrs:", paste(used_chrs, collapse = ", "), "\n")

## ========== 9) Two-panel link plot (Kronos on top, Svevo bottom) ==========
col_query_bar <- "#007C84"  # Kronos
col_target_bar <- "#CC5146" # Svevo
col_link      <- "#222222"
col_grid_v    <- "grey88"
col_grid_h    <- "grey92"

t_y1 <- 0.60; t_y2 <- 0.68; t_y_mid <- (t_y1 + t_y2) / 2   # Svevo (bottom)
q_y1 <- 1.32; q_y2 <- 1.40; q_y_mid <- (q_y1 + q_y2) / 2   # Kronos (top)

n_chr <- length(used_chrs)
ncol  <- 3
nrow  <- ceiling(n_chr / ncol)

pdf_path <- file.path(out_dir, paste0(base_prefix, "_liftover_connections.pdf"))
pdf(pdf_path, width = 18, height = max(8, nrow * 2.8), useDingbats = FALSE)
op <- par(no.readonly = TRUE)
on.exit({par(op); dev.off()}, add = TRUE)

par(oma = c(0.6, 0.6, 1.8, 0.6))
par(mfrow = c(nrow, ncol), mar = c(3.1, 5.2, 2.9, 1.2), mgp = c(2.2, 0.7, 0), xaxs = "i", yaxs = "i")

for (ch in used_chrs) {
  q_len <- fai_q$length [match(ch, fai_q$chr)]
  t_len <- fai_t$length[match(ch, fai_t$chr)]
  sub   <- merged[merged$chr_query_std == ch | merged$chr_target_std == ch, , drop = FALSE]

  x_max <- suppressWarnings(max(c(q_len, t_len), na.rm = TRUE))
  if (!is.finite(x_max)) {
    x_max <- suppressWarnings(max(c(sub$midpoint_query, sub$midpoint_target), na.rm = TRUE))
    if (!is.finite(x_max)) x_max <- 1
  }

  plot(NA, xlim = c(0, x_max), ylim = c(0.4, 1.55), xaxt = "n", yaxt = "n",
       xlab = "Genomic position (Mb)", ylab = "",
       main = paste0("Chr", ch), cex.main = 1.35, font.main = 2, bty = "n")

  at <- pretty(c(0, x_max), n = 6)
  abline(v = at, col = col_grid_v, lwd = 0.9)
  abline(h = c(t_y1, t_y2, q_y1, q_y2), col = col_grid_h, lwd = 0.9)
  axis(1, at = at, labels = formatC(at / 1e6, format = "f", digits = 1), cex.axis = 0.95)
  axis(2, at = c(t_y_mid, q_y_mid), labels = c(target_label, query_label), las = 1, cex.axis = 1.05)

  if (is.finite(t_len)) rect(0, t_y1, t_len, t_y2, col = adjustcolor(col_target_bar, 0.55), border = "black", lwd = 0.9)
  else                 rect(0, t_y1, x_max, t_y2, col = adjustcolor(col_target_bar, 0.20), border = "black", lty = 3, lwd = 0.9)
  if (is.finite(q_len)) rect(0, q_y1, q_len, q_y2, col = adjustcolor(col_query_bar, 0.55), border = "black", lwd = 0.9)
  else                 rect(0, q_y1, x_max, q_y2, col = adjustcolor(col_query_bar, 0.20), border = "black", lty = 3, lwd = 0.9)

  ok <- is.finite(sub$midpoint_target) & is.finite(sub$midpoint_query)
  if (any(ok)) {
    sub_ok <- sub[ok, , drop = FALSE]
    if (is.finite(t_len)) sub_ok$midpoint_target <- pmin(pmax(sub_ok$midpoint_target, 0), t_len)
    if (is.finite(q_len)) sub_ok$midpoint_query  <- pmin(pmax(sub_ok$midpoint_query,  0), q_len)
    segments(sub_ok$midpoint_target, t_y_mid, sub_ok$midpoint_query,  q_y_mid,
             col = adjustcolor(col_link, 0.85), lwd = 1.2)
    points(sub_ok$midpoint_target, rep(t_y_mid, nrow(sub_ok)), pch = 16, cex = 0.95, col = col_target_bar)
    points(sub_ok$midpoint_query,  rep(q_y_mid, nrow(sub_ok)), pch = 16, cex = 0.95, col = col_query_bar)
  }
}
mtext("PHAS cross-genome mapping (Kronos → Svevo)", side = 3, line = 0.6, outer = TRUE, cex = 1.22, font = 2)
dev.off()
cat("Wrote:", pdf_path, "\n")


## ========== 10) Export liftover table ==========
lift_out <- subset(merged,
                   !is.na(chr_query_std) & chr_query_std != "UN" &
                   !is.na(chr_target_std) & chr_target_std != "UN" &
                   is.finite(start_query) & is.finite(end_query) &
                   is.finite(start_target) & is.finite(end_target))

ord <- order(
  as.integer(sub("^([0-9]+).*", "\\1", lift_out$chr_target_std)),
  sub("^[0-9]+", "", lift_out$chr_target_std),
  lift_out$start_target,
  na.last = TRUE
)
lift_out <- lift_out[ord, , drop = FALSE]

out_tsv <- file.path(out_dir, paste0(base_prefix, "_KronosPHAS_lifted_to_Svevo.tsv"))
tab <- data.frame(
  PHAS_id          = lift_out$name_base,
  Kronos_chr       = lift_out$chr_query_std,
  Kronos_start     = as.integer(lift_out$start_query),
  Kronos_end       = as.integer(lift_out$end_query),
  Kronos_midpoint  = as.integer(lift_out$midpoint_query),
  Svevo_chr        = lift_out$chr_target_std,
  Svevo_start      = as.integer(lift_out$start_target),
  Svevo_end        = as.integer(lift_out$end_target),
  Svevo_midpoint   = as.integer(lift_out$midpoint_target),
  stringsAsFactors = FALSE
)
write.table(tab, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote liftover table:\n  ", out_tsv, "\n")
