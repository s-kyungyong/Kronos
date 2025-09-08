# ========== 1. Retrieve and lift ==========
system("bedtools getfasta -fi Ttu_Svevov1.fa -fo Ttu_Svevov1_PHAS.fa -bed Ttu_Svevov1_PHAS.bed -name")

base_prefix="Ttu_Svevov1_PHAS_to_Kronos"
genome_db="Ttu_Kronos_genome_db"
genome_fa="Ttu_Kronos_genome_clean.fa"
PHAS="Ttu_Svevov1_PHAS.fa"
query_bed_path="Ttu_Svevov1_PHAS.bed"
query_fai_path="Ttu_Svevov1.fa.fai"
query_label="Svevo"
target_label="Kronos"
out_dir="results"

mkdir -p "$out_dir"

blast_result="${out_dir}/${base_prefix}_blast.tab"
final_bed="${blast_result%.tab}.bed"
extracted_fa="${final_bed%.bed}.fa"
upper_fa="${final_bed%.bed}_Up.fa"

system("blastn -query Ttu_Svevov1_PHAS.fa -db Ttu_Kronos_genome_db -outfmt 6 -out $blast_result -num_threads 12")

# ========== 2. Read and filter ==========
## ---- Inputs (from your shell vars) ----
blast_result <- "results/Ttu_Svevov1_PHAS_to_Kronos_blast.tab"
query_bed_path <- "Ttu_Svevov1_PHAS.bed"

## ---- 1) Read BLAST + BED (base R) ----
df <- read.table(blast_result, sep="\t", header=FALSE, stringsAsFactors=FALSE,
                 col.names=c("query","subject","identity","alignment_length","mismatches",
                             "gap_opens","q_start","q_end","s_start","s_end","evalue","bit_score"))

bed <- read.table(query_bed_path, sep="\t", header=FALSE, stringsAsFactors=FALSE,col.names=c("chrom","start","end","name"))

## ---- 2) Try joining query->chrom from BED ----
df1 <- merge(df, bed[, c("name","chrom")], by.x="query", by.y="name", all.x=TRUE)

## If that didn’t match (chrom all NA), normalize key by stripping '::...'
if (all(is.na(df1$chrom))) {
    print("Uhum....")
  df$query_base <- sub("::.*$", "", df$query)
  bed$name_base <- sub("::.*$", "", bed$name)
  df1 <- merge(df, bed[, c("name_base","chrom")],
               by.x="query_base", by.y="name_base", all.x=TRUE)
}

## ---- 3) Coalesce chrom.x / chrom.y -> chrom_join; fallback parse from query ----
cx <- if ("chrom.x" %in% names(df1)) df1$chrom.x else rep(NA_character_, nrow(df1))
cy <- if ("chrom.y" %in% names(df1)) df1$chrom.y else if ("chrom" %in% names(df1)) df1$chrom else rep(NA_character_, nrow(df1))
df1$chrom_join <- ifelse(is.na(cx), cy, cx)

# If still NA, parse like '...::Ttu_1A:123..456' -> 'Ttu_1A'
need <- is.na(df1$chrom_join)
if (any(need)) {
  parsed <- sub("^.*::([^:]+):.*$", "\\1", df1$query[need])
  parsed[!grepl("_", parsed)] <- NA_character_  # require prefix_like 'Ttu_1A'
  df1$chrom_join[need] <- parsed
}
df1$chrom_join <- trimws(df1$chrom_join)

## ---- 4) Standardize subject to same prefix as chrom_join (e.g., 'Ttu_1A') ----
# infer prefix from first non-NA chrom_join (e.g., 'Ttu' from 'Ttu_1A')
prefix <- ""
nz <- which(!is.na(df1$chrom_join))
if (length(nz) > 0) {
  ex <- df1$chrom_join[nz[1]]
  pr <- sub("_.*$", "", ex)
  if (nzchar(pr)) prefix <- pr
}
df1$subject <- trimws(df1$subject)
df1$subject_std <- paste0(prefix, "_", df1$subject)

## ---- 5) Homology flag (same chromosome) ----
df1$Chr_same <- ifelse(!is.na(df1$chrom_join) & df1$subject_std == df1$chrom_join,
                       "Same", "Different")

## ---- 6) Pick top hit per query, preferring Same-chr; fallback to best overall ----
ord <- order(df1$query, -(as.numeric(df1$bit_score)), -(as.numeric(df1$identity)))
df1_sorted <- df1[ord, , drop=FALSE]

# First, best Same-chr per query
same <- df1_sorted[df1_sorted$Chr_same == "Same", , drop=FALSE]
top_same <- same[!duplicated(same$query), , drop=FALSE]

# For queries without a Same-chr hit, fill with best overall
have_same <- unique(top_same$query)
missing_q <- setdiff(unique(df1_sorted$query), have_same)
fallback <- df1_sorted[df1_sorted$query %in% missing_q, , drop=FALSE]
top_fallback <- fallback[!duplicated(fallback$query), , drop=FALSE]

top_hits <- rbind(top_same, top_fallback)

## ---- 7) Write BED of target loci (subject coords), strand-aware ----
start0 <- pmin(as.numeric(top_hits$s_start), as.numeric(top_hits$s_end)) - 1L
end1   <- pmax(as.numeric(top_hits$s_start), as.numeric(top_hits$s_end))
strand <- ifelse(as.numeric(top_hits$s_start) <= as.numeric(top_hits$s_end), "+", "-")

bed_out <- data.frame(
  chrom  = top_hits$subject_std,   # Kronos chromosome standardized to 'Ttu_#A/#B'
  start  = start0,
  end    = end1,
  name   = top_hits$query,
  score  = ".",
  strand = strand,
  stringsAsFactors = FALSE
)

out_bed <- sub("\\.tab$", ".bed", blast_result)
write.table(bed_out, file=out_bed, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

## ---- 8) Quick report ----
cat("Queries total:", length(unique(df$query)), "\n")
cat("chrom_join known for:", sum(!is.na(df1$chrom_join)), "rows\n")
cat("Top hits chosen:", nrow(top_hits), "\n")
cat("  — Same-chromosome:", sum(top_hits$Chr_same == "Same", na.rm=TRUE), "\n")
cat("BED written:", out_bed, "\n")

#Queries total: 9104 
#chrom_join known for: 20750492 rows
#Top hits chosen: 9104 
#  — Same-chromosome: 9080 
#BED written: results/Ttu_Svevov1_PHAS_to_Kronos_blast.bed 

# ========== 2.5 Extract sequence ==========

out_dir <- "results"
base_prefix <- "Ttu_Svevov1_PHAS_to_Kronos"

final_bed   <- "results/Ttu_Svevov1_PHAS_to_Kronos_blast.bed"  # your BED path
genome_fa   <- "Ttu_Kronos_genome_clean.fa"
extracted_fa <- sub("\\.bed$", ".fa", final_bed)

# 1) Read BED (6 columns: chrom, start, end, name, score, strand)
bed <- read.table(final_bed, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
if (ncol(bed) < 3) stop("BED has fewer than 3 columns")

# 2) Map chrom names to FASTA style
map_chrom_to_fasta <- function(x) {
  x <- trimws(x)
  if (grepl("^Chr", x)) return(x)               # already matches FASTA
  y <- sub("^Ttu_", "", x, ignore.case = TRUE)  # strip 'Ttu_' if present
  if (grepl("^[1-7][AB]$", y)) return(paste0("Chr", y))
  if (tolower(y) == "un") return("ChrUn")
  # Fall back: if it looks like "1A" without Ttu_, still fix it
  if (grepl("^[1-7][AB]$", x)) return(paste0("Chr", x))
  # Otherwise return original (will likely be absent from FASTA)
  x
}

old_chr <- bed[[1]]
new_chr <- vapply(old_chr, map_chrom_to_fasta, character(1))
bed[[1]] <- new_chr

# 3) Optional: sanity check against FASTA headers
fa_hdrs <- system(sprintf("grep '^>' %s | sed 's/^>//'", shQuote(genome_fa)),
                  intern = TRUE)
fa_hdrs <- sub("\\s.*$", "", fa_hdrs)  # keep only the first token on header line

missing <- setdiff(unique(bed[[1]]), fa_hdrs)
if (length(missing)) {
  warning(sprintf("These chromosomes are not present in FASTA: %s",
                  paste(missing, collapse = ", ")))
}

# 4) Write a temp BED and run bedtools
tmp_bed <- tempfile(fileext = ".bed")
write.table(bed, file = tmp_bed, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

cmd <- sprintf("bedtools getfasta -fi %s -bed %s -fo %s -s",
               shQuote(genome_fa), shQuote(tmp_bed), shQuote(extracted_fa))
cat(cmd, "\n")
system(cmd)

blast_result <- file.path(out_dir, paste0(base_prefix, "_blast.tab"))
final_bed <- sub("\\.tab$", ".bed", blast_result)
extracted_fa <- sub("\\.bed$", ".fa", final_bed)
upper_fa <- sub("\\.bed$", "_Up.fa", final_bed)


system(sprintf("awk '{if ($0 ~ /^>/) print $0; else print toupper($0)}' %s > %s", extracted_fa, upper_fa))
# ========== 3. Plotting (refactored) ==========
## -------- 1) Read inputs --------
read_bed4_or_6 <- function(path) {
  x <- read.table(path, header = FALSE, sep = "\t", quote = "", comment.char = "",
                  stringsAsFactors = FALSE, fill = TRUE, check.names = FALSE)
  if (!ncol(x) %in% c(4, 6)) stop(sprintf("Expected 4 or 6 columns in %s, found %d.", path, ncol(x)))
  if (ncol(x) == 4) {
    colnames(x) <- c("chr","start","end","name")
  } else {
    colnames(x) <- c("chr","start","end","name","score","strand")
    x <- x[, c("chr","start","end","name")]
  }
  x$chr   <- trimws(x$chr)
  x$name  <- trimws(x$name)
  x$start <- as.integer(x$start)
  x$end   <- as.integer(x$end)
  x
}

query_bed  <- read_bed4_or_6(query_bed_path)
mapped_bed <- read_bed4_or_6(final_bed)

## normalize names by stripping "::..." suffix on both sides (join key)
normalize_names <- TRUE
if (normalize_names) {
  query_bed$name_base  <- sub("::.*$", "", query_bed$name)
  mapped_bed$name_base <- sub("::.*$", "", mapped_bed$name)
  by_key <- "name_base"
} else {
  by_key <- "name"
}

## Add metadata & midpoints
query_bed$genome    <- query_label
mapped_bed$genome   <- target_label
query_bed$midpoint  <- (query_bed$start  + query_bed$end)  / 2
mapped_bed$midpoint <- (mapped_bed$start + mapped_bed$end) / 2

## Inner join by name
merged <- merge(query_bed, mapped_bed,
                by = by_key, suffixes = c("_query","_target"),
                all = FALSE, sort = FALSE)

cat("Rows in query_bed:", nrow(query_bed), "\n")
cat("Rows in mapped_bed:", nrow(mapped_bed), "\n")
cat("Inner-join rows:", nrow(merged), "\n")

if (nrow(merged) == 0) {
  qn <- unique(query_bed[[by_key]]); mn <- unique(mapped_bed[[by_key]])
  cat("Example only in query:\n");  print(utils::head(setdiff(qn, mn), 5))
  cat("Example only in mapped:\n"); print(utils::head(setdiff(mn, qn), 5))
  stop("No rows after inner-join; cannot plot.")
}

## -------- 2) FAI and chromosome normalization --------
query_fai  <- read.table(query_fai_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
target_fai <- read.table(paste0(genome_fa, ".fai"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(query_fai)  <- c("chr","length","offset","linebases","linewidth")
colnames(target_fai) <- c("chr","length","offset","linebases","linewidth")
query_fai$genome  <- query_label
target_fai$genome <- target_label

## Normalize chromosome labels to 1A/1B/... and drop UN/Unplaced
norm_chr <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y <- sub("^CHR",  "", y)       # "Chr1A" -> "1A"
  y <- sub("^TTU_", "", y)       # "Ttu_1A" -> "1A"
  y <- gsub("[^0-9A-Z]+", "", y) # keep only digits and A–Z
  y <- ifelse(grepl("^(UN|UNPLACED)$", y, ignore.case = TRUE), "UN", y)
  y
}

merged$chr_query_std  <- norm_chr(merged$chr_query)
merged$chr_target_std <- norm_chr(merged$chr_target)
query_fai$chr_std     <- norm_chr(query_fai$chr)
target_fai$chr_std    <- norm_chr(target_fai$chr)

## Drop UN everywhere
merged     <- subset(merged, chr_query_std != "UN" & chr_target_std != "UN")
query_fai  <- subset(query_fai,  chr_std != "UN")
target_fai <- subset(target_fai, chr_std != "UN")

## Keep only chromosomes that actually appear (std labels)
used_chrs <- sort(unique(c(merged$chr_query_std, merged$chr_target_std)))

## Helper to order as 1A,1B,2A,2B,...
order_chr <- function(x) {
  num <- as.integer(sub("^([0-9]+).*", "\\1", x))
  arm <- sub("^[0-9]+", "", x)
  o <- order(num, arm, na.last = NA)
  x[o]
}
used_chrs <- order_chr(used_chrs)

## Trim FAIs to used chrs and promote std labels as operative 'chr'
query_fai  <- query_fai [query_fai$chr_std  %in% used_chrs, , drop = FALSE]
target_fai <- target_fai[target_fai$chr_std %in% used_chrs, , drop = FALSE]
query_fai$chr  <- query_fai$chr_std;   query_fai$chr_std  <- NULL
target_fai$chr <- target_fai$chr_std;  target_fai$chr_std <- NULL

cat("Used chrs (std):", paste(used_chrs, collapse = ", "), "\n")
cat("FAI kept -> query:", nrow(query_fai), " target:", nrow(target_fai), "\n")

## -------- 3) Panel plotting (two-column layout) --------
# y-positions for the two genome “bars”
## -------- Plot liftover connections (3 cols; darker colors; thinner bars; bigger titles) --------
# Palette (darker teal/salmon)
col_query_bar <- "#007C84"  # Svevo bar
col_target_bar <- "#CC5146" # Kronos bar
col_query_pt  <- col_query_bar
col_target_pt <- col_target_bar
col_link      <- "#222222"
col_grid_v    <- "grey88"
col_grid_h    <- "grey92"

# Bar geometry (thinner)
t_y1 <- 0.60; t_y2 <- 0.68; t_y_mid <- (t_y1 + t_y2) / 2   # Kronos (bottom)
q_y1 <- 1.32; q_y2 <- 1.40; q_y_mid <- (q_y1 + q_y2) / 2   # Svevo  (top)

n_chr <- length(used_chrs)
ncol  <- 3
nrow  <- ceiling(n_chr / ncol)

pdf_path <- file.path(out_dir, paste0(base_prefix, "_liftover_connections.pdf"))
pdf(pdf_path, width = 18, height = max(8, nrow * 2.8), useDingbats = FALSE)

op <- par(no.readonly = TRUE)
on.exit({par(op); dev.off()}, add = TRUE)

# Outer margins for safe global title placement
par(oma = c(0.6, 0.6, 1.8, 0.6))  # bottom, left, top, right
par(mfrow = c(nrow, ncol), mar = c(3.1, 5.2, 2.9, 1.2), mgp = c(2.2, 0.7, 0), xaxs = "i", yaxs = "i")

for (ch in used_chrs) {
  q_len <- query_fai$length [match(ch, query_fai$chr)]
  t_len <- target_fai$length[match(ch, target_fai$chr)]

  # subset links for this chromosome (std labels)
  sub <- merged[merged$chr_query_std == ch | merged$chr_target_std == ch, , drop = FALSE]

  # x-range
  x_max <- suppressWarnings(max(c(q_len, t_len), na.rm = TRUE))
  if (!is.finite(x_max)) {
    x_max <- suppressWarnings(max(c(sub$midpoint_query, sub$midpoint_target), na.rm = TRUE))
    if (!is.finite(x_max)) x_max <- 1
  }

  # panel
  plot(NA, xlim = c(0, x_max), ylim = c(0.4, 1.55), xaxt = "n", yaxt = "n",
       xlab = "Genomic position (Mb)", ylab = "",
       main = paste0("Chr", ch), cex.main = 1.35, font.main = 2, bty = "n")

  # grid & axes
  at <- pretty(c(0, x_max), n = 6)
  abline(v = at, col = col_grid_v, lwd = 0.9)
  abline(h = c(t_y1, t_y2, q_y1, q_y2), col = col_grid_h, lwd = 0.9)
  axis(1, at = at, labels = formatC(at / 1e6, format = "f", digits = 1), cex.axis = 0.95)
  axis(2, at = c(t_y_mid, q_y_mid), labels = c(target_label, query_label), las = 1, cex.axis = 1.05)

  # bars (thinner; slightly more opaque)
  if (is.finite(t_len)) {
    rect(0, t_y1, t_len, t_y2, col = adjustcolor(col_target_bar, 0.55), border = "black", lwd = 0.9)
  } else {
    rect(0, t_y1, x_max, t_y2, col = adjustcolor(col_target_bar, 0.20), border = "black", lty = 3, lwd = 0.9)
  }
  if (is.finite(q_len)) {
    rect(0, q_y1, q_len, q_y2, col = adjustcolor(col_query_bar, 0.55), border = "black", lwd = 0.9)
  } else {
    rect(0, q_y1, x_max, q_y2, col = adjustcolor(col_query_bar, 0.20), border = "black", lty = 3, lwd = 0.9)
  }

  # links & colored points
  ok <- is.finite(sub$midpoint_target) & is.finite(sub$midpoint_query)
  if (any(ok)) {
    sub_ok <- sub[ok, , drop = FALSE]
    if (is.finite(t_len)) sub_ok$midpoint_target <- pmin(pmax(sub_ok$midpoint_target, 0), t_len)
    if (is.finite(q_len)) sub_ok$midpoint_query  <- pmin(pmax(sub_ok$midpoint_query,  0), q_len)

    segments(sub_ok$midpoint_target, t_y_mid,
             sub_ok$midpoint_query,  q_y_mid,
             col = adjustcolor(col_link, 0.85), lwd = 1.2)
    points(sub_ok$midpoint_target, rep(t_y_mid, nrow(sub_ok)), pch = 16, cex = 0.95, col = col_target_pt)
    points(sub_ok$midpoint_query,  rep(q_y_mid, nrow(sub_ok)), pch = 16, cex = 0.95, col = col_query_pt)
  }
}

# clean global title in the outer margin (visible)
mtext("PHAS cross-genome mapping", side = 3, line = 0.6, outer = TRUE, cex = 1.22, font = 2)

# close device (also guarded by on.exit)
dev.off()
cat("Wrote: ", pdf_path, "\n")

####
## ===== OmicCircos cross-genome links (segAnglePo-safe: 5-column seg.frame) =====
stopifnot(requireNamespace("OmicCircos", quietly = TRUE))
library(OmicCircos)

# 1) Choose chromosomes present in links & FAIs; order 1A,1B,2A,2B,...
order_chr <- function(x) x[order(as.integer(sub("^([0-9]+).*","\\1", x)),
                                 sub("^[0-9]+","", x), na.last = NA)]
used_t <- order_chr(intersect(unique(merged$chr_target_std), unique(target_fai$chr)))
used_q <- order_chr(intersect(unique(merged$chr_query_std),  unique(query_fai$chr)))
stopifnot(length(used_t) > 0, length(used_q) > 0)

# 2) Build 5-column segment frames (required 1–3, optional 4–5)
#    (col4 = display name; col5 = group label; content can be anything)
Tpref <- "K_"              # Kronos block prefix
Qpref <- "S_"              # Svevo  block prefix

seg_target <- data.frame(
  seg.name      = paste0(Tpref, used_t),                                   # 1
  seg.Start     = 0L,                                                      # 2
  seg.End       = as.integer(target_fai$length[match(used_t, target_fai$chr)]), # 3
  seg.name2     = used_t,                                                  # 4 (optional)
  seg.group     = target_label,                                            # 5 (optional)
  stringsAsFactors = FALSE
)
seg_query <- data.frame(
  seg.name      = paste0(Qpref, used_q),
  seg.Start     = 0L,
  seg.End       = as.integer(query_fai$length[match(used_q, query_fai$chr)]),
  seg.name2     = used_q,
  seg.group     = query_label,
  stringsAsFactors = FALSE
)

# Drop any NA-length segments (defensive)
seg_target <- seg_target[is.finite(seg_target$seg.End), , drop = FALSE]
seg_query  <- seg_query [is.finite(seg_query$seg.End),  , drop = FALSE]
stopifnot(nrow(seg_target) > 0, nrow(seg_query) > 0)

seg.frame <- rbind(seg_target, seg_query)
seg.order <- seg.frame$seg.name
db <- OmicCircos::segAnglePo(seg.frame, seg = seg.order)

# 3) Links (Kronos -> Svevo)
key_col <- if ("name_base" %in% names(merged)) {
  "name_base"
} else if ("name_query" %in% names(merged)) {
  "name_query"
} else {
  names(merged)[1]
}

link <- data.frame(
  seg1  = paste0(Tpref, merged$chr_target_std),
  po1   = as.integer(merged$midpoint_target),
  name1 = merged[[key_col]],
  seg2  = paste0(Qpref, merged$chr_query_std),
  po2   = as.integer(merged$midpoint_query),
  name2 = merged[[key_col]],
  stringsAsFactors = FALSE
)

# Keep only links with existing segments & finite positions
keep <- link$seg1 %in% seg_target$seg.name &
        link$seg2 %in% seg_query$seg.name &
        is.finite(link$po1) & is.finite(link$po2)
link <- link[keep, , drop = FALSE]
stopifnot(nrow(link) > 0)

# Clamp positions to segment ends
end_lookup <- setNames(seg.frame$seg.End, seg.frame$seg.name)
link$po1 <- pmin(pmax(link$po1, 0L), end_lookup[link$seg1])
link$po2 <- pmin(pmax(link$po2, 0L), end_lookup[link$seg2])

# 4) Plot
col_target <- "#CC5146"  # Kronos (dark salmon)
col_query  <- "#007C84"  # Svevo  (dark teal)
col_link   <- "#2B2B2B"

pdf_path <- file.path(out_dir, paste0(base_prefix, "_omicircos_links.pdf"))
pdf(pdf_path, width = 10, height = 10, useDingbats = FALSE)
par(mar = c(2,2,3,2))
plot(c(1,800), c(1,800), type = "n", axes = FALSE, xlab = "", ylab = "",
     main = paste0("PHAS cross-genome mapping: ", query_label, " ↔ ", target_label))

seg.cols <- c(
  rep(adjustcolor(col_target, 0.65), nrow(seg_target)),
  rep(adjustcolor(col_query,  0.65), nrow(seg_query))
)

OmicCircos::circos(R = 380, cir = db, type = "chr",
                   col = seg.cols, W = 22, scale = TRUE, print.chr.lab = TRUE)

OmicCircos::circos(R = 360, cir = db, W = 60, mapping = link, type = "link",
                   col = adjustcolor(col_link, 0.45), lwd = 0.9)

dev.off()
cat("Wrote: ", pdf_path, "\n")

###############################

## === Export Svevo→Kronos PHAS liftover table for Fig. B work ===
## Requires already in memory: merged, out_dir, base_prefix

stopifnot(exists("merged"))

# pick the best locus key we have
key_col <- if ("name_base" %in% names(merged)) "name_base" else if ("name_query" %in% names(merged)) "name_query" else names(merged)[1]

# keep only rows with both sides present (inner-join already did this),
# drop any UN/Unplaced if they linger
lift_out <- subset(merged,
                   !is.na(chr_query_std) & chr_query_std != "UN" &
                   !is.na(chr_target_std) & chr_target_std != "UN" &
                   is.finite(start_query) & is.finite(end_query) &
                   is.finite(start_target) & is.finite(end_target))

# select & rename columns for the PI
lift_out <- data.frame(
  PHAS_id          = lift_out[[key_col]],
  Svevo_chr        = lift_out$chr_query_std,
  Svevo_start      = as.integer(lift_out$start_query),
  Svevo_end        = as.integer(lift_out$end_query),
  Svevo_midpoint   = as.integer(lift_out$midpoint_query),
  Kronos_chr       = lift_out$chr_target_std,
  Kronos_start     = as.integer(lift_out$start_target),
  Kronos_end       = as.integer(lift_out$end_target),
  Kronos_midpoint  = as.integer(lift_out$midpoint_target),
  stringsAsFactors = FALSE
)

# (optional) sort by chromosome then start for readability
ord <- order(
  as.integer(sub("^([0-9]+).*", "\\1", lift_out$Kronos_chr)),
  sub("^[0-9]+", "", lift_out$Kronos_chr),
  lift_out$Kronos_start,
  na.last = TRUE
)
lift_out <- lift_out[ord, , drop = FALSE]

# write it
out_file <- file.path(out_dir, paste0(base_prefix, "_SvevoPHAS_lifted_to_Kronos.tsv"))
write.table(lift_out, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote liftover table for PI:\n  ", out_file, "\n")

# quick sanity prints
cat("Rows exported:", nrow(lift_out), "\n")
utils::head(lift_out, 3)
q()
