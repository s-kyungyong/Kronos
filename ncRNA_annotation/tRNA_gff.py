#!/usr/bin/env python3
from collections import defaultdict
from intervaltree import Interval, IntervalTree
import sys
import re

# ------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------
PREFIX   = "TrturKRN"          # species + cultivar + assembly tag
VERSION  = "01"                # used in column-2 (source)
TYPE_PREFIX = "TR"             # maps biotype -> prefix in ID
ctr = defaultdict(int)

# ------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------
def parse_rfam_hits(rfam_file):
    """Return list of (chrom, strand, start, end) for every tRNA hit."""
    hits = []
    with open(rfam_file) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            _, name, accession, chrom, _, _, _, _, _, start, end, strand, *rest = ln.split()
            if name not in ("tRNA", "tRNA-Sec"):
                continue
            s, e = int(start), int(end)
            if s > e:
                s, e = e, s
            hits.append((chrom, strand, s, e))
    return hits


def build_rfam_trees(hits):
    trees = defaultdict(IntervalTree)
    for chrom, strand, s, e in hits:
        trees[(chrom, strand)].add(Interval(s, e + 1))  # 1-based inclusive
    return trees


def enough_overlap(qs, qe, intervals, thr=0.7):
    qlen = qe - qs + 1
    for iv in intervals:
        ov = max(0, min(qe + 1, iv.end) - max(qs, iv.begin))
        if ov / qlen >= thr:
            return True
    return False


def make_new_id(chrom):
    ctr[chrom] += 1
    return f"{PREFIX}{chrom}01{TYPE_PREFIX}2{ctr[chrom]:05d}"


def integrate(tRNAscan_file, rfam_trees, out_path):
    records = []                       # collect lines for later sorting

    with open(tRNAscan_file) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue

            chrom, num, start, end, isotype, anticodon, int_s, int_e, *rest = ln.split()
            if isotype == "Undet":
                continue

            s, e = int(start), int(end)
            strand = "+" if s < e else "-"
            if s > e: s, e = e, s

            int_s, int_e = int(int_s), int(int_e)
            if int_s and int_s > int_e:
                int_s, int_e = int_e, int_s

            pseudo = "pseudo" in ln.lower()
            biotype = "tRNA_pseudogene" if pseudo else "tRNA"

            # Rfam support?
            tree_key = (chrom, strand)
            rfam_sup = tree_key in rfam_trees and \
                       enough_overlap(s, e, rfam_trees[tree_key])

            # assign stable ID
            gene_id = f"{chrom}-{num}"

            # --- build attributes ---
            g_attrs = f"ID={gene_id};Name={gene_id};biotype={biotype}"
            t_attrs = [
                f"ID={gene_id}.1",
                f"Parent={gene_id}",
                f"biotype={biotype}",
                f"anticodon={anticodon}",
                f"isotype={isotype}",
                "source=tRNAscan-SE"
            ]
            if pseudo:
                t_attrs.append("pseudo=true")
            if rfam_sup:
                t_attrs.append("rfam_supported=true")
            if int_s:
                t_attrs.append(f"Note=predicted intron {int_s}-{int_e}")

            # stash tuple for later sort: (chrom, start, strand, line)
            records.append((chrom, s, strand,
                f"{chrom}\tKRNncRNAv1.0\tncRNA_gene\t{s}\t{e}\t.\t{strand}\t.\t{g_attrs}\n"))
            records.append((chrom, s, strand,
                f"{chrom}\tKRNncRNAv1.0\t{biotype}\t{s}\t{e}\t.\t{strand}\t.\t{';'.join(t_attrs)}\n"))

    # ------------------------------------------------------------------
    # sort: chrom (lex), start (int), strand
    # ------------------------------------------------------------------
    records.sort(key=lambda x: (x[0], x[1], x[2]))
    id_map = {}
    with open(out_path, "w") as out:
        out.write("##gff-version 3\n")
        for chrom, start, strand, line in records:
            if "ncRNA_gene" in line.split()[2]:
                old_id   = re.search(r"ID=([^;]+)", line).group(1)
                new_id = make_new_id(chrom)
                id_map[old_id] = new_id
            else:
                old_id = re.search(r"Parent=([^;]+)", line).group(1)
                new_id = id_map[old_id]

            out.write(line.replace(old_id, new_id))


# ------------------------------------------------------------------
# main driver
# ------------------------------------------------------------------
if __name__ == "__main__":
    rfam_file   = "All.Rfam.tblout"
    trna_file   = "tRNAscan-SE.out"
    output_file = "tRNAs.final.gff3"

    trees = build_rfam_trees(parse_rfam_hits(rfam_file))
    integrate(trna_file, trees, output_file)
    print(f"✔ Wrote {output_file}")
