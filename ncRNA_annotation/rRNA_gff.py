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
TYPE_PREFIX = "RR"             # maps biotype -> prefix in ID
reserved = "3"
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
            if not 'rRNA' in ln :
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
    ctr[chrom] += 10
    return f"{PREFIX}{chrom}01{TYPE_PREFIX}{reserved}{ctr[chrom]:05d}"


def integrate(input_file, rfam_trees, out_path):
    records = []                       # collect lines for later sorting

    num = 0
    with open(input_file) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue

            num += 1
            fields = ln.strip().split("\t")
            chrom, s, e = fields[0], int(fields[3]), int(fields[4])
            partial = "partial" in ln.lower()
            name = fields[-1].split('Name=')[1].split(';')[0]
            product = fields[-1].split('product=')[1].split(' (')[0]
            strand = fields[6]

            gene_id = f"{chrom}-{num}"
            biotype = "rRNA_fragment" if partial else "rRNA"

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
                f"product={product}",
                "source=barrnap"
            ]
            if partial:
                 t_attrs.append("partial=true")
            if rfam_sup:
                 t_attrs.append("rfam_supported=true")

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
    rfam_file   = "../../Rfam_analysis/All.Rfam.tblout"
    input_file   = "barrnap.gff"
    output_file = "rRNAs.final.gff3"

    trees = build_rfam_trees(parse_rfam_hits(rfam_file))
    integrate(input_file, trees, output_file)
    print(f"✔ Wrote {output_file}")
