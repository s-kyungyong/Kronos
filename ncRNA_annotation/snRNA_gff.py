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
TYPE_PREFIX = "SN"             # maps biotype -> prefix in ID
reserved = "5"
ctr = defaultdict(int)

# ------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------

def make_new_id(chrom):
    ctr[chrom] += 10
    return f"{PREFIX}{chrom}01{TYPE_PREFIX}{reserved}{ctr[chrom]:05d}"

def select_rf(category_file):

    rfam = {}
    with open(category_file) as fh:
        for ln in fh:
            rfam[ ln.split()[0] ] = 'snRNA' if 'snoRNA' not in ln else 'snoRNA'
    return rfam

def integrate(input_file, rfam, out_path):
    records = []                       # collect lines for later sorting

    num = 0
    with open(input_file) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue

            num += 1
            fields = ln.strip().split()
            racc = fields[2]
            if racc in rfam:
                biotype = rfam[ racc ]
            else:
                continue

            chrom, s, e = fields[3], int(fields[9]), int(fields[10])
            if s > e: s, e = e, s

            match_length = int(fields[8]) - int(fields[7]) + 1
            cov = (e - s + 1 )/match_length
            product = fields[1] + ":" + " ".join(fields[28:])
            strand = fields[11]

            gene_id = f"{chrom}-{num}"
            # --- build attributes ---
            g_attrs = f"ID={gene_id};Name={gene_id};biotype={biotype}"
            t_attrs = [
                f"ID={gene_id}.1",
                f"Parent={gene_id}",
                f"biotype={biotype}",
                f"product={product}",
                f"Note=Rfam:{racc}",
                "source=Rfam"
            ]
            if cov < 0.7:
                 t_attrs.append("partial=true")

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
    category_file = 'snRNA_category.list'
    output_file = "snoRNAs.final.gff3"

    rfam = select_rf(category_file)
    integrate(rfam_file, rfam, output_file)
    print(f"✔ Wrote {output_file}")
