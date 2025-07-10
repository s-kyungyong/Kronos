#!/usr/bin/env python3
from collections import defaultdict
from intervaltree import Interval, IntervalTree
import pandas as pd
import sys
import re

# ------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------
PREFIX   = "TrturKRN"          # species + cultivar + assembly tag
VERSION  = "01"                # used in column-2 (source)
TYPE_PREFIX = "LN"             # maps biotype -> prefix in ID
ctr = defaultdict(int)

# ------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------

def make_new_id(chrom):
    ctr[chrom] += 10
    return f"{PREFIX}{chrom}{VERSION}{TYPE_PREFIX}{ctr[chrom]:06d}"

def classify_combined(path, distance_cutoff=2000):
    """Return {lncRNA_gene: composite_class_label} with multiple components."""
    classes = {}

    with open(path) as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or ln.startswith("isBest"):
                continue

            cols = ln.split("\t")
            is_best      = cols[0] == "1"
            lnc_gene     = cols[1]
            direction    = cols[5]      # sense / antisense
            rel_type     = cols[6]      # intergenic / exonic / intronic
            distance     = int(cols[7])
            location     = cols[9]      # upstream / downstream

            if not is_best:
                continue

            if direction == "antisense" and rel_type == "intergenic" and distance <= distance_cutoff and location == "upstream":
                direction_label = "bidirectional"
            else:
                direction_label = direction

            if rel_type == "intergenic":
                if distance > distance_cutoff:
                    type_label = "distal"
                else:
                    type_label = location  # upstream or downstream
            else:
                type_label = rel_type  # exonic or intronic

            composite_label = f"{direction_label}_{type_label}"
            classes[lnc_gene] = composite_label

    return classes

def TE_labels(input_file):
    df = pd.read_csv(input_file, sep="\t")
    best = df.sort_values("Overlap", ascending=False).drop_duplicates("LncRNA_ID")
    return {row.LncRNA_ID: (row.TE_ID, row.Overlap) for _, row in best.iterrows()}

def integrate(input_file, output_file, classified, te_info):
    records = []                       # collect lines for later sorting
    gene_dict = {}
    exon_count = defaultdict(int)

    with open(input_file) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue

            fields = ln.strip().split("\t")
            chrom = fields[0]
            feature_type = fields[2]
            s, e = int(fields[3]), int(fields[4])
            strand = fields[6]
            attributes = fields[8]

            # --- build gene attributes ---
            if fields[2] == "transcript":
                transcript_id = re.search(r"ID=([^;]+)", attributes).group(1)
                gene_id = transcript_id.rsplit('.', 1)[0]

                g_attrs = [f"ID={gene_id}", f"Name={gene_id}", "biotype=lncRNA"]
                if gene_id in classified:
                    g_attrs.append(f"class={classified[gene_id]}")

                if gene_id in te_info:
                    g_attrs.append("te_overlap=true")
                    g_attrs.append(f"te_family={te_info[gene_id][0]}")
                    g_attrs.append(f"overlap={te_info[gene_id][1]}")

                if gene_id not in gene_dict:
                    gene_dict[gene_id] = [s, e, strand, chrom, g_attrs]

                else:
                    gene_dict[gene_id][0] = min(gene_dict[gene_id][0], s)
                    gene_dict[gene_id][1] = max(gene_dict[gene_id][1], e)

                t_attrs = [
                    f"ID={transcript_id}",
                    f"Parent={gene_id}"
                ]

                records.append((chrom, s, strand,
                    f"{chrom}\tKRNncRNAv1.0\talncRNA\t{s}\t{e}\t.\t{strand}\t.\t{';'.join(t_attrs)}\n"))

            if feature_type == "exon":
                transcript_id = re.search(r"Parent=([^;]+)", attributes).group(1)
                exon_count[transcript_id] += 1

                e_attrs = [
                    f"ID={transcript_id}.exon{exon_count[transcript_id]}",
                    f"Parent={transcript_id}",
                ]

                records.append((chrom, s, strand,
                    f"{chrom}\tKRNncRNAv1.0\texon\t{s}\t{e}\t.\t{strand}\t.\t{';'.join(e_attrs)}\n"))


    # Sort: chrom (lex), start (int), strand
    records.sort(key=lambda x: (x[0], x[1], x[2]))
    id_map = {}
    written = {}

    with open(output_file, "w") as out:
        out.write("##gff-version 3\n")
        for chrom, start, strand, line in records:
            if "alncRNA" in line:
                old_id = re.search(r"Parent=([^;]+)", line).group(1).strip()

                if not old_id in id_map:
                    new_id = make_new_id(chrom)
                    id_map[old_id] = new_id
                new_id = id_map[old_id]

                if old_id not in written:
                    written[old_id] = ''
                    s, e, st, c, g_attrs = gene_dict[old_id]
                    g_components = map(str, [ c, "KRNncRNAv1.0" , "ncRNA_gene", s, e, ".", st, ".", ";".join(g_attrs).replace(old_id, new_id) ])
                    out.write("\t".join(g_components) + "\n")

                out.write(line.replace(old_id, new_id).replace(f"Name={old_id}", f"Name={new_id}").replace("alncRNA", "lncRNA"))

            else:
                tr_id = re.search(r"Parent=([^;]+)", line).group(1).strip()
                old_id = tr_id.rsplit('.', 1)[0]
                new_id = id_map.get(old_id, old_id)

                if strand == "-":
                    exon_id = re.search(r"ID=([^;]+)", line).group(1).strip()
                    old_ct = exon_id.split("exon")[1]
                    ex_ct = exon_count[ tr_id ] - int(old_ct) + 1
                    out.write(line.replace(old_id, new_id).replace(f"exon{old_ct}", f"exon{ex_ct}"))

                else:
                    out.write(line.replace(old_id, new_id))

# ------------------------------------------------------------------
# main driver
# ------------------------------------------------------------------
if __name__ == "__main__":
    te_path = "TE_lncRNA_intersect.txt"
    class_path = "lncRNA_classes.txt"
    input_file = "mikado_stringtie_merged_final_lncRNA.gff3"
    output_file = "lncRNAs.final.gff3"

    classified = classify_combined(class_path)
    te_info = TE_labels(te_path)
    integrate(input_file, output_file, classified, te_info)

    print(f"✔ Wrote {output_file}")
