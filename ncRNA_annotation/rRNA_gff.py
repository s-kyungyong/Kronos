#!/usr/bin/env python3

import sys
from collections import defaultdict
from intervaltree import Interval, IntervalTree

def parse_rfam_hits(rfam_file):
    """
    Parse Rfam hits and yield only tRNA-related hits.
    """
    hits = []
    with open(rfam_file) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            fields = ln.strip().split()
            _, name, accession, chrom, _, _, _, _, _, start, end, strand, *rest = fields
            if 'rRNA' not in name:
                continue
            start, end = int(start), int(end)
            if start > end:
                start, end = end, start
            hits.append((chrom, strand, start, end))
    return hits

def build_rfam_trees(hits):
    """
    Build an IntervalTree for each (chrom, strand).
    """
    trees = defaultdict(IntervalTree)
    for chrom, strand, start, end in hits:
        trees[(chrom, strand)].add(Interval(start, end + 1))  # 1-based inclusive
    return trees

def overlaps_enough(query_start, query_end, intervals, threshold=0.7):
    """
    Check if the query interval overlaps ≥ threshold with any interval in the tree.
    """
    query_len = query_end - query_start + 1
    for iv in intervals:
        ov_start = max(query_start, iv.begin)
        ov_end = min(query_end + 1, iv.end)
        overlap = max(0, ov_end - ov_start)
        if (overlap / query_len) >= threshold:
            return True
    return False

def integrate_rRNA_predictions(rRNA_gff, rfam_trees, output_file):
    """
    Read barrnap results, annotate with Rfam support, and write GFF3 output.
    """
    num = 0
    with open(output_file, 'w') as out:
        out.write("##gff-version 3\n")
        with open(rRNA_gff) as f:
            for ln in f:
                if ln.startswith("#") or not ln.strip():
                    continue
                num += 1
                fields = ln.strip().split("\t")
                chrom, start, end = fields[0], int(fields[3]), int(fields[4])
                partial = "partial" in ln.lower()
                name = fields[-1].split('Name=')[1].split(';')[0]
                product = fields[-1].split('product=')[1].split(' (')[0]
                strand = fields[6]

                tree_key = (chrom, strand)
                rfam_supported = False
                if tree_key in rfam_trees:
                    overlapping = rfam_trees[tree_key].overlap(start, end + 1)
                    rfam_supported = overlaps_enough(start, end, overlapping)

                gene_id = f"{chrom}-{num}-rRNA"
                transcript_type = "rRNA_fragment" if partial else "rRNA"

                # Gene line
                gene_attrs = [f"ID={gene_id}",
                    f"Name={name}",
                    f"biotype=rRNA",
                    f"product={product}",
                    "logic_name=barrnap"
                  ]
                if partial:
                      gene_attrs.append("partial=true")
                if rfam_supported:
                      gene_attrs.append("rfam_supported=true")

                attr_str = ";".join(gene_attrs)
                out.write(f"{chrom}\tKRN_ncRNA\tncRNA_gene\t{start}\t{end}\t.\t{strand}\t.\t{attr_str}\n")

                # Transcript line
                transcript_attrs = [
                    f"ID={gene_id}.1",
                    f"Parent={gene_id}",
                  ]

                attr_str = ";".join(transcript_attrs)
                out.write(f"{chrom}\tKRN_ncRNA\t{transcript_type}\t{start}\t{end}\t.\t{strand}\t.\t{attr_str}\n")

                # Exon line
                exon_attrs = [
                    f"ID={gene_id}.1.exon1",
                    f"Parent={gene_id}.1",
                ]

                attr_str = ";".join(exon_attrs)
                out.write(f"{chrom}\tKRN_ncRNA\texon\t{start}\t{end}\t.\t{strand}\t.\t{attr_str}\n")


def main():
    rfam_file = "All.Rfam.tblout"
    rRNA_file = "barrnap.gff"
    output_file = "rRNAs.final.gff3"

    rfam_hits = parse_rfam_hits(rfam_file)
    rfam_trees = build_rfam_trees(rfam_hits)
    integrate_rRNA_predictions(rRNA_file, rfam_trees, output_file)

if __name__ == "__main__":
    main()
