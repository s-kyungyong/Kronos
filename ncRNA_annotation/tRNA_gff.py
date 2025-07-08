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
            if name not in ['tRNA', 'tRNA-Sec']:
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

def integrate_tRNA_predictions(tRNAscan_file, rfam_trees, output_file):
    """
    Read tRNAscan results, annotate with Rfam support, and write GFF3 output.
    """
    with open(output_file, 'w') as out:
        out.write("##gff-version 3\n")
        with open(tRNAscan_file) as f:
            for ln in f:
                if ln.startswith("#") or not ln.strip():
                    continue
                chrom, num, start, end, anticodon, isotype, int_start, int_end, *rest = ln.strip().split()
                start, end = int(start), int(end)
                int_start, int_end = int(int_start), int(int_end)


                if isotype == 'Undet':
                    print(f'Ignoring {isotype}')
                    continue

                if start > end:
                    strand = "-"
                    start, end = end, start
                else:
                    strand = "+"  # adjust if available in input

                if int_start != 0 and int_start > int_end:
                    int_start, int_end = int_end, int_start

                pseudo = "pseudo" in ln.lower()

                tree_key = (chrom, strand)
                rfam_supported = False
                if tree_key in rfam_trees:
                    overlapping = rfam_trees[tree_key].overlap(start, end + 1)
                    rfam_supported = overlaps_enough(start, end, overlapping)

                gene_id = f"{chrom}-{num}-tRNA"
                transcript_type = "tRNA_pseudogene" if pseudo else "tRNA"
                # Gene line
                gene_attrs = [f"ID={gene_id}",
                    f"Name={gene_id}",
                    f"biotype={transcript_type}",
                    f"anticodon={anticodon}",
                    f"isotype={isotype}",
                    "logic_name=trnascan_gene"
                  ]
                if pseudo:
                      gene_attrs.append("pseudo=true")
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

                if int_start == 0:
                    out.write(f"{chrom}\tKRN_ncRNA\texon\t{start}\t{end}\t.\t{strand}\t.\t{attr_str}\n")
                else:
                    if strand == "+":
                      out.write(f"{chrom}\tKRN_ncRNA\texon\t{start}\t{int_start-1}\t.\t{strand}\t.\t{attr_str}\n")
                      out.write(f"{chrom}\tKRN_ncRNA\texon\t{int_end+1}\t{end}\t.\t{strand}\t.\t{attr_str.replace('.exon1', '.exon2')}\n")
                    else:
                      out.write(f"{chrom}\tKRN_ncRNA\texon\t{start}\t{int_start-1}\t.\t{strand}\t.\t{attr_str.replace('.exon1', '.exon.2')}\n")
                      out.write(f"{chrom}\tKRN_ncRNA\texon\t{int_end+1}\t{end}\t.\t{strand}\t.\t{attr_str}\n")


def main():
    rfam_file = "../../Rfam_analysis/All.Rfam.tblout"
    tRNA_file = "tRNAscan-SE.out"
    output_file = "tRNAs.final.gff3"

    rfam_hits = parse_rfam_hits(rfam_file)
    rfam_trees = build_rfam_trees(rfam_hits)
    integrate_tRNA_predictions(tRNA_file, rfam_trees, output_file)

if __name__ == "__main__":
    main()
