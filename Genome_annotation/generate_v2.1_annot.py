import sys
from intervaltree import Interval, IntervalTree

def check_overlap(intervals, chromosome, start, end, gene):
    """
    Check if the given interval overlaps with existing intervals on the chromosome,
    excluding the same gene.
    """
    if chromosome in intervals:
        overlapping_intervals = intervals[chromosome].overlap(start, end)
        filtered = [interval for interval in overlapping_intervals]
        if filtered:
            return True, filtered

    return False, []

def get_geneID(transcriptID):
    geneID = transcriptID.rsplit('.', 1)[0]
    return geneID

def load_coordinates(gff3_file):
    """
    Load high-confidence gene coordinates from the GFF3 file and create a dictionary
    mapping genes to their corresponding transcript IDs.
    """
    plus = {}
    minus = {}
    cds_entries = []
    gene_transcripts_dict = {}

    with open(gff3_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split()
            chromosome = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            direction = fields[6]

            if fields[2] == "CDS":
                transcript = fields[8].split('Parent=')[1].split(';')[0].rstrip()
                gene = get_geneID(transcript)
                cds_entry = (chromosome, start, end, direction, gene, transcript)
                cds_entries.append(cds_entry)

                if gene not in gene_transcripts_dict:
                    gene_transcripts_dict[gene] = []

                gene_transcripts_dict[gene].append(transcript)

                if direction == "+":
                    plus.setdefault(chromosome, []).append((start, end, transcript))
                elif direction == "-":
                    minus.setdefault(chromosome, []).append((start, end, transcript))

    return plus, minus, cds_entries, gene_transcripts_dict


def build_interval_trees(coordinates):
    """
    Convert high-confidence gene coordinates into interval trees per chromosome.
    """
    intervals = {}
    for chromosome, coords in coordinates.items():
        tree = IntervalTree()
        for start, end, gene in coords:
            if start != end:
                tree.add(Interval(start, end, gene))
        intervals[chromosome] = tree

    return intervals


def process_gff3(v2_annot, nlr_annot, nlr_confidence, output_gff3):
    """
    Process the GFF3 file and identify disqualified CDS entries.
    """
    plus, minus, cds_entries, gene_transcripts_dict = load_coordinates(v2_annot)
    plus_intervals = build_interval_trees(plus)
    minus_intervals = build_interval_trees(minus)

    nlr_confidence = {}
    to_be_removed = {}
    to_be_added   = {}
    for line in open(nlr_confidence, 'r'):
        fields = line.split()
        nlr_confidence[ fields[0] ] = fields[2]
        if fields[2] in ["High", "Medium"]:
            to_be_added[ fields[0] ] = ''

    _, _, nlr_entries, _ = load_coordinates(nlr_annot)

    for chromosome, start, end, direction, gene, transcript in nlr_entries:
        if direction == '+':
            overlap, overlapping_genes = check_overlap(plus_intervals, chromosome, start, end, gene)
            if overlap:
                for s, e, tr in overlapping_genes:
                    to_be_removed[tr.split('.')[0]] = ''

        else:
            overlap, overlapping_genes = check_overlap(minus_intervals, chromosome, start, end, gene)
            if overlap:
                for s, e, tr in overlapping_genes:
                    to_be_removed[tr.split('.')[0]] = ''


    print(f'{len(to_be_added.keys())} will be added')
    print(f'{len(to_be_removed.keys())} will be removed')

    removed = 0
    added   = 0
    with open(output_gff, 'w') as o:
        for line in open(v2_annot, 'r'):
            fields = line.strip().split()
            fields[1] = "KRNv2.1"

            if fields[-1].split('ID=')[1].split('.')[0].strip() in to_be_removed:
                if fields[2] == "gene":
                    removed += 1
                continue
            o.write('\t'.join(fields) + '\n')

        for line in open(nlr_annot, 'r'):
            fields = line.strip().split('\t')
            if fields[-1].split('ID=')[1].split('.')[0].split(';')[0].strip() in to_be_added:
                o.write('\t'.join(fields) + '\n')
                if fields[2] == "gene":
                    added += 1

    print(f'{added} was added')
    print(f'{removed} was removed')

if __name__ == "__main__":

    v2_annot  = 'Kronos.v2.0.gff3'
    v21_annot = "Kronos.v2.1.initial.gff3"
    nlr_confidence = 'Kronos_all.NLRs.final.conf.list'
    nlr_annot = 'Kronos_all.NLRs.final.gff3'
    process_gff3(v2_annot, nlr_annot, nlr_confidence, v21_annot)
