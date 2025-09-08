import argparse
import sys
from Bio import SeqIO
from intervaltree import Interval, IntervalTree

def merge_intervals(raw_coordinate_dict, flanking):
    """Merge overlapping intervals with flanking regions."""
    merged_coordinate_dict = {}

    for chromosome, intervals in raw_coordinate_dict.items():
        # Create an interval tree and add existing intervals with flanking regions
        tree = IntervalTree(Interval(start - flanking, end + flanking) for start, end in intervals)

        # Merge overlapping intervals
        tree.merge_overlaps()

        # Extract and sort merged intervals
        merged_intervals = sorted((interval.begin, interval.end) for interval in tree)
        merged_coordinate_dict[chromosome] = merged_intervals

    return merged_coordinate_dict

def read_hmmout(hmmout_file):
    """Read NLR coordinates from hmmsearch results."""
    coordinates = {}


    with open(hmmout_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue

            parts = line.split('\t')
            chromosome = parts[0].split('_')[0]
            start, end = map(int, parts[-1].split('[')[1].split(']')[0].split('-'))

            if chromosome not in coordinates:
                coordinates[chromosome] = []
            coordinates[chromosome].append((start, end))

    return coordinates

def read_nlrannot(nlrannot_file):
    """Read NLR coordinates from NLR-Annotator CSV results."""
    coordinates = {}

    with open(nlrannot_file, 'r') as file:
        next(file)  # Skip header
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            chromosome = parts[0]
            start, end = map(int, parts[3:5])

            if chromosome not in coordinates:
                coordinates[chromosome] = []
            coordinates[chromosome].append((start, end))

    return coordinates

def load_existing_contigs(contig_file):
    """Load existing contig intervals."""
    contig_intervals = {}

    with open(contig_file, "r") as file:
        for line in file:
            if line.startswith(">"):
                parts = line.strip()[1:].split("_")
                if len(parts) == 2:
                    chromosome, positions = parts
                    start, end = map(int, positions.split("-"))
                    if chromosome not in contig_intervals:
                        contig_intervals[chromosome] = IntervalTree()
                    contig_intervals[chromosome].add(Interval(start, end))

    return contig_intervals

def remove_overlapping_intervals(merged_intervals, contig_intervals):
    """Remove intervals that overlap with existing contig fragments."""
    filtered_intervals = {}

    for chromosome, intervals in merged_intervals.items():
        if chromosome not in contig_intervals:
            filtered_intervals[chromosome] = intervals
            continue

        filtered_intervals[chromosome] = []
        for start, end in intervals:
            # Check for overlap
            if not contig_intervals[chromosome].overlap(start, end):
                filtered_intervals[chromosome].append((start, end))
            else:
                print(f'Overlap detected: {chromosome} {start} {end}')

    return filtered_intervals

def crop_genome(intervals, genome_fa, output_name):
    """Extract genome regions that do not overlap with existing contigs."""
    with open(output_name, 'w') as output_file:
        for record in SeqIO.parse(genome_fa, 'fasta'):
            chromosome = str(record.id)
            chrLen = len(record.seq)
            if chromosome in intervals:
                for start, end in intervals[chromosome]:
                    start = max(1, start)  # Ensure start is at least 1
                    end = min(chrLen, end)  # Ensure end does not exceed chromosome length
                    output_file.write(f'>{chromosome}_{start}-{end}\n{record.seq[start-1:end]}\n')

def main():
    parser = argparse.ArgumentParser(description="Process NLR loci from HMM and NLR-Annotator outputs.")
    parser.add_argument("--hmm", help="Path to hmmsearch output file", required=False)
    parser.add_argument("--nlrannot", help="Path to NLR-Annotator CSV file", required=False)
    parser.add_argument("--genome", help="Path to genome FASTA file", required=True)
    parser.add_argument("--contigs", help="Path to existing NLR contigs file", required=False)
    parser.add_argument("--flanking", type=int, default=15000, help="Flanking region size (default: 15000 bp)")

    args = parser.parse_args()

    if not args.hmm and not args.nlrannot:
        print("Error: You must provide at least one input file (HMM or NLR-Annotator).", file=sys.stderr)
        sys.exit(1)

    raw_intervals = {}

    # Process HMM output if provided
    if args.hmm:
        hmm_intervals = read_hmmout(args.hmm)
        for chrom, intervals in hmm_intervals.items():
            if chrom not in raw_intervals:
                raw_intervals[chrom] = []
            raw_intervals[chrom].extend(intervals)

    # Process NLR-Annotator output if provided
    if args.nlrannot:
        nlrannot_intervals = read_nlrannot(args.nlrannot)
        for chrom, intervals in nlrannot_intervals.items():
            if chrom not in raw_intervals:
                raw_intervals[chrom] = []
            raw_intervals[chrom].extend(intervals)

    # Merge overlapping intervals
    merged_intervals = merge_intervals(raw_intervals, args.flanking)

    # Load existing contig regions if provided
    if args.contigs:
        existing_contigs = load_existing_contigs(args.contigs)
        non_overlapping_intervals = remove_overlapping_intervals(merged_intervals, existing_contigs)
    else:
        non_overlapping_intervals = merged_intervals

    # Output file name
    output_name = f"NLR_loci.fa"

    # Crop genome and save non-overlapping sequences
    crop_genome(non_overlapping_intervals, args.genome, output_name)

    print(f"Filtered extra NLR loci saved to {output_name}")

if __name__ == "__main__":
    main()
