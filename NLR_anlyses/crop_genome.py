from Bio import SeqIO
from intervaltree import Interval, IntervalTree
import sys

def merge_intervals(raw_coordinate_dict, flanking):
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

def crop_genome(intervals, genome_fa, output_name):
    with open(output_name, 'w') as output_file:
        for record in SeqIO.parse(genome_fa, 'fasta'):
            chromosome = str(record.id)
            chrLen     = len(record.seq)
            if chromosome in intervals:
                for start, end in intervals[chromosome]:
                    start = start if start > 0 else 1
                    end   = end if end < chrLen else chrLen
                    output_file.write(f'>{chromosome}_{start}-{end}\n{record.seq[start-1:end]}\n')

def main():
    if len(sys.argv) != 3:
        print("Usage: script.py <hmmout_file> <genome_fasta>")
        sys.exit(1)

    hmmfile = sys.argv[1]
    genome = sys.argv[2]
    output_name = f"{genome.split('/')[-1]}.NLR_loci.fa"
    flanking = 15000

    raw_intervals = read_hmmout(hmmfile)
    merged_intervals = merge_intervals(raw_intervals, flanking)
    crop_genome(merged_intervals, genome, output_name)

if __name__ == "__main__":
    main()
