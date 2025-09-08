import re
import sys

# Input files
gff_file = sys.argv[1]
scaffold_file = "coordinates.list"
output_file = gff_file.replace('gff', 'recoordinated.gff')

# Load scaffold information into a dictionary
def load_scaffolds(scaffold_file):
    scaffolds = {}
    with open(scaffold_file, 'r') as f:
        for line in f:
            chrom, coords = line.strip().split('_', 1)
            start, end = map(int, coords.split('-'))
            scaffolds[chrom] = scaffolds.get(chrom, []) + [(start, end)]
    return scaffolds

# Adjust coordinates based on scaffold
def adjust_coordinates(scaffolds, chrom, start, end):
    if chrom not in scaffolds:
        return start, end, '', False  # Return unchanged if chromosome not found

    for scaffold_start, scaffold_end in scaffolds[chrom]:
        if scaffold_start <= start <= scaffold_end and scaffold_start <= end <= scaffold_end:
            new_start = start - scaffold_start + 1
            new_end = end - scaffold_start + 1
            new_chromosome = f"{chrom}_{scaffold_start}-{scaffold_end}"
            return new_start, new_end, new_chromosome, True
    return start, end, '', False  # Return unchanged if no matching scaffold found

# Process the GFF file
def recoordinate_gff(gff_file, scaffolds, output_file):
    candidates = {}
    with open(gff_file, 'r') as gff, open(output_file, 'w') as out:
        for line in gff:
            if line.startswith("#") or not line.strip() :
                continue

            fields = line.strip().split()
            chrom, source, feature, start, end, score, strand, phase, attributes = fields[:9]
            start, end = int(start), int(end)

            if fields[2] == "mRNA":
                gID = attributes.split('ID=')[1].rstrip().split(';')[0]
                new_start, new_end, new_chromosome, adjusted = adjust_coordinates(scaffolds, chrom, start, end)

                if adjusted:
                    candidates[gID] = [ new_start, new_end, new_chromosome, adjusted ]

            elif fields[2] == "exon" or fields[2] == "CDS":
                gID = attributes.split('Parent=')[1].rstrip().split(';')[0]

            else:
                continue

            if gID in candidates:
                new_start, new_end, new_chromosome, adjusted = adjust_coordinates(scaffolds, chrom, start, end)
                fields[3], fields[4] = str(new_start), str(new_end)
                fields[0] = new_chromosome
                fields[1] = "v1Annot" #fields[1]
                out.write('\t'.join(fields[:9]) + '\n')

# Main function
def main():
    scaffolds = load_scaffolds(scaffold_file)
    recoordinate_gff(gff_file, scaffolds, output_file)
    print(f"Re-coordinated GFF file saved to {output_file}")

if __name__ == "__main__":
    main()
