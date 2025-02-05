import sys

# Ensure input file is provided
if len(sys.argv) < 2:
    sys.exit("Usage: python script.py <input_vcf.gz>")

input_file = sys.argv[1]
accessions_file = "MAPS_groups.list"

# Read accession mappings into a dictionary
accession2kr = {}
with open(accessions_file, 'r') as acc_file:
    for line in acc_file:
        parts = line.strip().split()
        accession2kr[parts[3]] = parts[1]

# Extract sample name (SRR ID) from filename
srr = input_file.split('.')[0]

# Ensure SRR ID is found in the mapping
if srr not in accession2kr:
    sys.exit(f"Error: {srr} not found in {accessions_file}")

# Update output filename using the accession mapping
output_file = output_file.replace(srr, accession2kr[srr])

# Define chromosome sizes only once
CHROMOSOME_SIZES = """\
##contig=<ID=1A,length=600443981>
##contig=<ID=1B,length=708842986>
##contig=<ID=2A,length=795820389>
##contig=<ID=2B,length=828541533>
##contig=<ID=3A,length=759128228>
##contig=<ID=3B,length=864152387>
##contig=<ID=4A,length=767865717>
##contig=<ID=4B,length=699696956>
##contig=<ID=5A,length=720280059>
##contig=<ID=5B,length=731153026>
##contig=<ID=6A,length=624303373>
##contig=<ID=6B,length=733599645>
##contig=<ID=7A,length=753476766>
##contig=<ID=7B,length=766026795>
##contig=<ID=Un,length=211250944>
"""

# Process the VCF file and reformat as needed
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    chromosome_size_written = False

    for line in infile:
        line = line.strip()  # Remove trailing whitespace/newlines

        # Insert chromosome sizes if encountering the first contig header
        if line.startswith('##contig='):
            if not chromosome_size_written:
                outfile.write(CHROMOSOME_SIZES)
                chromosome_size_written = True
            continue  # Skip the original contig lines

        # Write header lines as is
        if line.startswith('#'):
            outfile.write(line + '\n')
            continue

        # Process variant lines
        items = line.split('\t')

        # Handle chromosome renaming for "chr_pos" format
        if '_' in items[0]:
            chromosome, pos = items[0].split("_", 1)
            items[0] = chromosome
           
