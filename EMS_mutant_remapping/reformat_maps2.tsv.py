import sys

input_file = sys.argv[1]
output_file = input_file.replace('.tsv', '.reformatted.tsv')

accessions = "accessions.list"
suffix     = ".sorted.rmdup.bam"

accession2kr = {}
used         = {}
for line in open(accessions, 'r'):
  accession2kr[ line.split()[1] + suffix ] = line.split()[0]

with open(output_file, 'w') as o:
  for i, line in enumerate(open(input_file, 'r')):
    if i == 0:
      o.write(line)
      continue

    if 'Chrom' in line:
        continue

    items = line.split()
    items[6] = accession2kr[ line.split()[6] ]
    used[ line.split()[6] ] = ''
    if '_' in items[0]:
      chr, pos = items[0].split("_")
      items[0] = chr
      items[1] = str(int(pos) + int(items[1]))

    o.write("\t".join(items) + "\n")

for accession in accession2kr:
    if accession not in used:
        print(accession)
