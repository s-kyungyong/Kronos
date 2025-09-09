
from Bio import SeqIO
import sys

fasta = sys.argv[1]
prefix = sys.argv[2]

seqID = { str(record.id):'' for record in SeqIO.parse(fasta, 'fasta') }
chromosome = { '1A': '1', '1B': '2', 'Chr1': '3',
               '2A': '4', '2B': '5', 'Chr2': '6',
               '3A': '7', '3B': '8', 'Chr3': '9',
               '4A': '10', '4B': '11', 'Chr4': '12',
               '5A': '13', '5B': '14', 'Chr5': '15',
               '6A': '16', '6B': '17', 'Chr6': '18',
               '7A': '19', '7B': '20', 'Chr7': '21' }

with open(f'{prefix}.bed', 'w') as o1, open(f'{prefix}.gff', 'w') as o2:
  for line in open(gff, 'r'):
      if line.startswith('#'):
          continue
  
      fields = line.split("\t")
      if not fields[2] == "mRNA":
          continue
  
      if fields[0] not in chromosome:
          continue
  
      tid = fields[-1].split('ID=')[1].split(';')[0].replace('transcript:', '')
  
      if tid in seqID:
        chrom = f'{prefix}{chromosome[fields[0]]}'
        o1.write(f'{chrom}\t{tid}\t{fields[3]}\t{fields[4]}')
        o2.write(f'{chrom}\t{fields[4]}\t{fields[4]}\t{tid}')
