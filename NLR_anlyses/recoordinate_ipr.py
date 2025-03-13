import sys
from Bio import SeqIO

def process_orf(orf_output):
  #get original genomic coordinates
  orf_coordinates = {}
  for record in SeqIO.parse(orf_output, 'fasta'):
    seqid = record.id
    start, end = map(int, record.description.replace('[', '').replace(']', ' ').split()[1].split('-'))
    direction = record.description.split('(')[1].split(')')[0]
    orf_coordinates[seqid] = [start, end, direction]

  return orf_coordinates

def process_iprscan(iprscan_output, orf_coordinates, continue):
  #reformat interproscan results
  for line in open(iprscan_gff3, 'r'):
    items = line.split("\t")
    orf, start, stop, description = items[0], items[6], items[7], items[5]
    if not orf.startswith(chromosome):
      continue
      
    gff_fields = [''] * 9

    gff_fields[0] = orf.rsplit('_', 1)[0] #nlr loci
    gff_fields[1] = 'iprscan'
    gff_fields[2] = 'match'

    #adjust coordinates
    if orf_coordinates[orf][2] == "+":
      gff_fields[3] = str( int(start) * 3 - 3 + orf_coordinates[orf][0] + 1 )# orf 0 index
      gff_fields[4] = str( int(stop) * 3 - 3 + orf_coordinates[orf][0] + 1 )

    elif orf_coordinates[orf][2] == "-":
      gff_fields[3] = str(orf_coordinates[orf][1] - (int(stop) * 3 - 3))  # orf 0 index
      gff_fields[4] = str(orf_coordinates[orf][1] - (int(start) * 3 - 3))

    gff_fields[5] = '.'
    gff_fields[6] = orf_coordinates[orf][2]
    gff_fields[7] = '.'
    gff_fields[8] = f'Name={description};ID=match_{ct}'
    ct += 1

    print('\t'.join(gff_fields))

iprscan = sys.argv[1]
orfs_fa = sys.argv[2]
chromosome = sys.argv[3]

orf_coordinates = process_orf(orfs_fa)
process_iprscan(iprscan, orf_coordinates, chromosome)
