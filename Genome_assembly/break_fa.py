from Bio import SeqIO
import sys

genome = sys.argv[1]
broken = genome.replace('.fa', '.broken.fa')
window = 100000000 #100 Mb

with open(broken, 'w') as o:
  for record in SeqIO.parse(genome, 'fasta'):
    seqID = record.id
    sequence = record.seq
    seqLen = len(sequence)

    for k in range(seqLen//window + 1):
      # create 25 Mb overlaping regions on the left arm
      left_ = k * window - window//4
      left  = left_ if left_ > 0 else 0
      right = (k + 1) * window
      o.write(f'>{seqID}_{k}\n{sequence[ left: right ]}\n')
