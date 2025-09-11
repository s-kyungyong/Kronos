from Bio import SeqIO
import sys 

input = sys.argv[1]
output = sys.argv[2]

with open(output, 'w') as o:
    for record in SeqIO.parse(input, 'fasta'):

        description = record.description
        start = description.split('start:')[1].split()[0]
        end   = description.split('stop:')[1].rstrip()

        if end in ['TAG', 'TGA', 'TAA'] and start == 'ATG':
            o.write(f'>{record.id.split()[0]}\n{record.seq}\n')
