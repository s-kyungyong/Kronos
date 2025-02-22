from Bio import SeqIO

infile  = 'confident_TE.cons.fa.classified'
outfile = 'confident_TE.cons.fa.classified.filtered.fa'

with open(outfile, 'w') as o:
    for record in SeqIO.parse(infile, 'fasta'):
        if 'Unknown' not in record.id:
            seqid = record.id
            if seqid.count('/') == 0:
                seqid += '/unknown'

            o.write(f'>{seqid}\n{record.seq}\n')
