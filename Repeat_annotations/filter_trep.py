infile  = 'trep-db_complete_Rel-19.fasta'
outfile = 'trep-db_complete_Rel-19.triticum.filtered.fa'

with open(outfile, 'w') as o:
    for record in SeqIO.parse(infile, 'fasta'):
        seqDesc = record.description
        elements = seqDesc.split(';')
        repTypes = [ x.strip() for x in elements[1].split(',') ]

        if 'Triticum' in seqDesc and repTypes[0] != "unknown":
            if 'consensus' in seqDesc or 'complete' in seqDesc:
              subclass = repTypes[1].split('(')[1][:-1] if '(' in repTypes[1] else repTypes[1]
              subclass = subclass if subclass != 'Helitron' else 'DNA'
              superfamily = repTypes[2]

              seqid = f'{record.id}#{subclass}/{superfamily}'
              o.write(f'>{seqid}\n{record.seq}\n')
