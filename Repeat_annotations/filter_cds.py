from Bio import SeqIO

infile = 'trep-db_complete_Rel-19.triticum.filtered.against.Kronosv2.dmnd.out'
infa  = 'Kronos.v2.0.cds.fa'
outfa = 'Kronos.v2.0.cds.filtered.fa'


discard = {}
for line in open(infile, 'r'):
    fields = line.split()

    hstart, hend, hlen = map(int, [ fields[8], fields[9], fields[-1] ] )
    hcov = (hend - hstart + 1)/float(hlen)

    if hcov > 0.4 and float(fields[10]) < 1e-4:
        discard[ fields[1].split('.')[0] ] = ''

print(f'{len(discard.keys())} detected')
with open(outfa, 'w') as o:
    for record in SeqIO.parse(infa, 'fasta'):
        if record.id.split('.')[0] in discard:
            continue
        SeqIO.write(record, o, 'fasta')
