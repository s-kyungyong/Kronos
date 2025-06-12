from Bio import SeqIO

refseq = 'NCBI_refseq.aa.fa'
dmnd_out = 'Kronos.v2.1.against.NCBI.refseq.dmnd.out'
eggnog_out = 'Kronos.v2.1.pep.eggnog.tsv'
output = 'Kronos.v2.1.description'

seqdict = SeqIO.to_dict(SeqIO.parse(refseq, 'fasta'))
description = {}

for line in open(dmnd_out, 'r'):
    fields = line.split()
    query, hit = fields[0], fields[1]
    gid = query.split('.')[0]

    # Corrected tuple unpacking
    qstart, qend, qlen, hstart, hend, hlen = map(int, (fields[6], fields[7], fields[-2], fields[8], fields[9], fields[-1]))
    perid, evalue = map(float, (fields[2], fields[10]))


    qcov = (qend - qstart + 1)/qlen
    hcov = (hend - hstart + 1)/hlen

    if qcov > 0.90 and hcov > 0.90 and perid > 0.90 and evalue < 1e-10:
        if query not in description:
            annot = seqdict[hit].description.split(' ', 1)[1].split('[')[0]

            #remove reelavant isoform tag
            if 'isoform' in annot:
                annot = annot.rsplit(' isoform', 1)[0]

            #remove irrelavant locus tacg
            if 'uncharacterized protein' in annot:
                annot = 'uncharacterized protein'

            description[gid] = annot

print(f'{len(description.keys())} annotations were transferred from the NCBI RefSeq')
ct = 0
for line in open(eggnog_out, 'r'):
    if line.startswith('#'):
        continue

    fields = line.split('\t')
    query, annot = fields[0], fields[7]
    gid = query.split('.')[0]

    if query in description:
        continue

    if gid not in description or description[gid] == '-':
        description[gid] = annot
        if annot != '-':
            ct += 1
print(f'{ct} annotations were transferred from the Eggnog output')

#replace not useful terms from eggnog
replace = [ '-' , 'source UniProtKB', 'protein Sb04g004280 source' ]

with open(output, 'w') as o:
    for record in sorted(description.keys()):
        annot = 'Hypothetical protein' if description[record] in replace else description[record]
        o.write(f'{record}\t{annot}\n')
