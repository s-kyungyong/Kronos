from Bio import SeqIO

input = "Putative_NLRs_in_Triticum_IsoSeq.aa.complete.fa"
output = "Putative_NLRs_in_Triticum_IsoSeq.aa.complete.hc.aa.fa"
dmnd  = "Putative_NLRs_in_Triticum_IsoSeq.aa.complete.against.dmnd.out"

with open(output, 'w') as o:
    complete = {}

    for line in open(dmnd, 'r'):
        fields = line.split("\t")
        qlen, hlen = int(fields[-2]), int(fields[-1])

        if qlen > hlen * 0.9:
            complete[fields[0]] = ''

    for record in SeqIO.parse(input, 'fasta'):
        if record.id in complete:
            SeqIO.write(record, o, 'fasta')
