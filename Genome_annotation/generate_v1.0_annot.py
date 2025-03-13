from Bio import SeqIO

stop_codons = {'TAA', 'TAG', 'TGA'}
fasta_records = {}
geneCount, gname, exonCount, currentgeneCount, tr2gene, gene2tr = {}, {}, {}, {}, {}, {}

cds_fa   = 'Kronos.EVM.pasa.cds.fa'
gff3     = 'Kronos.EVM.pasa.gff3'
pasa_out = 'Kronos.EVM.pasa.against.all.pasa_orf.max10'
prot_out = 'Kronos.EVM.pasa.against.all.all_prot.max10'
description = 'Kronos.EVM.pep.descriptions'

def classify_cds(sequence):
    start, end = sequence[:3].upper(), sequence[-3:].upper()
    length = len(sequence) // 3

    #this may actually pick up CDS not divisible by 3 and contains stop
    #should have checked whether sequence.translate() startswith 'M' and ends with '*')
    return ("complete", length - 1) if start == "ATG" and end in stop_codons else ("partial", length - (1 if end in stop_codons else 0))

# Parse CDS fasta file
for record in SeqIO.parse(cds_fa, 'fasta'):
    cds_type, length = classify_cds(str(record.seq))
    fasta_records[record.id] = f'CDS_type={cds_type};start={record.seq[:3]};end={record.seq[-3:]};length={length};'

# Parse transcript evidence
done = set()
with open(pasa_out) as file:
    for line in file:
        cols = line.split()
        q, perid, q_start, q_end, q_len, h_start, h_end, h_len, eval_score = cols[0], cols[2], cols[6], cols[7], cols[-2], cols[8], cols[9], cols[-1], cols[10]
        if q not in done:
            t_type = 'complete' if q_len == h_len and q_start == h_start == '1' and q_end == h_end else 'incomplete'
            done.add(q)
        fasta_records[q] += f'transcript={t_type};transcript_match={q_start}-{q_end}_{perid}_{eval_score};'

# Ensure all records have transcript annotations
for key in fasta_records:
    if 'transcript_match' not in fasta_records[key]:
        fasta_records[key] += 'transcript=None;transcript_match=None;'

# Parse protein evidence
prot_evidence = {}
with open(prot_out) as file:
    for line in file:
        cols = line.split()
        q, h, perid, q_start, q_end, q_len, h_start, h_end, h_len, eval_score = cols[0], cols[1], cols[2], cols[6], cols[7], cols[-2], cols[8], cols[9], cols[-1], cols[10]
        qcov, hcov = round((int(q_end) - int(q_start) + 1) / float(q_len), 2), round((int(h_end) - int(h_start) + 1) / float(h_len), 2)
        p_type = 'complete' if float(perid) > 80.0 and qcov >= 0.97 and hcov >= 0.97 else 'incomplete'
        formatted = f'protein={p_type};protein_match={q_start}-{q_end}_{perid}_{eval_score}_{h}'
        if q not in prot_evidence or (p_type == 'complete' and '=complete' not in prot_evidence[q]):
            prot_evidence[q] = formatted

# Assign protein evidence to records
for key in fasta_records:
    fasta_records[key] += prot_evidence.get(key, 'protein=None;protein_match=None')

# store description
with open(description, 'w') as o:
    for key in fasta_records:
      desc = fasta_records[key]

      if 'protein=complete' in desc and 'CDS_type=complete' in desc:
        o.write(f'{key}\tversion=1.0;quality=high;{fasta_records[key]}\n')
      else:
        o.write(f'{key}\tversion=1.0;quality=low;{fasta_records[key]}\n')

# Store gene counts for each chromosome
for line in open(gff3):

    if line.startswith('#'):
        continue

    items = line.split()

    if items[2] == "gene":
      if items[0] not in geneCount:
        geneCount[ items[0] ] = 0
      geneCount[items[0]] += 1

    elif items[2] == "mRNA":
        parent = items[-1].strip().split('Parent=')[1].split(';')[0]
        rnaid  = items[-1].strip().split('ID=')[1].split(';')[0]
        tr2gene[rnaid] = parent
        gene2tr[parent] = rnaid

# Parse GFF file for gene renaming
for line in open(gff3):
    if line.startswith('#'):
        continue

    items = line.split()
    if not items:
        continue
    chr, attr = items[0], items[2]
    if attr == "gene":
        currentgeneCount[chr] = currentgeneCount.get(chr, 0) + 1
        raw_gNum = currentgeneCount[chr] * 10
        gNum = '0' * (6 - len(str(raw_gNum))) + str(raw_gNum)
        gname[items[-1].split("ID=")[1].split(';')[0]] = f'TrturKRN{chr}01G{gNum}'
    elif attr == "mRNA":
        exonCount[items[-1].split("ID=")[1].split(";")[0]] = 0
    elif attr == "exon":
        parent = items[-1].split("Parent=")[1].strip()
        exonCount[parent] += 1

# Write Name conversion file
with open('Name_conversion.list', 'w') as o:
    for item in gname:
        o.write(f'{item}\t{gname[item]}\n')

# Read descriptions
gdict = {line.split()[0]: line.split()[1] for line in open(description)}
hq    = { tr2gene[key]: gdict[key] for key in gdict if 'quality=high' in gdict[key] }

tr_counts = {}
tr_name   = {}
with open(f'{description}.w_converted_names', 'w') as o2, \
     open('Gene_confidnece.list', 'w') as o4, \
     open('Kronos.v1.0.high.gff3', 'w') as o, \
     open('Kronos.v1.0.low.gff3', 'w') as o3:
    for line in open(gff3):
        if line.startswith('#'):
            continue

        items = line.split()
        if not items:
            o.write('\n')
            continue
        chr, attr = items[0], items[2]
        items[1] = 'KRNv1.0'

        if attr == "gene":
            gid = items[-1].split("ID=")[1].split(';')[0]
            rid = gene2tr[gid]
            newName = gname[gid]
            items[-1] = 'ID=' + newName
            orgname = gid.replace('TU', 'model').replace('gene', 'model')
            quality = 'high' if gid in hq else 'low'
            o4.write(f'{orgname}\t{newName}\t{quality}\n')

        elif attr == "mRNA":
            gid = items[-1].split("Parent=")[1].split(';')[0]
            rid = items[-1].split('ID=')[1].split(';')[0]

            if gid not in tr_counts:
                tr_counts[gid] = 1

            newName = f'{gname[gid]}.{tr_counts[gid]}'
            tr_name[rid] = newName
            tr_counts[gid] += 1
            items[-1] = f'ID={newName};Parent={gname[gid]}'
            o2.write(f'{gid}\t{rid}\t{newName}\t{gdict[rid]}\n')

        else:
            rid = items[-1].split("Parent=")[1].split(';')[0]
            newName = tr_name[rid]
            items[-1] = f'Parent={newName}'

        (o if quality == 'high' else o3).write("\t".join(items) + "\n")
