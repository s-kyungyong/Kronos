from Bio import SeqIO
import subprocess
import random
import os
import sys

MIN_SEQ_LENGTH = 350
SEQID=99.5
num_gene=int(sys.argv[1])

# Constants
def parse_gff_attributes(attributes):
    """
    Parse GFF attributes and return them as a dictionary.
    """
    attr_dict = {}
    for attribute in attributes.split(';'):
        if '=' in attribute:
            key, value = attribute.split('=')
            attr_dict[key] = value
    return attr_dict

def select_braker_models(blast_out, braker_aa, braker_gff3):
    """
    Select Braker gene models based on BLAST results.
    """
    fasta_dict = SeqIO.to_dict(SeqIO.parse(braker_aa, 'fasta'))
    selected = {}

    stop_codon = {}
    for line in open(braker_gff3, 'r'):
      if line.split()[2] == 'stop_codon':
        attrs = parse_gff_attributes( line.split()[-1] )
        stop_codon[ attrs['Parent'] ] = ''

    with open(blast_out, "r") as blast_file:
        for line in blast_file:
            fields = line.split()
            q_start, h_start, qend, hend, qlen, hlen = fields[6], fields[8], fields[7], fields[9], fields[-2], fields[-1]
            perid = float( fields[2] )
            transcript = fields[0]

            if perid >= SEQID and q_start == h_start == "1" and qlen == hlen == qend == hend and int(qlen) >= MIN_SEQ_LENGTH:
                seq = fasta_dict.get(fields[0], None)
                seqlen = len(seq.seq)
                gene = transcript.split(".")[0]
                if seq and seq.seq.startswith('M') and fields[0] in stop_codon:
                    if transcript not in selected:
                      selected[ gene ] = [ transcript, seqlen ]
                    else:
                      if seqlen > selected[gene][1]:
                        selected[ gene ] = [ transript, seqlen ]

    return [selected[key][0] for key in selected]

def generate_gff_from_ids(input_gff, output_gff, selected_ids):
    """
    Generate a GFF file containing only the selected gene models.
    """
    target_ids = list(set(selected_ids))
    mRNA_coordinates = {}
    gids = [ t.split(".")[0] for t in target_ids ]

    with open(output_gff, 'w') as output_file, open(input_gff, "r") as gff_file:
        for line in gff_file:
            if line.startswith("#"):
              continue

            columns = line.split('\t')
            feature = columns[2]
            attributes = parse_gff_attributes(columns[-1])

            if feature == "gene":
              gname = attributes.get('ID', '')
              if gname in gids:
                 output_file.write(line)

            elif feature == "mRNA":
              gname = ".".join(attributes.get('ID', '').split(".")[:2])
              if gname in target_ids:
                mRNA_coordinates[gname] = [ int(columns[3]), int(columns[4]) ]
                output_file.write(line)

            elif feature == "exon" or feature == "CDS":
              gname = ".".join(attributes.get('ID', '').split(".")[:2])
              parent = attributes.get('Parent', '')

              if gname in target_ids:
                if columns[6] == "-":
                  exon_num = attributes.get('ID', '').split(".")[-1]
                  exon_num = exon_num.replace('exon', '').replace('CDS', '')
                  if exon_num == "1":
                   start = int(columns[3])
                   if start -3 == mRNA_coordinates[parent][0]:
                     columns[3] = str( start - 3 )

                elif columns[6] == '+':
                   end = int(columns[4])
                   if end + 3 == mRNA_coordinates[parent][1]:
                     columns[4] = str( end + 3)

                output_file.write("\t".join(columns))

    print(f'{output_gff} has {len(selected_ids)} genes')

# Define input and output file names
braker_aa = 'braker.aa'
braker_gff = 'braker.gff3'
blast_out = 'braker.blast.out'

augustus_gff = 'augustus.gff3'
snap_gff = 'snap.gff3'

# Select the Braker gene models for training
gene_models = select_braker_models('braker.blast.out', braker_aa, braker_gff)

print(f"Total number of selected gene models: {len(gene_models)}")

# Split the gene models for Augustus and SNAP
random.shuffle(gene_models)
augustus = gene_models[0:num_gene]
rest = gene_models[num_gene:]
random.shuffle(rest)
snap = rest[0:num_gene]

# Generate two GFF files
generate_gff_from_ids(braker_gff, augustus_gff, augustus)
generate_gff_from_ids(braker_gff, snap_gff, snap)
