from Bio import SeqIO

# Process low masking region
mask = {}
masked_seqs = 0
previous_seq_id = None

with open("Triticum_kronos.proteins.segmakser.out", "r") as mask_file:
    for line in mask_file:
        if line.startswith(">"):
            if previous_seq_id is not None:
                mask[previous_seq_id] = masked_seqs
            masked_seqs = 0
            previous_seq_id = line.split()[1].replace(">", "")
        else:
            start, end = map(int, line.split()[:2])
            masked_seqs += end - start + 1

# Store last sequence if it exists
if previous_seq_id is not None:
    mask[previous_seq_id] = masked_seqs

# Check start and stop codons
stop_codons = {"TAA", "TAG", "TGA"}
aaseq = {}

with open("Triticum_kronos.cds-transcripts.fa", "r") as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        if seq.startswith("ATG") and seq[-3:] in stop_codons:
            gene_id = str(record.id).rsplit("-", 1)[0]
            aaseq[gene_id] = len(seq) // 3  # Store length in codons

# Identify genes with low complexity regions
remove = {
    seq_id: ""
    for seq_id, masked_length in mask.items()
    if seq_id in aaseq and (masked_length / float(aaseq[seq_id])) > 0.70
}

# Function to extract gene ID from a GFF3 line
def get_gene_id(gff_line):
    fields = gff_line.strip().split("\t")
    attributes = fields[-1]

    if fields[2] == "gene":
        return attributes.split("ID=")[-1].split(";")[0]
    return attributes.split("ID=")[-1].split("-")[0]

# Identify tRNA genes
tRNA_genes = {}

with open("Triticum_kronos.gff3", "r") as gff_file:
    for line in gff_file:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) > 2 and fields[2] == "tRNA":
            tRNA_genes[get_gene_id(line)] = ""

# Filter GFF3 file and count removed genes
with open("Triticum_kronos.filtered.gff3", "w") as output_file:
    with open("Triticum_kronos.gff3", "r") as gff_file:
        for line in gff_file:
            if line.startswith("#"):
                output_file.write(line)
                continue

            gene_id = get_gene_id(line)

            if gene_id in remove and gene_id not in tRNA_genes:
                output_file.write(line)

