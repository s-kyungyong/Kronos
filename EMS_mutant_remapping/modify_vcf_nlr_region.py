import os
import sys
from Bio import SeqIO
from intervaltree import IntervalTree

# --- Inputs ---
chromosome = sys.argv[1]  # e.g. "1A"
vcf_dir = "Exome-capture_GATK_vcf_annot_v2.1/" #a folder with GATK derived vcfs
nlrs = "nlr_regions.fa" # target regions annotated as "1A:12345-45678"

# --- Load NLR regions and build interval tree ---
trees = {}
mutations = {}  # { region_id: { sample_id: [vcf_lines] } }

print(f"Building interval tree for {chromosome}...")
trees[chromosome] = IntervalTree()

for record in SeqIO.parse(nlrs, "fasta"):
    chrom_region = record.id  # e.g. "1A:12345-45678"
    chrom, coords = chrom_region.split(":")
    if chrom != chromosome:
        continue

    start, end = map(int, coords.split("-"))
    start += 1
    region_id = chrom + "_" + str(start) + "-" + str(end)

    trees[chrom][start:end] = region_id
    mutations[region_id] = {}  # Init per region

    # Save reference FASTA
    os.makedirs(region_id, exist_ok=True)
    ref_path = os.path.join(region_id, f"{region_id}.fa")

    with open(ref_path, 'w') as o:
        o.write(f'>{region_id}\n{record.seq}')

# --- Parse VCF files and assign lines to regions ---
print("Scanning VCF files...")
headers = []

for vcf_name in os.listdir(vcf_dir):
    if not vcf_name.endswith(".vcf") or vcf_name.startswith("Kronos0"):
        continue

    sample_id = vcf_name.split('.')[0]
    vcf_path = os.path.join(vcf_dir, vcf_name)

    with open(vcf_path) as vcf_file:
        for line in vcf_file:
            if line.startswith("#"):
                if not headers:
                    headers.append(line)
                continue

            fields = line.strip().split("\t")
            vcf_chrom = fields[0]
            pos = int(fields[1])
            quality = float( fields[5] )

            if vcf_chrom != chromosome:
                continue

            if not quality > 30:
                continue

            dp = int( fields[7].split('DP=')[1].split(';')[0] )

            if not dp > 5:
                continue

            hits = trees[chromosome][pos]
            for hit in hits:
                region_key = hit.data
                if sample_id not in mutations[region_key]:
                    mutations[region_key][sample_id] = []
                mutations[region_key][sample_id].append(line)

# --- Write region-specific VCFs ---
print("Writing region-specific VCFs...")

for region in mutations:
    start = int(region.split('_')[1].split('-')[0])
    for mutant in mutations[region]:
        vcf_lines = mutations[region][mutant]
        if vcf_lines:
            with open(f"{region}/{mutant}.vcf", 'w') as o:
                o.write("##fileformat=VCFv4.2\n")
                o.write("\t".join("#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Kronos".split()))
                for line in vcf_lines:
                    fields = line.split("\t")
                    fields[0] = region
                    fields[1] = str( int(fields[1]) - start + 1)

                    o.write("\n" + "\t".join(fields) )

    vcf_files = [ v for v in os.listdir(region) if v.endswith('vcf') ]
    if len(vcf_files) == 0:
        os.system(f'rm -r {region}')

print("Done.")
