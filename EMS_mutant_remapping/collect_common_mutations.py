import gzip
import glob
from collections import defaultdict
import sys
import os

ems = sys.argv[1]
if ems not in ("True","False"):
    sys.exit("ems must be True/False")
is_ems = ems == "True"

min_rate = 0.05 # 0.05 for exome capture and 0.15 for promoter capture 
common = {}

for line in open('exome-promoter.common.list', 'r'):
    common[ line.strip() ] = ''

vcf_files = [ x for x in os.listdir() if x.endswith('filtered.vcf.gz') and x.split('.')[0] in common ]
site_mutant_count = defaultdict(int)

#count how many mutants carry a mutation at each site
valid_gts = {"0/1","1/0","0|1","1|0","1/1","1|1"}
for vcf_path in vcf_files:
    seen_sites = set()
    prefix = vcf_path.split('.')[0]
  
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            chrom, pos, fmt, sample_data = parts[0], parts[1], parts[8], parts[9]
            wt, mut = parts[3], parts[4]

            is_ems_mut = (wt == "G" and mut == "A") or (wt == "C" and mut == "T")
            if is_ems != is_ems_mut:   # XOR logic
                continue               # skip site
            
          site = f"{chrom}_{pos}"

            fmt_fields = fmt.split(":")
            sample_fields = sample_data.split(":")
            if "GT" not in fmt_fields:
                continue
            gt = sample_fields[fmt_fields.index("GT")]

            #only consider biallelic sites
            #take zygosity for the mutation
            #hetrozygous: 0/1, 1/0 0|1 1|0
            #homozygous for alternative allel: 1/1 1|1
            if gt not in valid_gts:
                continue

            seen_sites.add(site)

    for site in seen_sites:
        site_mutant_count[site] += 1


min_required = len(vcf_files) * min_rate
selected_sites = {s for s, c in site_mutant_count.items() if c >= min_required}
print(f"Selected {len(selected_sites)} variant sites present in ≥{min_required} mutants.")

# Pass 2 — output relevant lines
with open(f"common_EMS_{ems}_per{min_rate}_{len(selected_sites)}_sites_{min_required}_mutants.tsv", "w") as out:
    out.write("sample\tsite\tgenotype\n")
    for vcf_path in vcf_files:
        sample = vcf_path.split("/")[-1].split(".")[0]
        with gzip.open(vcf_path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                chrom, pos, fmt, sample_data = parts[0], parts[1], parts[8], parts[9]
                site = f"{chrom}_{pos}"
                if site not in selected_sites:
                    continue

                fmt_fields = fmt.split(":")
                sample_fields = sample_data.split(":")
                if "GT" not in fmt_fields:
                    continue
                gt = sample_fields[fmt_fields.index("GT")]

                if gt in ("0/1", "1/0", "0|1", "1|0", "1/1", "1|1"):
                    out.write(f"{sample}\t{site}\t{gt}\n")
