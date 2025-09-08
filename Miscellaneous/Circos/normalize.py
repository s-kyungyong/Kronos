# normalize_snp_tracks_separately.py
def normalize_file(input_path, output_path, flip_sign=False):
    total = 0
    data = []

    with open(input_path) as f:
        for line in f:
            if line.strip() == "":
                continue
            chrom, start, end, count = line.strip().split()
            count = int(count)
            total += count
            data.append((chrom, int(start), int(end), count))

    with open(output_path, "w") as out:
        for chrom, start, end, count in data:
            norm = count / total * 10000
            if flip_sign:
                norm *= -1
            out.write(f"{chrom}\t{start}\t{end}\t{norm:.6f}\n")

# Example usage:
# Normalize exome data (positive values)
normalize_file("exome_density.bed", "exome_normalized.bed", flip_sign=False)

# Normalize promoter data (negative values)
normalize_file("promoter_density.bed", "promoter_normalized.bed", flip_sign=True)
