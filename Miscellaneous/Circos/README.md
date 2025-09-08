# Circos

## Software version

```
circos v0.69-8
```

------


### 1. Chromosomal synteny 

```
#for each chromosome pairs, $chr1 and $chr2
minimap2 -f 0.05 -t 40 -x asm10 Kronos.collapsed.chromosomes.masked.v1.1.${chr1}.fa Kronos.collapsed.chromosomes.masked.v1.1.${chr2}.fa > ${chr1}_vs_${chr2}.minimap.paf
```

```
#change paf files to colored links 
#!/bin/bash

# Define color map based on query chromosomes
declare -A colors
colors["1A"]="red"
colors["1B"]="darkred"
colors["2A"]="orange"
colors["2B"]="darkorange"
colors["3A"]="yellow"
colors["3B"]="olive"
colors["4A"]="green"
colors["4B"]="darkgreen"
colors["5A"]="blue"
colors["5B"]="steelblue"
colors["6A"]="indigo"
colors["6B"]="slateblue"
colors["7A"]="violet"
colors["7B"]="purple"

# Output dir
mkdir -p circos_links

# Loop through all paf files
for paf in *.minimap.paf; do
    # Extract filename components
    fname=$(basename "$paf")
    query_chr=$(echo "$fname" | cut -d'_' -f1) # e.g. 1A
    color=${colors[$query_chr]}

    # Output file
    out="circos_links/${fname%.minimap.paf}.links.txt"

    # Filter and format
    awk -v col="$color" '$13 == "tp:A:P" && $11 > 150000 {
        print $1, $3, $4, $6, $8, $9, "color="col
    }' OFS='\t' "$paf" > "$out"
done
```

### 2. Repeat contents
```
#hard mask genomes using annotated repeats in Kronos.collapsed.chromosomes.masked.v1.1.fa.mod.EDTA.TEanno.gff3
bedtools makewindows -g Kronos.collapsed.chromosomes.masked.v1.1.fa.fai -w 1000000 > windows.bed
bedtools nuc -fi Kronos.collapsed.chromosomes.masked.v1.1.hard_masked_with_X.fa -pattern X -bed windows.bed | awk 'BEGIN{OFS="\t"} NR>1 { print $1, $2, $3, $13/$12 }'> repeat_density.txt
```

### 3. Gene contents
```
#use annotation version 2.1 (gff to bed)
bedtools intersect -a windows.bed -b genes.v21.bed -c | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > gene_density.bed
```
