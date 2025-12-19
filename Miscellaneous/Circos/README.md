# Circos

## Software version

```
circos v0.69-8
minimap v2.28-r1209
bedtools v2.31.1
```

------


### 1. Chromosomal synteny 

```
#DNA-based
#for each chromosome pairs, $chr1 and $chr2
minimap2 -f 0.05 -t 40 -x asm10 Kronos.collapsed.chromosomes.masked.v1.1.${chr1}.fa Kronos.collapsed.chromosomes.masked.v1.1.${chr2}.fa > ${chr1}_vs_${chr2}.minimap.paf
```

```
#protein-based
diamond makedb --in Kronos.v2.1.pep.longest.fa -d Kronos.v2.1.pep.longest.dmnd

diamond blastp -q Kronos.v2.1.pep.longest.fa -d Kronos.v2.1.pep.longest.dmnd \
               -o Kronos.blast -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
               --masking 0 --max-target-seqs 5 --evalue 1e-10 --threads 56
MCScanX -b 0 -a Kronos

##reformat collinearity -> anchor: 3 fields, with block sepearted by #
head Kronos.anchor
#
TrturKRN1A02G023390.2   TrturKRN4A02G051510.1   0
TrturKRN1A02G023410.1   TrturKRN4A02G051580.1   9e-52
TrturKRN1A02G023550.1   TrturKRN4A02G051780.1   9e-117
TrturKRN1A02G023900.1   TrturKRN4A02G051980.1   7e-126
TrturKRN1A02G023920.1   TrturKRN4A02G052210.1   2e-27
TrturKRN1A02G023940.1   TrturKRN4A02G052270.1   1e-169
TrturKRN1A02G024020.1   TrturKRN4A02G052410.1   0
TrturKRN1A02G024030.1   TrturKRN4A02G052440.1   5e-266
#

#minimum 10 gene blocks
python -m jcvi.compara.synteny screen Kronos.anchor --minspan=10 Kronos.jcvi.anchors --qbed=Kronos.bed --sbed=Kronos.bed
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
    out="${fname%.minimap.paf}.links.txt"

    # Filter and format
    awk -v col="$color" '$13 == "tp:A:P" && $11 > 150000 {
        print $1, $3, $4, $6, $8, $9, "color="col
    }' OFS='\t' "$paf" > "$out"
done

cat *.links.txt > 00.synteny_links.txt
```

### 2. Repeat density
```
#hard mask genomes using annotated repeats in Kronos.collapsed.chromosomes.masked.v1.1.fa.mod.EDTA.TEanno.gff3
bedtools makewindows -g Kronos.collapsed.chromosomes.masked.v1.1.fa.fai -w 1000000 > windows.bed
bedtools nuc -fi Kronos.collapsed.chromosomes.masked.v1.1.hard_masked_with_X.fa -pattern X -bed windows.bed | awk 'BEGIN{OFS="\t"} NR>1 { print $1, $2, $3, $13/$12 }'> 01.repeat_density.txt
```

### 3. Gene density
```
#use annotation version 2.1 (gff to bed)
bedtools intersect -a windows.bed -b genes.v21.bed -c | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > 02.gene_density.txt
```

### 4. NLR density
```
#use high and medium confidence NLRs (gff to bed)
bedtools intersect -a windows.bed -b NLRs.reliable.bed -c | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > 03.NLR_density.txt
```

### 5. EMS mutation density
```
#use MAPS-drivin EMS mutation density
awk '{print $1 "\t" $2 "\t" $2+1 "\t" 1}' Kronos_v1.1.Exom-capture.corrected.deduped.10kb_bins.RH.byContig.MI.No_RH.maps.substitutions.vcf > exom.maps.ems.bed
bedtools intersect -a windows.bed -b exom.maps.ems.bed -c | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > exome_density.bed
awk '{print $1 "\t" $2 "\t" $2+1 "\t" 1}' Kronos_v1.1.Promoter-capture.corrected.deduped.10kb_bins.RH.byContig.MI.No_RH.maps.substitutions.vcf > promoter.maps.ems.bed
bedtools intersect -a windows.bed -b promoter.maps.ems.bed -c | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > promoter_density.bed
python normalize.py
```

### 6. PHAS and MIR loci
```
#use sRNA datasets 
awk '$3 == "21PHAS" { print $1, $4, $5, "green" } $3 == "24PHAS" { print $1, $4, $5, "red" }'phas.gff3 > 05.PHAS_loci.txt
awk '$3 == "MIRNA_hairpin" {print $1 "\t" $4 "\t" $5 }' miRNAs.gff3  | sed 's/Chr//g' > MIR_loci.bed
bedtools intersect -a windows.bed -b MIR_loci.bed -c | awk '$4 > 0 {print $1 "\t" $2 "\t" $3 "\t" $4}' > 05.MIR_loci.txt
```


### 7. Circos
```
#with *.conf files in the current folder
#check out the output circos.png
#we further edited this backbone in Illustrator to create Fig. 1. 
circos
```

