# Transcriptome analyses





## Software version

```
fastp v0.24.0
salmon v1.10.3
```

---
### 1. Quality Control

**📥 Inputs**  
• `*.fastq`: raw RNA-seq data from NCBI

**📥 Outputs**  
• `*.filtered.fastq`: filtered, trimmed reads

⚙️ **Trim with fastp**  
```
#for paired-end 
for fq1 in *_1.fastq; do
  fq2=$(echo $fq1  | sed 's/_1.fastq/_2.fastq/g')
  out1=$(echo $fq1 | sed 's/fastq/filtered.fastq/g')
  out2=$(echo $fq2 | sed 's/fastq/filtered.fastq/g')
  fastp --in1 $fq1 --out1 $out1 --in2 $fq2 --out2 $out2 -q 20 --length_required 50 --detect_adapter_for_pe -w 16
done

#for single end
for fq1 in *.fastq; do
  out1=$(echo $fq1 | sed 's/fastq/filtered.fastq/g')
  fastp --in1 $fq1 --out1 $out1 -q 20 --length_required 50 -w 16
done
```
---
### 2. Quantification
**📥 Inputs**  
• `*.filtered.fastq`: filtered, trimmed reads  
• `Kronos.collapsed.chromosomes.masked.v1.1.fa`: reference genome  
• `Kronos.v2.1.gff3`: v2.1 annotations  

**📥 Outputs**  
• `quant.genes.sf`: Gene-level quantifications    

⚙️ **Indexing**  
```
#prepare database
grep ">" Kronos.collapsed.chromosomes.masked.v1.1.fa | cut -d " " -f 1 | cut -d ">" -f 2 > decoys.txt
gffread -w Kronos.v2.1.transcripts.fa -g Kronos.collapsed.chromosomes.masked.v1.1.fa Kronos.v2.1.gff3
cat Kronos.v2.1.transcripts.fa Kronos.collapsed.chromosomes.masked.v1.1.fa > Kronos.gentrome.fa

#index gentrome
salmon index -t Kronos.gentrome.fa -d decoys.txt -p 30 -i salmon_index
```
⚙️ **Quantification**  
```
#for paired-end data:
indir=$1
for fq1 in ${indir}/*_1.filtered.fastq; do
  prefix=$(basename "$fq1" | cut -d "_" -f 1) 
  fq2=$(echo $fq | sed 's/_1.filtered/_2.filtered/g')

  salmon quant -l A -1 "$fq1" -2 "$fq2" -p 40 \
    -g Kronos.v2.1.gtf \
    -i salmon_index/ \
    -o "${prefix}" --validateMappings
done
```
```
#for single_end
indir=$1

for fq1 in ${indir}/*.filtered.fastq; do
  prefix=$(basename "$fq1" | cut -d "." -f 1) 

  salmon quant -l A -r "$fq1" -p 40 \
    -g Kronos.v2.1.gtf \
    -i salmon_index/ \
    -o "${prefix}" --validateMappings
done
```


