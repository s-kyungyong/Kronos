# Repeat Annotation
## Data Availability
Repeat annotations are available in Zenodo.   
• `Repeat annotation [✨Final✨]`: https://zenodo.org/records/15399687  

## Software version
```
HiTE v3.0
EDTA v2.2.2
RepeatMasker v4.1.5
```

---

## 1. Repeat Identification with HiTE
Repeative elements annotated in the Kronos reference genome v1.0 and v1.1 were predicted by [HiTE](https://github.com/CSU-KangHu/HiTE). 

📥 Inputs  

• `Kronos.collapsed.chromosomes.fa`: Haplotype-collapsed chromosomes (v1.0)  

📥 Outputs  

• `Kronos.collapsed.chromosomes.masked.fa`: Haplotype-collapsed masked chromosomes (v1.0)  

⚙️ **Run HiTE**  
```
singularity run HiTE.sif python main.py --genome Kronos.collapsed.chromosomes.fa \
     --thread 56 --outdir HiTE --recover 1 --annotate 1 \
     --plant 1 --classified 1 --domain 1 --recover 1 --debug 1
```

⚙️ **Soft-mask with RepeatMasker**  
```
RepeatMasker -xsmall -e ncbi -pa 56 -q -no_is -norna -nolow -div 40 -gff -lib confident_TE.cons.fa.classified -cutoff 225 Kronos.collapsed.chromosomes.fa
```

---


## 2. Repeat Identification with EDTA
The next version of repeat annotations was produced using EDTA v2.2.2. 

📥 Inputs
• `Kronos.collapsed.chromosomes.v1.1.fa`: Haplotype-collapsed chromosomes (v1.1)
• `trep-db_complete_Rel-19.fasta`: Curated repeat database
• `confident_TE.cons.fa.classified`: HiTE-derived repeat annotation
• `Kronos.v2.0.pep.fa Kronos.v2.0.cds.fa`: Kronos v2.0 annotaion


⚙️ **Generate Databases**  
• Annotated TREP Databases 
```
wget https://trep-db.uzh.ch/downloads/trep-db_complete_Rel-19.fasta.gz
gunzip trep-db_complete_Rel-19.fasta.gz
python filter_trep.py #this creates trep-db_complete_Rel-19.triticum.filtered.fa

#add telomeric repeats to the file
>Generic_telomere#telomere/telomere 
TTTAGGGTTTAGGGTTTAGGGTTTAGGG
```

• Pre-curated HiTE annotations 
```
python filter_hite.py #get confident_TE.cons.fa.classified from HiTE
```

• Pre-curated HiTE annotations 
```
diamond makedb -d Kronos.v2.0.pep.dmnd --in Kronos.v2.0.pep.fa
diamond blastx -d Kronos.v2.0.pep.dmnd --out trep-db_complete_Rel-19.triticum.filtered.against.Kronosv2.dmnd.out \
     --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
     -q trep-db_complete_Rel-19.triticum.filtered.fa
python filter_cds.py #get clean cds 
```

⚙️ **Run EDTA**  

• Annotate non-TIR elements 
Note: EDTA has a really long runtime as discussed [here](https://github.com/oushujun/EDTA/issues/61). 
```
#Let's divde and conquer. For LTR, LINE, SINE and helitron:
EDTA_raw.pl --genome Kronos.collapsed.chromosomes.masked.v1.1.fa --species others --type ${target} -t 40 --overwrite 0 --rmlib confident_TE.cons.fa.classified.filtered.fa
```

• Annotate TIR elements 
Note: TIR prediction is extremely slow. The runtime for the bread wheat was 4 weeks!
```
#run for each chromosome to speed up
EDTA_raw.pl --genome Kronos.v1.1.${target}.fa --species others --type tir -t 20 --overwrite 0 --rmlib confident_TE.cons.fa.classified.filtered.fa

#merge all EDTA outputs
cat genome_split_for_TIR/*.raw/*mod.TIR.intact.raw.fa > Kronos.collapsed.chromosomes.masked.v1.1.fa.mod.TIR.intact.raw.fa
cat genome_split_for_TIR/*.raw/*mod.TIR.intact.raw.bed > Kronos.collapsed.chromosomes.masked.v1.1.fa.mod.TIR.intact.raw.bed
cat .genome_split_for_TIR/*.raw/*mod.TIR.intact.raw.gff3 > Kronos.collapsed.chromosomes.masked.v1.1.fa.mod.TIR.intact.raw.gff3
cat genome_split_for_TIR/*.raw/*mod.TIR.intact.raw.fa.anno.list > Kronos.collapsed.chromosomes.masked.v1.1.fa.mod.TIR.intact.raw.fa.anno.list
```

• Finish EDTA  
```
EDTA.pl --genome Kronos.collapsed.chromosomes.v1.1.fa --species others --step filter \
        --cds Kronos.v2.0.cds.filtered.fa --curatedlib trep-db_complete_Rel-19.triticum.filtered.fa \
        --rmlib confident_TE.cons.fa.classified.filtered.fa --sensitive 1 --anno 1 \
        --evaluate 1 --overwrite 0 --threads 56
```

