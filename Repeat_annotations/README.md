# Repeat Annotation
## Data Availability


## Method


Repetitive elements were initially annotated with HiTE v3.0.0 (Hu et al., 2024) and used to soft-mask the reference genomes v1.0 and v1.1. After generating the reference annotation v2.0, we re-annotated repetitive elements with EDTA v2.2.2 (Ou et al., 2019). To enhance repeat prediction and classification, complete and consensus repeats for Triticum were retrieved from the TREP database and included as curated libraries (Schlagenhauf and Wicker et al., 2016). Additionally, classified repeats from HiTE were integrated as RepeatModeler libraries. To prevent over-masking, the coding sequences of the v2.0 annotations were also provided. 


---

## 1. Repeat Identification with HiTE
Repeative elements annotated in the Kronos reference genome v1.0 and v1.1 are predicted by [HiTE](https://github.com/CSU-KangHu/HiTE) v3.0.0. 
```
singularity run HiTE.sif python main.py --genome Kronos.collapsed.chromosomes.fa \
     --thread 56 --outdir HiTE --recover 1 --annotate 1 \
     --plant 1 --classified 1 --domain 1 --recover 1 --debug 1

#mask v1.0 genome
RepeatMasker -xsmall -e ncbi -pa 56 -q -no_is -norna -nolow -div 40 -gff -lib confident_TE.cons.fa.classified -cutoff 225 Kronos.collapsed.chromosomes.fa
```

## 2. Repeat Annotations with EDTA
The next version of repeat annotations is produced by EDTA v2.2.2. Some input files need to be prepared. 

### Annotated TREP Databases 
Download TREP database and filter out some sequences. We will only use complete or consensus annotations from *Triticum* and exlcudes any unknown classes. The output will be used as --curatedlib.
```
wget https://trep-db.uzh.ch/downloads/trep-db_complete_Rel-19.fasta.gz
gunzip trep-db_complete_Rel-19.fasta.gz
python filter_trep.py #this creates trep-db_complete_Rel-19.triticum.filtered.fa

#add telomeric repeats to the file
>Generic_telomere#telomere/telomere 
TTTAGGGTTTAGGGTTTAGGGTTTAGGG
```

Let's create an inputfile that would be used for --rmlib. We will use the repeat annotation file from HiTE. Similarly, let's exclude any unknown classes. 
```
#get confident_TE.cons.fa.classified from HiTE
python filter_hite.py
```

Lastly, CDS of annotated genes will be used for --cds. We will remove matches to TEs from v2.0 annotation sets. The criteria will be set very lenient with E-value < 1e-04 and 40% coverage for the Kronos genes. This does not necessarily mean that these genes are TEs. We are simply giving EDTA options. 
```
#run diamond v2.1.9.163
diamond makedb -d Kronos.v2.0.pep.dmnd --in Kronos.v2.0.pep.fa
diamond blastx -d Kronos.v2.0.pep.dmnd --out trep-db_complete_Rel-19.triticum.filtered.against.Kronosv2.dmnd.out --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen -q trep-db_complete_Rel-19.triticum.filtered.fa
python filter_cds.py
```

Now, let's run EDTA. EDTA has a really long runtime as discussed [here](https://github.com/oushujun/EDTA/issues/61). Let's divde and conquer. For LTR, LINE, SINE and helitron:
```
#target is eitheir ltr, line, sine or helitron 
EDTA_raw.pl --genome Kronos.collapsed.chromosomes.masked.v1.1.fa --species others --type ${target} -t 40 --overwrite 0 --rmlib confident_TE.cons.fa.classified.filtered.fa
```

TIR prediction is extremely slow. The runtime for the bread wheat was 4 weeks, based on the post above. After running this step for 8 days, We also decided that we need to split genome and run EDTA individually to speed this process up. 
```
#target is each chromosome: 1A, 1B ... 7A, 7B, and Un
EDTA_raw.pl --genome Kronos.v1.1.${target}.fa --species others --type tir -t 20 --overwrite 0 --rmlib confident_TE.cons.fa.classified.filtered.fa

#merge all EDTA outputs
cat ../genome_split_for_TIR/*.raw/*mod.TIR.intact.raw.fa > Kronos.collapsed.chromosomes.masked.v1.1.fa.mod.TIR.intact.raw.fa
cat ../genome_split_for_TIR/*.raw/*mod.TIR.intact.raw.bed > Kronos.collapsed.chromosomes.masked.v1.1.fa.mod.TIR.intact.raw.bed
cat ../genome_split_for_TIR/*.raw/*mod.TIR.intact.raw.gff3 > Kronos.collapsed.chromosomes.masked.v1.1.fa.mod.TIR.intact.raw.gff3
cat ../genome_split_for_TIR/*.raw/*mod.TIR.intact.raw.fa.anno.list > Kronos.collapsed.chromosomes.masked.v1.1.fa.mod.TIR.intact.raw.fa.anno.list
```

Finish EDTA.
```
EDTA.pl --genome Kronos.collapsed.chromosomes.masked.v1.1.fa --species others --step filter \
        --cds Kronos.v2.0.cds.filtered.fa --curatedlib trep-db_complete_Rel-19.triticum.filtered.fa \
        --rmlib confident_TE.cons.fa.classified.filtered.fa --sensitive 1 --anno 1 \
        --evaluate 1 --overwrite 0 --threads 56
```

