# Repeat Annotation

### Repeat Identification with HiTE
Repeative elements annotated in the Kronos reference genome v1.0 and v1.1 are predicted by [HiTE](https://github.com/CSU-KangHu/HiTE) v3.0.0. 
```
singularity run HiTE.sif python main.py --genome Kronos.collapsed.chromosomes.fa \
     --thread 56 --outdir HiTE --recover 1 --annotate 1 \
     --plant 1 --classified 1 --domain 1 --recover 1 --debug 1
```

### Repeat Annotations with EDTA
The next version repeat annotation is produced by EDTA v2.2.2. Some input files need to be prepared. 

### Annotated TREP Databases 
Download TREP database and filter out some sequences. We will only use annotations from *Triticum* and exlcudes any unknown classes. The output will be used as --curatedlib.
```
wget https://trep-db.uzh.ch/downloads/trep-db_complete_Rel-19.fasta.gz
gunzip trep-db_complete_Rel-19.fasta.gz
python filter_trep.py #this creates trep-db_complete_Rel-19.triticum.filtered.fa

#add telomeric repeats
#to the file
>Generic_telomere#telomere/telomere 
TTTAGGGTTTAGGGTTTAGGGTTTAGGG
```

Let's create an inputfile that would be used for --rmlib. We will use the repeat annotation file from HiTE. 
```
#get confident_TE.cons.fa.classified from HiTE
python filter_hite.py
```

Lastly, CDS of annotated genes will be used for --cds. We will remove matches to TEs from v2.0 annotation sets based on some criteria and let EDTA decide what to do for those loci.
```
diamond makedb -d Kronos.v2.0.pep.dmnd --in Kronos.v2.0.pep.fa
diamond blastx -d Kronos.v2.0.pep.dmnd --out trep-db_complete_Rel-19.triticum.filtered.against.Kronosv2.dmnd.out --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen -q trep-db_complete_Rel-19.triticum.filtered.fa
python filter_cds.py
```

Now, let's run EDTA. 
```
singularity exec -B $(pwd) EDTA.sif EDTA.pl --genome Kronos.collapsed.chromosomes.masked.v1.1.fa --species others --step all \
                    --cds Kronos.v2.0.cds.filtered.fa --curatedlib trep-db_complete_Rel-19.triticum.filtered.fa \
                    --rmlib confident_TE.cons.fa.classified.filtered.fa --sensitive 1 --annot 1 \
                    --evaluate 1 --threads 56
```

