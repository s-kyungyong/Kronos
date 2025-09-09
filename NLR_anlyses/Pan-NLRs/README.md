

## NLR Prediction in Wheat Genomes

```
augustus v3.3.3
snap v2006-07-28
seqkit v2.8.2
cd-hit v4.8.1
```


---
### 1. Trainning Augustus  
To predict wheat NLRs, we first trained ab initio parameters using Kronos NLRs.  

**📥 Inputs**  
• `Kronos_NLRs.final.cds.fa`: Kronos NLR CDS  
• `Kronos_all.NLRs.final.gff3`: Kronos NLR gff  


**📥 Outputs**    
• `Wheat_NLR`: optimized AUGUSTUS parameters for Kronos NLRs  
  
⚙️ **Get high-confidence NLRs**    
```
#get high-confidence NLRs and reduce redundancy 
seqkit grep -f <(awk '$2 == "High" {print $1".1"}' NLR_confidence.list) Kronos_NLRs.final.cds.fa > Kronos_hc.cds.fa
cd-hit-est -c 0.9 -i Kronos_hc.cds.fa -o Kronos_hc.cds.cd-hit.c_0.9.fa -T 30
grep ">" Kronos_hc.cds.cd-hit.c_0.9.fa | sed 's/>//g' > geneIDs.list #this step generates 828 high-confidence nlrs

#create a separate gtf file for the selected high-confidence nlr genes
gffread -T --ids geneIDs.list Kronos_all.NLRs.final.gff3  > 828_hc_nlrs.gtf
awk '$3 == "CDS" {print}' 828_hc_nlrs.gtf  > 828_hc_nlrs.cds.gtf
```
⚙️ **Train AUGUSTUS**    
```
#create a genbank file with 1,000 bp flanking regions
perl gff2gbSmallDNA.pl 828_hc_nlrs.cds.gtf Kronos.collapsed.chromosomes.masked.v1.1.fa 1000 first.gb

#run quality check
#here, everything passes, but TrturKRN7A02G026020 got removed from the previous possibly as this gene is within an intron of another gene. 
etraining --species=generic first.gb > train.err

#separate the genes into 126 test set, and 700 trainning set
randomSplit.pl first.gb 126

#train augustus 
perl new_species.pl --species=Wheat_NLR
etraining --species=Wheat_NLR first.gb.train
optimize_augustus.pl --cpus=40 --kfold=40 --species=Wheat_NLR first.gb.train
```
---
We can compare how well the parameters do in predicting the NLRs in the test set. 

```
#this parameter is before the optimization step 
augustus --species=Wheat_NLR first.gb.test | tee firsttest.out
*******      Evaluation of gene prediction     *******
---------------------------------------------\
                 | sensitivity | specificity |
---------------------------------------------|
nucleotide level |       0.983 |       0.964 |
---------------------------------------------/

----------------------------------------------------------------------------------------------------------\
           |  #pred |  #anno |      |    FP = false pos. |    FN = false neg. |             |             |
           | total/ | total/ |   TP |--------------------|--------------------| sensitivity | specificity |
           | unique | unique |      | part | ovlp | wrng | part | ovlp | wrng |             |             |
----------------------------------------------------------------------------------------------------------|
           |        |        |      |                 82 |                106 |             |             |
exon level |    328 |    352 |  246 | ------------------ | ------------------ |       0.699 |        0.75 |
           |    328 |    352 |      |   44 |    2 |   36 |   46 |    4 |   56 |             |             |
----------------------------------------------------------------------------------------------------------/

----------------------------------------------------------------------------\
transcript | #pred | #anno |   TP |   FP |   FN | sensitivity | specificity |
----------------------------------------------------------------------------|
gene level |   137 |   126 |   72 |   65 |   54 |       0.571 |       0.526 |
----------------------------------------------------------------------------/
```

```
#this parameter is after the optimization step 
augustus --species=Wheat_NLR first.gb.test | tee firsttest.out
*******      Evaluation of gene prediction     *******
---------------------------------------------\
                 | sensitivity | specificity |
---------------------------------------------|
nucleotide level |       0.985 |        0.96 |
---------------------------------------------/

----------------------------------------------------------------------------------------------------------\
           |  #pred |  #anno |      |    FP = false pos. |    FN = false neg. |             |             |
           | total/ | total/ |   TP |--------------------|--------------------| sensitivity | specificity |
           | unique | unique |      | part | ovlp | wrng | part | ovlp | wrng |             |             |
----------------------------------------------------------------------------------------------------------|
           |        |        |      |                 89 |                 97 |             |             |
exon level |    344 |    352 |  255 | ------------------ | ------------------ |       0.724 |       0.741 |
           |    344 |    352 |      |   42 |    3 |   44 |   45 |    4 |   48 |             |             |
----------------------------------------------------------------------------------------------------------/

----------------------------------------------------------------------------\
transcript | #pred | #anno |   TP |   FP |   FN | sensitivity | specificity |
----------------------------------------------------------------------------|
gene level |   141 |   126 |   75 |   66 |   51 |       0.595 |       0.532 |
----------------------------------------------------------------------------/
```

Let's compare this to the other parameters we used in whole genome annotation. 


```
#auto trained parameters by braker in the version 1 annotation 
augustus --species=Kronos_collapsed first.gb.test | tee firsttest.out
*******      Evaluation of gene prediction     *******
---------------------------------------------\
                 | sensitivity | specificity |
---------------------------------------------|
nucleotide level |       0.861 |       0.919 |
---------------------------------------------/

----------------------------------------------------------------------------------------------------------\
           |  #pred |  #anno |      |    FP = false pos. |    FN = false neg. |             |             |
           | total/ | total/ |   TP |--------------------|--------------------| sensitivity | specificity |
           | unique | unique |      | part | ovlp | wrng | part | ovlp | wrng |             |             |
----------------------------------------------------------------------------------------------------------|
           |        |        |      |                369 |                225 |             |             |
exon level |    496 |    352 |  127 | ------------------ | ------------------ |       0.361 |       0.256 |
           |    496 |    352 |      |  174 |   75 |  120 |  161 |   21 |   43 |             |             |
----------------------------------------------------------------------------------------------------------/

----------------------------------------------------------------------------\
transcript | #pred | #anno |   TP |   FP |   FN | sensitivity | specificity |
----------------------------------------------------------------------------|
gene level |   174 |   126 |    0 |  174 |  126 |           0 |           0 |
----------------------------------------------------------------------------/
```

```
#manually trained parameters for ginger
augustus --species=Kronos_manual first.gb.test | tee firsttest.out

---------------------------------------------\
                 | sensitivity | specificity |
---------------------------------------------|
nucleotide level |       0.836 |       0.921 |
---------------------------------------------/

----------------------------------------------------------------------------------------------------------\
           |  #pred |  #anno |      |    FP = false pos. |    FN = false neg. |             |             |
           | total/ | total/ |   TP |--------------------|--------------------| sensitivity | specificity |
           | unique | unique |      | part | ovlp | wrng | part | ovlp | wrng |             |             |
----------------------------------------------------------------------------------------------------------|
           |        |        |      |                298 |                161 |             |             |
exon level |    489 |    352 |  191 | ------------------ | ------------------ |       0.543 |       0.391 |
           |    489 |    352 |      |  147 |   39 |  112 |  114 |    4 |   43 |             |             |
----------------------------------------------------------------------------------------------------------/

----------------------------------------------------------------------------\
transcript | #pred | #anno |   TP |   FP |   FN | sensitivity | specificity |
----------------------------------------------------------------------------|
gene level |   160 |   126 |   33 |  127 |   93 |       0.262 |       0.206 |
----------------------------------------------------------------------------/
```

```
#this is manually trained parameter for funannotate
augustus --AUGUSTUS_CONFIG_PATH /global/scratch/users/skyungyong/Software/anaconda3/envs/funannotate/config/ --species=Kronos_maker first.gb.test | tee firsttest.out
---------------------------------------------\
                 | sensitivity | specificity |
---------------------------------------------|
nucleotide level |       0.879 |       0.928 |
---------------------------------------------/

----------------------------------------------------------------------------------------------------------\
           |  #pred |  #anno |      |    FP = false pos. |    FN = false neg. |             |             |
           | total/ | total/ |   TP |--------------------|--------------------| sensitivity | specificity |
           | unique | unique |      | part | ovlp | wrng | part | ovlp | wrng |             |             |
----------------------------------------------------------------------------------------------------------|
           |        |        |      |                258 |                139 |             |             |
exon level |    471 |    352 |  213 | ------------------ | ------------------ |       0.605 |       0.452 |
           |    471 |    352 |      |  120 |   39 |   99 |   93 |    7 |   39 |             |             |
----------------------------------------------------------------------------------------------------------/

----------------------------------------------------------------------------\
transcript | #pred | #anno |   TP |   FP |   FN | sensitivity | specificity |
----------------------------------------------------------------------------|
gene level |   160 |   126 |   42 |  118 |   84 |       0.333 |       0.263 |
----------------------------------------------------------------------------/
```
----
### 2. Trainning SNAP 
**📥 Inputs**  
• `Kronos_hc.cds.cd-hit.c_0.9.fa`: filtered high-confidence NLRs from the previous step  
• `828_hc_nlrs.gtf`: filtered high-confidence NLRs in gff from the previous step  

**📥 Outputs**    
• `Wheat_NLR.hmm`: SNAP parameters for Kronos NLRs  

⚙️ **Train SNAP**  
```
#use randomly selected 700 genes for training 
grep ">"  Kronos_hc.cds.cd-hit.c_0.9.fa  | cut -d ">" -f 2 | shuf | head -n 700 > selected.hc.list
gffread --ids selected.hc.list  828_hc_nlrs.gtf > selected.hc.gff
awk '$3 != "exon" {print}' selected.hc.gff  > selected.hc.cds.gff
agat_convert_sp_gff2zff.pl --gff selected.hc.cds.gff  --fasta Kronos.collapsed.chromosomes.masked.v1.1.folded.fa -o genome.ann
#this step is to re-order the genome sequences to match the trainning sequences
grep ">" genome.ann  | cut -d ">" -f 2 | while read line; do awk -v seq=$line -v RS=">" '$1 == seq {print RS $0; exit}' Kronos.collapsed.chromosomes.masked.v1.1.fa; done > genome.dna

fathom -validate genome.ann genome.dna
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl Kronos_NLRs . > Wheat_NLR.hmm
```
---
### 3. NLR Prediction in Wheat Genomes
**📥 Inputs**  
• `genome.fa`: genome data available in **genome.list**. 

⚙️ **Donwload genomes** 
```
#download genomes
while IFS=$'\t' read -r col1 prefix url; do
    mkdir -p "$prefix" && cd "$prefix" || continue
    wget "$url"
    cd ..
done < genome.list
```

⚙️ **Isolate putative NLR loci** 
```
for prefix in */; do
  cd "${prefix}" || exit
  
  gz=$(ls *.gz) #some genomes are compressed differently. Change accordingly
  genome="${gz%.gz}"
  gunzip "$gz"

  #orf prediction and hmmsearch
  orfipy --procs 20 --bed orfs.bed --pep orfs.aa.fa --max 45000 --ignore-case --partial-3 --partial-5 "${genome}"
  hmmsearch --cpu 20 --domE 1e-4 -E 1e-4 --domtblout orfs.against.NBARC.out \
              PF00931.hmm "${genome}_out/orfs.aa.fa"

  #nlr-annotator
  java -jar ./NLR-Annotator/NLR-Annotator-v2.1b.jar \
      -t 20 -x ./NLR-Annotator/src/mot.txt \
      -y ./NLR-Annotator/src/store.txt \
      -i "${genome}" -o NLRannotator.whole-genome.out \
      -g NLRannotator.whole-genome.gff3
  done

  #crop genomes 
  python crop_genome.py --hmm orfs.aa.fa.against.NBARC.out --nlrannot NLRannotator.whole-genome.gff3 --genome 
  cd ..
done
```


Then, gene models were predicted with maker. The evidence and ab initio predictors used is different from the Kronos NLR prediction.
```
est evidence: transcript assemblies from stringtie that combined short and long-read transcriptome alignments. This was produced in the v2.0 annotation.
protein evidence: reliably curated Kronos NLRs, cloned functional NLRs, and NLRs from RefPlantNLR
ab initio prediction: Wheat_NLR, trained above.
```

Let's first create the evidence datasets. We can capture NLR sequences from long-read transcriptome data in **long-read_transripts.list**.
```
while read -r accession; do
    fa=${accession}.fasta
    
    #detect all orfs 
    orfipy -procs 56 --ignore-case --pep prot.fa ${fa}
    cd orfipy_${fa}_out
    
    #serch for nbarc domains and extract the hits
    hmmsearch --domtblout prot.fa.against.Kronos_NBARC.hmm.out --cpu 56 -E 1E-4 --domE 1e-4 orfipy_${fa}_out Kronos_NBARC.hmm prot.fa
    seqkit grep -f <(awk '!/^#/ {print $1}' prot.fa.against.Kronos_NBARC.hmm.out | sort -u) prot.fa > hits.fasta
    
    #remove the exact matches 
    cd-hit -c 1 -T 40 -i hits.fasta -o hits.reduced.fasta
    cd ..
done < long-read_transcripts.list
```

Combine all hits and remove exact matches again
```
#concatnate all NBARC-containging sequences into a single file
cat *.fasta_out/hits.reduced.fasta > all_hits/hits.reduced.combined.fa

#remove redundancy and reextract the hits 
cd-hit -c 1.0 -M 6000000 -T 40 -i hits.reduced.combined.fa -o hits.reduced.combined.reduced.fa
seqkit grep -f <( grep ">" hits.reduced.combined.reduced.fa | cut -d "_" -f 1 | sort -u) ../*.fasta > hits.reduced.combined.reduced.est.fa
grep ">" hits.reduced.combined.reduced.fa | cut -d "_" -f 1 | sort -u | sed 's/>//g' > hit_ids.txt
awk -F"[_\\.]" '{gsub(/^>/, "", $1); print $1"."$2 > "hit_ids_"$1".txt"}' hit_ids.txt

for fa in ../*.fasta; do
    prefix=$(basename "$fa" .fasta)
    if [ -f "hit_ids_${prefix}.txt" ]; then
        seqkit grep -f "hit_ids_${prefix}.txt" "$fa" >> hits.reduced.combined.reduced.est.fa
    fi
done
```
 
For each locus, a separate directory was made to run maker in parallel. We ran > 500 jobs at the same time to speed up this step. This still took over 3 days. 
```
#for each genome
dir=$1

cd "$dir" || exit
mkdir -p split_genome

# Split the genome sequences into individual FASTA files
seqkit split -i -O split_genome NLR_loci.fa

cd split_genome || exit
for fa in *.fa; do
  prefix="${fa%.fa}"  # Extract prefix by removing ".fa"
  mkdir -p "$prefix"
  mv "$fa" "$prefix"/
  cp /global/scratch/users/skyungyong/Kronos/NLR_annotations/Pan-NLRome/Evidence/maker* "$prefix"/ #copy control files
#  maker -RM_off -genome ${dir}.fa 
done
```

Once all jobs were finished, gff3 files were collected and genes were extracted. 
```
#create gff3 outputs
ls -d * | parallel -j 40 'gff3_merge -d "{}/{}.maker.output/{}_master_datastore_index.log" -o "{}/{}.gff3" && echo "Processed: {}"'

#collect raw outputs into a separate folder
find split_genome -type f -name "*.gff3" -print0 | xargs -0 mv -t raw_gff/

#combine gene annotations into a single gff3 file
find raw_gff -type f -name "*.gff3" -exec grep -E 'gene|exon|mRNA|CDS' {} + | awk '$3 == "gene" || $3 == "exon" || $3 == "mRNA" || $3 == "CDS"' | cut -d ":" -f 2- > NLR_loci.maker.gff3

awk -F ':' '{print $2}' NLR_loci.maker.gff3 > cleaned_NLR_loci.maker.gff3

#get protein sequences
for dir in */; do
    if [[ -d "$dir" ]]; then
        pushd "$dir" || continue
        if [[ -f "NLR_loci.maker.gff3" && -f "NLR_loci.fa" ]]; then
            cut -d ":" -f 2- NLR_loci.maker.gff3 > cleaned_NLR_loci.maker.gff3
            gffread -x NLR_loci.maker.cds.fa -y NLR_loci.maker.pep.fa -g NLR_loci.fa cleaned_NLR_loci.maker.gff3
            hmmsearch --domtblout NLR_loci.maker.pep.against.Kronos_NBARC.hmm.out -E 1e-4 --domE 1e-4 --cpu 56 /global/scratch/users/skyungyong/Kronos/NLR_annotations/HMM/Kronos_NBARC.hmm NLR_loci.maker.pep.fa
            diamond blastp -q NLR_loci.maker.pep.fa -d /global/scratch/users/skyungyong/Kronos/NLR_annotations/HMM/blastdb/Kronos.NLRs.reliable.dmnd --evalue 1e-4 --masking 0 --max-target-seqs 1 --max-hsps 1 --threads 56 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen -o NLR_loci.maker.pep.against.Kronos_NBARC.dmnd.out
        else
            echo "Skipping $dir: Required files not found."
        fi
        popd
    fi
done

```
