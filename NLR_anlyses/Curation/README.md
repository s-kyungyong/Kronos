# NLR Curation

TrturKRN2B02G099990 is a low-confidence NLR

## Data Availability
Curated NLRs were deposited in [Zenodo](https://zenodo.org/records/15539721).

## Software Versions
```
orfipy v0.0.4
hmmsearch v3.4
nlr-annotator v2.1b
maker v3.01.03
seqkit v
star v2.7.11b
deeptools v3.5.5
samtools v
interproscan v5.68-100.0
apllo v2.0.6
agat v0.8.0
```

---

### 1. Putative NLR Loci Detection
The Kronos genome is large, and for targeted manual curation, NLR loci need to be first extracted from the genome. NLRs were loosely defined as NB-ARC domain-containing genes or proteins. Genomic regions containing NB-ARC domains were identified and extracted with 15,000 flanking sequences from both ends. You may choose to increase the flanking size, as a small number of genes could not be fully contained in this region.

**📥 Inputs**  
• `Kronos.collapsed.chromosomes.masked.v1.1.fa`: Kronos genome  

**📥 Outputs**  
• `Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa`: NLR loci  


⚙️ **identify open reading frames**  
```
orfipy --procs 56 --bed orfs.bed --pep orfs.aa.fa --max 45000 --ignore-case --partial-3 --partial-5 Kronos.collapsed.chromosomes.masked.v1.1.fa
```
⚙️ **Search for NB-ARC domains**  
```
hmmsearch --cpu 56 --domE 1e-4 -E 1e-4 --domtblout orfs.against.NBARC.out PF00931.hmm orfs.aa.fa
```
⚙️ **Extract NLR loci**  
We used 15,000 flanking regions, but some NLRs had really long introns could could not be contained in extracted sequences. It may better to increase this threshold. 
```
python crop_genome.py --hmm orfs.aa.fa.against.NBARC.out --genome Kronos.collapsed.chromosomes.masked.v1.1.fa
```

⚙️ **NLR annotator**  
We later learned that some divergent NB-ARC domains cannot be properly detected by the HMM approach and additionally incoporated NLR-Annotator.
```
java -jar NLR-Annotator-v2.1b.jar -t 40 -x ./NLR-Annotator/src/mot.txt -y ./NLR-Annotator/src/store.txt -i Kronos.collapsed.chromosomes.masked.v1.1.fa -o NLRannotator.whole-genome.out -g NLRannotator.whole-genome.gff3

#crop genome to NB-ARC-containing loci
python crop_genome.py --nlrannot NLRannotator.whole-genome.gff3 --genome Kronos.collapsed.chromosomes.masked.v1.1.fa
```

Alternatively, the two outputs can be combined.
```
python crop_genome.py --hmm orfs.aa.fa.against.NBARC.out --nlrannot NLRannotator.whole-genome.gff3 --genome Kronos.collapsed.chromosomes.masked.v1.1.fa
```
---
### 2. Gene model prediction
Initial gene models were predicted with MAKER. 

**📥 Inputs**  
• `Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa`: NLR loci  
• `proteins.fa`: NLRs from 18 Poaceae species [(source)](https://zenodo.org/records/13627395) and 415 reference NLRs from [RefPlantNLR](https://zenodo.org/records/3936022).
• `est.fa`: Stringtie transcripts assememblies (v1.0 annotation (short reads) and v2.0 annotation (short/long-reads). See [this record](https://github.com/s-kyungyong/Kronos/tree/main/Genome_annotation)
• `Kronos_collapsed` & `Kronos.hmm`: The two ab initio parameters obtained in the v1 annotation: Augustus from BRAKER and SANP from GINGER. See [this folder](https://github.com/s-kyungyong/Kronos/tree/main/Genome_annotation/abinitio_parameters))

**📥 Outputs**  
• `Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.maker_out.gff3`: Initial NLR prediction results

⚙️ **Run MAKER**  
Default control files from MAKER will be used. In **maker_opts.ctl**, some parameters were modified as below. 
```
est=est.fa                        #est evidence
protein=proteins.fa               #protein evidence
snaphmm=Kronos.hmm                #paramters trained as part of ginger, v1.0 annotation
augustus_species=Kronos_collapsed #parameters trained as part of braker, v1.0 annotation
```

Due to some incompatibility between OpenMPI in our computer cluster and MAKER, each contig was separated and MAKER was run in parallel. 
```
#separate genome
mkdir split_genome
seqkit split -i -O split_genome Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa

#run maker for each directory
maker -RM_off -genome ${dir}.fa

#collect all outputs
for dir in $(ls -d *); do 
    gff3_merge -d ${dir}/${dir}.maker.output/${dir}_master_datastore_index.log -o ${dir}/${dir}.gff3
done

#combine into one file
cat */*.gff3 > Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.maker_out.gff3
```

---
### 3. Additional gene models 
Other than MAKER annotations, intermediate and final annotation files produced during the version 1 and 2 annotations were also included.  

**📥 Inputs**   
• `Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa`: NLR loci  
• `Kronos.v1.0.all.gff3`: version 1.0 annotation
• `Kronos.v2.0.all.gff3`: version 2.0 annotation
• `v1_abinitio.gff3`: annotations from BRAKER, Ginger and Funannotate produced as part of v 1.0 annotation. Re-coordinated file is available in this folder. 

**📥 Outputs**  
• `all_models.recoordinated.gff3`: Additional NLR gene models


⚙️ **Re-coordinate to match NLR loci**  
```
grep ">" Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa | cut -d ">" -f 2 | sort -u > coordinates.list 
python recoordinate_gff3.py Kronos.v1.0.all.gff3
python recoordinate_gff3.py Kronos.v2.0.gff3
python recoordinate_gff3.py v1_abinitio.gff3
cat Kronos.v1.0.all.recoordinated.gff3 Kronos.v2.0.recoordinated.gff3 v1_abinitio.recoordinated.gff3  > all_models.recoordinated.gff3
```

---
### 4. Transcriptome evidence

**📥 Inputs**   
• `Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa`: NLR loci  

**📥 Outputs**  
• `NLR.merged.bam`: transcriptome alignments
• `NLR.merged.bigwig`: transcriptome coverage

⚙️ **Run STAR**  
```
#index
STAR --runMode genomeGenerate \
     --genomeFastaFiles Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa \
     --genomeDir GenomeDir

#align transcriptome data from Kronos to the putative NLR loci
#Nnte that if coverage is high, the alignment may not be visuzlied in Apollo genome browser.
#use publicly available RNAseq data (reads.list)

while read -r read1 read2; do
    # Extract the prefix from the first read filename
    prefix=$(basename "$read1" | cut -d "_" -f 1)

    # Run STAR with the specified parameters
    STAR --runThreadN 56 \
        --genomeDir GenomeDir \
        --readFilesIn "$read1" "$read2" \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterMultimapNmax 3 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --outFileNamePrefix "${prefix}." \
        --readFilesCommand zcat
done < reads.list

#filter alignments
for bam in *.bam; do
        samtools view -F 260 -q 20 -@ 56 -b ${bam} > ${bam}.filtered.bam
done

#merge and sort
samtools merge -@ 56 NLR.merged.bam *.filtered.bam
samtools index -@ 56 NLR.merged.bam

#collapse the bam file to a bigwig file.
bamCoverage -b NLR.merged.bam -o NLR.merged.bigwig --binSize 1 --normalizeUsing None
```

---
### 4. Domain annotation

Lastly, domain prediction results are also useful. This information was generated as below.
```
#identify open reading frames (ORFs)
orfipy --procs 56 --bed orfs.bed --pep orfs.aa.fa --max 45000 --ignore-case --partial-3 --partial-5 Kronos.collapsed.chromosomes.masked.v1.1.fa

#predict pfam domains
global/scratch/users/skyungyong/Software/interproscan-5.68-100.0/interproscan.sh \
    -i orfs.aa.fa \
    --appl PFAM-37.0 \
    --disable-precalc \
    --output-dir ./output \
    --tempdir ./temp \
    --cpu 40 \
```

### 4. Manual curation

For manual curation, all these datasets need to be loaded into the Apollo Genome Browser. We conducted annotations for each chromosome. It took about 3 weeks to annotate all loci. After curation, there was a lot of downstream QC and manual examination over multiple iterations.

```
#create reference sequence
perl ./Apollo/bin/prepare-refseqs.pl --fasta ${chromosome}.nlr_loci.fa --out .

#load annotations
#first, extract genes in that chromosome
cat all_models.recoordinated.gff3 Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.maker_out.gff3 | \
    awk '$1 == ${chromosome} && ($3 == "gene" || $3 == "mRNA" || $3 == "exon" || $3 == "CDS") {print}' > ${chromosome}.genes.gff3

#get rid of models with in-frame stop codon using agat and then load
agat_sp_flag_premature_stop_codons.pl --gff ${chromosome}.genes.gff3 --fa ${chromosome}.nlr_loci.fa --out ${chromosome}.genes.fixed.gff3
perl ./Apollo/bin/flatfile-to-json.pl --trackLabel models --type mRNA --className mRNA --out . --gff ${chromosome}.genes.fixed.gff3

#change and load interproscan coordinates
python recoordinate_ipr.py orfs.aa.fa.gff3 ${chromosome}.nlr_loci.fa ${chromosome} > ${chromosome}.ipr.gff3
perl ./Apollo/bin/flatfile-to-json.pl --trackLabel iprscan --type match:iprscan --out . --gff ${chromosome}.ipr.gff3

#load transcriptome data
perl ./Apollo/bin/add-bam-track.pl --bam_url ${chromosome}.bam --label bam --in ./trackList.json
perl ./Apollo/bin/add-bw-track.pl --bw_url ${chromosome}.bam.bigwig --label bigwig --in ./trackList.json

#run apollo
./Apollo/bin/apollo run-local 7070
```



```
seqkit grep -f <(awk '$2 == "High" || $2 == "Medium" {print $1".1"}' NLR_confidence.list) Kronos.NLRs.fa > Kronos.NLRs.reliable.fa
[INFO] 1089 patterns loaded from file
 hmmsearch --domtblout Kronos.NLRs.reliable.against.PF00931.hmm.out --cpu 56 -E 1e-04 --domE 1e-04 PF00931.hmm Kronos.NLRs.reliable.fa

awk '{print $1}' Kronos.NLRs.reliable.against.PF00931.hmm.out | grep 'KRN' | sort -u | wc -l
1067

hmmalign --trim -o Kronos.NLRs.reliable.hmmalign.sto PF00931.hmm Kronos.NLRs.reliable.fa
esl-reformat fasta Kronos.NLRs.reliable.hmmalign.sto > Kronos.NLRs.reliable.hmmalign.fa
cd-hit -T 2 -c 0.8 -i Kronos.NLRs.reliable.hmmalign.fa -o Kronos.NLRs.reliable.hmmalign.cd-hit.filtered.fa
mafft --maxiterate 1000 --globalpair --thread 56 Kronos.NLRs.reliable.hmmalign.cd-hit.filtered.fa > Kronos.NLRs.reliable.hmmalign.cd-hit.filtered.mafft.msa.fa
trimal -gt 0.3 -in Kronos.NLRs.reliable.hmmalign.cd-hit.filtered.mafft.msa.fa -out Kronos.NLRs.reliable.hmmalign.cd-hit.filtered.mafft.msa.filtered.fa #246 positions left

python remove_gappy_seqs.py Kronos.NLRs.reliable.hmmalign.cd-hit.filtered.mafft.msa.filtered.fa Kronos.NLRs.reliable.hmmalign.cd-hit.filtered.mafft.msa.filtered.clean.fa
22 gappy sequences removed
hmmbuild -n Kronos_NBARC Kronos_NBARC.hmm Kronos.NLRs.reliable.hmmalign.cd-hit.filtered.mafft.msa.filtered.clean.fa

hmmsearch --domtblout Kronos.NLRs.reliable.against.Kronos_NBARC.hmm.out --cpu 56 -E 1e-04 --domE 1e-04 Kronos_NBARC.hmm Kronos.NLRs.reliable.fa
awk '{print $1}' Kronos.NLRs.reliable.against.Kronos_NBARC.hmm.out | grep 'KRN' | sort -u | wc -l
1077 

#there are small number of very divergent NBARC domain that cannot be detected using this method. We will likeley skip these

```

Phylogeny:
hmmsearch --domE 1e-4 -E 1e-4 --domtblout Kronos_and_known_NLRs.against.Kronos_NBARC.hmm.out --cpu 40 ../NLR_annotations/HMM/Kronos_NBARC.hmm Kronos_and_known_NLRs.fa
hmmalign --trim -o Kronos_and_known_NLRs.hmmalign.sto ../NLR_annotations/HMM/Kronos_NBARC.hmm Kronos_and_known_NLRs.fa
esl-reformat fasta Kronos_and_known_NLRs.hmmalign.sto >Kronos_and_known_NLRs.hmmalign.fasta
mafft --maxiterate 1000 --globalpair --thread 40 Kronos_and_known_NLRs.hmmalign.fasta > Kronos_and_known_NLRs.hmmalign.msa.fasta
trimal -gt 0.3 -in Kronos_and_known_NLRs.hmmalign.msa.fasta -out Kronos_and_known_NLRs.hmmalign.msa.filtered.fasta
python remove_gappy_seqs.py Kronos_and_known_NLRs.hmmalign.msa.filtered.fasta Kronos_and_known_NLRs.hmmalign.msa.filtered.clean.fasta
15 gappy sequences removed

less Kronos_and_known_NLRs.hmmalign.msa.filtered.clean.fasta | gre
p ">" | grep -c "KRN"
1074 #1121 total: non

raxml-ng --all \
  --msa Kronos_and_known_NLRs.hmmalign.msa.filtered.clean.fasta \
  --search 20 \
  --model LG+G8+F \
  --tree pars{10} \
  --bs-trees 1000 \
  --threads 56


## NLR Prediction in Wheat Genomes

```
augustus v3.3.3

```

In order to identify highly variable NLRs, enough coverages of homologous sequences are needed. We will collect them by predicting NLR genes in other published wheat genomes. Let's first train NLR specific prediction parameters. 

### 1. Trainning Augustus
After curating and labeling about 2,400 loci, we now know which genes are reliable and which genes have been pseudoginized. Using high-confidence genes, we will train prediction parameters.

```
#convert the nlr gff file to cds sequences. 
gffread -x Kronos_NLRs.final.cds.fa -g Kronos.collapsed.chromosomes.masked.v1.1.fa Kronos_NLRs.final.gff3

#get high-confidence NLRs and reduce redundancy 
#this step generates 828 high-confidence nlrs
seqkit grep -f <(awk '$2 == "High" {print $1".1"}' NLR_confidence.list) Kronos_NLRs.final.cds.fa > Kronos_hc.cds.fa
cd-hit-est -c 0.9 -i Kronos_hc.cds.fa -o Kronos_hc.cds.cd-hit.c_0.9.fa -T 30
grep ">" Kronos_hc.cds.cd-hit.c_0.9.fa | sed 's/>//g' > geneIDs.list

#create a separate gtf file for the selected high-confidence nlr genes
gffread -T --ids geneIDs.list ../NLR_final_datasets_curation/Kronos_all.NLRs.final.gff3  > 828_hc_nlrs.gtf
awk '$3 == "CDS" {print}' 828_hc_nlrs.gtf  > 828_hc_nlrs.cds.gtf

#create a genbank file with 1,000 bp flanking regions
perl gff2gbSmallDNA.pl 828_hc_nlrs.cds.gtf Kronos.collapsed.chromosomes.masked.v1.1.fa 1000 first.gb

#run quality check
#here, everything passes, but TrturKRN7A02G026020 got removed from the previous possibly as this gene is within an intron of another gene. 
etraining --species=generic first.gb > train.err

#separate the genes into 126 test set, and 700 trainning set
/global/scratch/users/skyungyong/Software/Augustus-3.3.3/augustus-3.3.3/scripts/perl randomSplit.pl first.gb 126

#train augustus 
perl /global/scratch/users/skyungyong/Software/Augustus-3.3.3/augustus-3.3.3/scripts/new_species.pl --species=Wheat_NLR #this is in maker environemnt 
etraining --species=Wheat_NLR first.gb.train
optimize_augustus.pl --cpus=40 --kfold=40 --species=Wheat_NLR first.gb.train
```

### 2. Performance Evaluation

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

So the NLR-specific parameters are much better!

We can also train NLR-specific parameters for SNAP as below. 

```
#use randomly selected 700 genes for training 
grep ">"  Kronos_hc.cds.cd-hit.c_0.9.fa  | cut -d ">" -f 2 | shuf | head -n 700 > selected.hc.list
gffread --ids selected.hc.list  828_hc_nlrs.gtf > selected.hc.gff
awk '$3 != "exon" {print}' selected.hc.gff  > selected.hc.cds.gff
agat_convert_sp_gff2zff.pl --gff selected.hc.cds.gff  --fasta Kronos.collapsed.chromosomes.masked.v1.1.folded.fa -o genome.ann
grep ">" genome.ann  | cut -d ">" -f 2 | while read line; do awk -v seq=$line -v RS=">" '$1 == seq {print RS $0; exit}' Kronos.collapsed.chromosomes.masked.v1.1.fa; done > genome.dna

fathom -validate genome.ann genome.dna
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl Kronos_NLRs . > Wheat_NLR.hmm
```

### 3. NLR Prediction in Wheat Genomes

NLRs were then predicted in wheat genomes. Genomes are listed in **genome.list**. 
```
#download genomes
while IFS=$'\t' read -r col1 prefix url; do
    mkdir -p "$prefix" && cd "$prefix" || continue
    wget "$url"
    cd ..
done < genome.list
```

Following the previous workflow, the NLR loci were detected
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
  python crop_genome.py --hmm orfs.aa.fa.against.NBARC.out --nlrannot NLRannotator.whole-genome.gff3 --genome Kronos.collapsed.chromosomes.masked.v1.1.fa

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
