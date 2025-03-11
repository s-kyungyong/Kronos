# NLR Curation and Analyses

## Data Availability

## Method

**Putative NLR loci detection**: to facilitate targeted NLR curation, we first identified, following a modified version of the previous workflow (Seong et al., 2020). The Kronos genome was translated in six frames, and open reading frames (ORFs) were identified using orfpy v0.0.4 (Singh and Wurtele, 2021). The predicted ORFs were searched for the NBARC domain, using the Hidden Markov Model obtained from PFAM (PF00931) and hmmsearch v3.4 (--domE 1e-4 -E 1e-4) (Eddy 2011; Mistry et al., 2021). While this approach captured the majority of putative NLR loci, some divergent NB-ARC domains remained undetected. To capture these additional loci, we also employed NLR-Annotator v2.1b (Steuernagel et al., 2020). All putative NLR loci detected by either approach were extracted from the Kronos genome with 15,000 flanking sequences on both sides. 

**Gene model prediction**: to support manual curation, initial gene models were predicted using MAKER v3.01.03 (Cantarel et al., 2008). Protein evidence was incorporated from NLR sequences of 18 Poaceae species (Toghani and Kamoun, 2024) and 415 reference NLRs from RefPlantNLR (Kourelis et al., 2021). EST evidence was derived from the transcripts assembled by Stringtie from short-read or long-read data, produced to create reference annotations v1.0 and v2.0, respectively. Additionally, two pre-trained ab initio prediction models developed during the version 1 annotation were used: Augustus parameters from BRAKER and SNAP parameters from GINGER. 

**Additional evidence**: To establish consensus, all intermediate and final gene models produced during reference genome annotations were incorporated. To support splicing site annotations, selected RNA-seq data were mapped to putative NLR loci using STAR v2.7.11b (Dobin et al., 2013), with coverage profiles generated using bamCoverage from deepTools v3.5.5 (Ramírez et al., 2014). To assist domain curation, PFAM domains (v36.0) were predicted using InterProScan v5.68.100 (Jones et al., 2014; Mistry et al., 2021). 

**Manual curation**: All predicted gene models, transcriptomic data, and domain annotations were loaded into Apollo Genome Browser v2.0.6 for manual curation (Dunn et al., 2019). Gene models containing an NB-ARC domain were systematically curated through visual inspection of gene structures and domain architectures. Each gene product was searched against the non-redundant database in The National Center for Biotechnology Information (NCBI) to assess sequence homology. During this process, each gene was assigned two classification labels: intact, partial or interrupted based on domain architecture and integrity, as well as consistent or inconsistent based on homology consistency (Fig. SX).

**Quality control**: In some cases, 30,000 flanking regions were insufficient to capture the full gene structure. To improve annotation accuracy, the draft NLR annotation sets were compared to the reference annotations to identify and rescue more complete gene models. The NLR annotations were loaded into Integrative Genome Browser v2.17.0 for further evaluation (Robinson et al., 2011). During the re-examination of gene models, each gene was assigned a confidence level. Genes with splicing sites supported by transcriptome data were classified as high-confidence. When transcriptome support was absent due to a lack of gene expression, gene models were evaluated based on structural conservation with close homologs in the NCBI database with available transcriptome evidence. Homology was reassessed to confirm that annotated splicing sites were not associated with gaps in pairwise alignments. Genes that contain mutations that disrupted splicing sites, introduced frameshifts, or caused premature termination were classified as low-confidence. Genes without strong supporting evidence but also without contradictions to their annotations were classified as medium-confidence.

---


## NLR Curation

### 1. Putative NLR loci detection
The Kronos genome is large, and for manual curation, NLR loci need to be first extracted from the genome. Locate the region of genomes in which NB-ARC domains are detected, with 15,000 flanking sequences from both ends. You may choose to increase the flanking size, as a small number of genes could not be fully contained in this region.

```
# Identify open reading frames (ORFs)
orfipy --procs 56 --bed orfs.bed --pep orfs.aa.fa --max 45000 --ignore-case --partial-3 --partial-5 Kronos.collapsed.chromosomes.masked.v1.1.fa

# Search for NB-ARC domains using HMMER
hmmsearch --cpu 56 --domE 1e-4 -E 1e-4 --domtblout orfs.against.NBARC.out PF00931.hmm orfs.aa.fa

# Crop genome to NB-ARC-containing loci
python crop_genome.py --hmm orfs.aa.fa.against.NBARC.out --genome Kronos.collapsed.chromosomes.masked.v1.1.fa
```


We later learned that some divergent NB-ARC domains in Kronos cannot be properly detected by this approach and additionally incoporated NLR-Annotator.
```
java -jar NLR-Annotator-v2.1b.jar -t 40 -x ./NLR-Annotator/src/mot.txt -y ./NLR-Annotator/src/store.txt -i Kronos.collapsed.chromosomes.masked.v1.1.fa -o NLRannotator.whole-genome.out -g NLRannotator.whole-genome.gff3

# Crop genome to NB-ARC-containing loci
python crop_genome.py --nlrannot NLRannotator.whole-genome.gff3 --genome Kronos.collapsed.chromosomes.masked.v1.1.fa
```

Alternatively, the two outputs can be combined.
```
python crop_genome.py --hmm orfs.aa.fa.against.NBARC.out --nlrannot NLRannotator.whole-genome.gff3 --genome Kronos.collapsed.chromosomes.masked.v1.1.fa
```

### 2. Gene model prediction
Initial gene models are needed to faciliate curation. Typically, if NLRs do not have any mutations that disrupt their gene structures (e.g. framshift mutations or mutations in splicing sites), evidence-based annotators and even ab initio annotators can correctly predict their gene structures (most of the time). Gene structures will be predicted with MAKER v3.01.03. For protein evidence, NLR sequences for 18 Poaceae species were collected from [this repository](https://zenodo.org/records/13627395) and 415 reference NLRs from [RefPlantNLR](https://zenodo.org/records/3936022). For EST evidence, the transcripts assemembled by Stringtie with short-reads (v1.0 annotation) and long-reads (v2.0 annotation) were used. The two ab initio parameters obtained in the v1 annotation were used: Augustus from BRAKER and SANP from GINGER. 

Default control files from MAKER will be used. In **maker_opts.ctl**, some parameters were modified as below. 
```
est=est.fa                        #est evidence
protein=proteins.fa               #protein evidence
snaphmm=Kronos.hmm                #paramters trained as part of ginger, v1.0 annotation
augustus_species=Kronos_collapsed #parameters trained as part of braker, v1.0 annotation
```

Due to some incompatibility between OpenMPI in our computer cluster and MAKER, each contig will be separated into a file and MAKER will run in parallel. 
```
#separate genome
mkdir split_genome
grep ">" Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa | cut -d ">" -f 2 | \
while read tig; do
  mkdir split_genome/${tig}
  awk -v seq=$tig -v RS=">" '$1 == seq {print RS $0; exit}' \
  Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa > split_genome/${tig}/${tig}.fa
done
```
Run MAKER. 
```
cd split_genome
for dir in $(ls -d *); do
  cd $dir
  maker -RM_off -genome ${dir}.fa 
  cd ..
done
```

After MAKER completes, collect all the annotations as gff files. 
```
for dir in $(ls -d *); do 
    gff3_merge -d ${dir}/${dir}.maker.output/${dir}_master_datastore_index.log -o ${dir}/${dir}.gff3
done
cat */*.gff3 > Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.maker_out.gff3
```

### 3. Additional evidence
Other than MAKER annotations, intermediate and final annotation files produced during the version 1 and 2 annotations are also included.
```
ls
Kronos.v1.0.all.gff3  #version 1.0 annotation: available through Zenodo
Kronos.v2.0.gff3      #version 2.0 annotation: available through Zenodo
v1_abinitio.gff3      #includes annotations from BRAKER, Ginger and Funannotate produced as part of the version 1 annotation
```

The coordinates in these GFF files need to be adjusted from the whole genome to the NLR-loci
```
less Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa | grep ">" | cut -d ">" -f 2 | sort -u > coordinates.list 
python recoordinate_gff3.py Kronos.v1.0.all.gff3
python recoordinate_gff3.py Kronos.v2.0.gff3
python recoordinate_gff3.py v1_abinitio.gff3
cat Kronos.v1.0.all.recoordinated.gff3 Kronos.v2.0.recoordinated.gff3 v1_abinitio.recoordinated.gff3  > all_models.recoordinated.gff3
```

Not all NLRs are well expressed. When they are, however, the transcriptome evidence can be useuful. Especially, when exon-intron structures are complicated, or when there are integrated domains nearby, this evidence can help improve the annotation. Let's generate this evidene track.

```
#index
STAR --runMode genomeGenerate \
     --genomeFastaFiles *Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa \
     --genomeDir GenomeDir

#align transcriptome data from Kronos to the putative NLR loci 
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

#merge all bam files
for bam in *.bam; do
        samtools view -F 260 -q 20 -@ 56 -b ${bam} > ${bam}.filtered.bam
done

#merge and sort
samtools merge -@ 56 merged.bam *.filtered.bam
samtools index -@ 56 merged.bam

#collapse the bam file to a bigwig file.
bamCoverage -b merged.bam -o output.bigwig --binSize 1 --normalizeUsing None
```


### 4. Manual curation

perl Apollo/bin/prepare-refseqs.pl --fasta 1A.fa --out .
for feature in {augustus,snap,maker,KRNv1.0,KRNv2.0,v1Annot}; do perl Apollo/bin/flatfile-to-json.pl --trackLabel ${feature} --type mRNA --className mRNA --out . --gff 1A.gff3; done

### 5. NLR-Annotator

After the initial curation, we run out first QC. In this step, we ensure that all putative NLR loci are detected and annotated. 

```
java -jar /global/scratch/users/skyungyong/Software/NLR-Annotator/NLR-Annotator-v2.1b.jar -t 40 -x /global/scratch/users/skyungyong/Software/NLR-Annotator/src/mot.txt -y /global/scratch/users/skyungyong/Software/NLR-Annotator/src/store.txt -i  /global/scratch/projects/vector_kvklab/KS-Kronos_remapping/Kronos.collapsed.chromosomes.masked.v1.1.fa -o NLRannotator.whole-genome.out -g NLRannotator.whole-genome.gff3
```

## NLR Prediction in Wheat Genomes

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
nucleotide level |       0.984 |        0.97 |
---------------------------------------------/

----------------------------------------------------------------------------------------------------------\
           |  #pred |  #anno |      |    FP = false pos. |    FN = false neg. |             |             |
           | total/ | total/ |   TP |--------------------|--------------------| sensitivity | specificity |
           | unique | unique |      | part | ovlp | wrng | part | ovlp | wrng |             |             |
----------------------------------------------------------------------------------------------------------|
           |        |        |      |                 79 |                 88 |             |             |
exon level |    331 |    340 |  252 | ------------------ | ------------------ |       0.741 |       0.761 |
           |    331 |    340 |      |   40 |    2 |   37 |   41 |    2 |   45 |             |             |
----------------------------------------------------------------------------------------------------------/

----------------------------------------------------------------------------\
transcript | #pred | #anno |   TP |   FP |   FN | sensitivity | specificity |
----------------------------------------------------------------------------|
gene level |   139 |   126 |   76 |   63 |   50 |       0.603 |       0.547 |
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

So the NLR-specific parameters are much better!




 
Here, our aim is to define highly variable NLR group within Kronos. Other wheat species may have divergent NLRs, whose the close homologs may be missing in Kronos. These genes may not be correctly predicted. This is OK. These genes will be filtered out anyways, as they will not offer any evolutionary information in the Shannon Entropy analyses of Kronos NLRs. 

while IFS=$'\t' read -r col1 prefix url; do
    mkdir -p "$prefix" && cd "$prefix" || continue
    wget "$url"
    cd ..
done < genome.list

```
for prefix in */; do
  cd "${prefix}" || exit
  
  gz=$(ls *.gz)
  genome="${gz%.gz}"
  gunzip "$gz"

  orfipy --procs 20 --bed orfs.bed --pep orfs.aa.fa --max 45000 --ignore-case --partial-3 --partial-5 "${genome}"

  hmmsearch --cpu 20 --domE 1e-4 -E 1e-4 --domtblout orfs.against.NBARC.out \
      /global/scratch/users/skyungyong/Kronos/NLR_annotations/Pan-NLRome/PF00931.hmm "${genome}_out/orfs.aa.fa"

  java -jar /global/scratch/users/skyungyong/Software/NLR-Annotator/NLR-Annotator-v2.1b.jar \
      -t 20 -x /global/scratch/users/skyungyong/Software/NLR-Annotator/src/mot.txt \
      -y /global/scratch/users/skyungyong/Software/NLR-Annotator/src/store.txt \
      -i "${genome}" -o NLRannotator.whole-genome.out -g NLRannotator.whole-genome.gff3
  done
  
  cd ..
done


python crop_genome.py --hmm orfs.aa.fa.against.NBARC.out --nlrannot NLRannotator.whole-genome.gff3 --genome Kronos.collapsed.chromosomes.masked.v1.1.fa
```


Due to some incompatibility between OpenMPI in our computer cluster and MAKER, each contig will be separated into a file and MAKER will run in parallel. 
```
#!/bin/bash

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
  cp /global/scratch/users/skyungyong/Kronos/NLR_annotations/Pan-NLRome/Evidence/maker* "$prefix"/
done

cd ../..

 
```
Run MAKER. 
```
cd split_genome
for dir in $(ls -d *); do
  cd $dir
  maker -RM_off -genome ${dir}.fa 
  cd ..
done
```
