# NLR Curation and Analyses


## NLR Curation

### 1. NLR Loci Isolation
The Kronos genome is large, and for manual curation, NLR loci need to be first extracted from the genome. Locate the region of genomes in which NB-ARC domains are detected, with 15,000 flanking sequences from both ends. These steps will generate an output named *Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa*. We later learned that some divergent NB-ARC domains in Kronos cannot be properly detected by this approach and additionally incoporated NLR-Annotator (step 7). Bringing NLR-Annotator early in this step can enhance the accuracy of loci detection.  

```
# Identify open reading frames (ORFs)
orfipy --procs 56 --bed orfs.bed --pep orfs.aa.fa --max 45000 --ignore-case --partial-3 --partial-5 Kronos.collapsed.chromosomes.masked.v1.1.fa

# Search for NB-ARC domains using HMMER
hmmsearch --cpu 56 --domE 1e-4 -E 1e-4 --domtblout orfs.against.NBARC.out PF00931.hmm orfs.aa.fa

# Crop genome to NB-ARC-containing loci
python crop_genome.py orfs.against.NBARC.out Kronos.collapsed.chromosomes.masked.v1.1.fa
```

### 2. Initial Annotations
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

### 3. Collecting Existing Annotations
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

### 4. Transcriptome Evidence
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
bamCoverage -b merged.bam -o output.bigwig --binSize 10 --normalizeUsing CPM
```

### 5

### 4. Manual curation

perl Apollo/bin/prepare-refseqs.pl --fasta 1A.fa --out .
for feature in {augustus,snap,maker,KRNv1.0,KRNv2.0,v1Annot}; do perl Apollo/bin/flatfile-to-json.pl --trackLabel ${feature} --type mRNA --className mRNA --out . --gff 1A.gff3; done

### 5. NLR-Annotator

After the initial curation, we run out first QC. In this step, we ensure that all putative NLR loci are detected and annotated. 

```
java -jar /global/scratch/users/skyungyong/Software/NLR-Annotator/NLR-Annotator-v2.1b.jar -t 40 -x /global/scratch/users/skyungyong/Software/NLR-Annotator/src/mot.txt -y /global/scratch/users/skyungyong/Software/NLR-Annotator/src/store.txt -i  /global/scratch/projects/vector_kvklab/KS-Kronos_remapping/Kronos.collapsed.chromosomes.masked.v1.1.fa -o NLRannotator.whole-genome.out -g NLRannotator.whole-genome.gff3
```


NLR prediction in other wheat species

gffread -x Kronos_NLRs.final.cds.fa -g ../../5.Annotations/Final/Final_Final_for_release/Kronos.collapsed.chromosomes.masked.v1.1.fa ../NLR_final_datasets_curation/Kronos_all.NLRs.final.gff3

seqkit grep -f <(awk '$2 == "High" {print $1".1"}' /global/scratch/projects/vector_kvklab/KS-Kronos_Final_datasets/03.NLRs/03.Final_NLR_datasets/NLR_confidence.list) Kronos_NLRs.final.cds.fa > Kronos_hc.cds.fa
cd-hit -c 0.9 -i Kronos_hc.cds.fa -o Kronos_hc.cds.cd-hit.c_0.9.fa -T 30
grep ">" Kronos_hc.cds.cd-hit.c_0.9.fa | sed 's/>//g' > geneIDs.list #827 nr-NLRs

gffread -T --ids geneIDs.list ../NLR_final_datasets_curation/Kronos_all.NLRs.final.gff3  > 827_hc_nlrs.gtf
awk '$3 == "CDS" {print}' 827_hc_nlrs.gtf  > 827_hc_nlrs.cds.gtf
perl /global/scratch/users/skyungyong/Software/Augustus-3.3.3/augustus-3.3.3/scripts/gff2gbSmallDNA.pl 827_hc_nlrs.cds.gtf ../../5.Annotations/Final/Final_Final_for_release/Kronos.collapsed.chromosomes.masked.v1.1.fa 1000 first.gb


etraining --species=generic first.gb > train.err #all genes are good; but TrturKRN7A02G026020 got removed from the previous step possibly as this gene is within an intron of another gene. 
perl /global/scratch/users/skyungyong/Software/Augustus-3.3.3/augustus-3.3.3/scripts/randomSplit.pl first.gb 126 #700 tranning 126 testing 

# Training 1
perl /global/scratch/users/skyungyong/Software/Augustus-3.3.3/augustus-3.3.3/scripts/new_species.pl --species=Wheat_NLR #this is in maker environemnt 
etraining --species=Wheat_NLR first.gb.train 
	
augustus --species=Wheat_NLR first.gb.test | tee firsttest.out
*******      Evaluation of gene prediction     *******

---------------------------------------------\
                 | sensitivity | specificity |
---------------------------------------------|
nucleotide level |       0.985 |       0.972 |
---------------------------------------------/

----------------------------------------------------------------------------------------------------------\
           |  #pred |  #anno |      |    FP = false pos. |    FN = false neg. |             |             |
           | total/ | total/ |   TP |--------------------|--------------------| sensitivity | specificity |
           | unique | unique |      | part | ovlp | wrng | part | ovlp | wrng |             |             |
----------------------------------------------------------------------------------------------------------|
           |        |        |      |                 78 |                 93 |             |             |
exon level |    325 |    340 |  247 | ------------------ | ------------------ |       0.726 |        0.76 |
           |    325 |    340 |      |   43 |    2 |   33 |   43 |    3 |   47 |             |             |
----------------------------------------------------------------------------------------------------------/

----------------------------------------------------------------------------\
transcript | #pred | #anno |   TP |   FP |   FN | sensitivity | specificity |
----------------------------------------------------------------------------|
gene level |   138 |   126 |   74 |   64 |   52 |       0.587 |       0.536 |
----------------------------------------------------------------------------/

------------------------------------------------------------------------\
            UTR | total pred | CDS bnd. corr. |   meanDiff | medianDiff |
------------------------------------------------------------------------|
            TSS |          4 |              0 |         -1 |         -1 |
            TTS |          9 |              0 |         -1 |         -1 |
------------------------------------------------------------------------|
            UTR | uniq. pred |    unique anno |      sens. |      spec. |
------------------------------------------------------------------------|
                |  true positive = 1 bound. exact, 1 bound. <= 20bp off |
 UTR exon level |          0 |              0 |       -nan |       -nan |
------------------------------------------------------------------------|
 UTR base level |          0 |              0 |       -nan |       -nan |
------------------------------------------------------------------------/

#braker-trained parameters 
augustus --species=Kronos_collapsed first.gb.test | tee firsttest.out
*******      Evaluation of gene prediction     *******

---------------------------------------------\
                 | sensitivity | specificity |
---------------------------------------------|
nucleotide level |       0.865 |       0.932 |
---------------------------------------------/

----------------------------------------------------------------------------------------------------------\
           |  #pred |  #anno |      |    FP = false pos. |    FN = false neg. |             |             |
           | total/ | total/ |   TP |--------------------|--------------------| sensitivity | specificity |
           | unique | unique |      | part | ovlp | wrng | part | ovlp | wrng |             |             |
----------------------------------------------------------------------------------------------------------|
           |        |        |      |                352 |                213 |             |             |
exon level |    479 |    340 |  127 | ------------------ | ------------------ |       0.374 |       0.265 |
           |    479 |    340 |      |  171 |   75 |  106 |  159 |   20 |   34 |             |             |
----------------------------------------------------------------------------------------------------------/

----------------------------------------------------------------------------\
transcript | #pred | #anno |   TP |   FP |   FN | sensitivity | specificity |
----------------------------------------------------------------------------|
gene level |   170 |   126 |    0 |  170 |  126 |           0 |           0 |
----------------------------------------------------------------------------/

------------------------------------------------------------------------\
            UTR | total pred | CDS bnd. corr. |   meanDiff | medianDiff |
------------------------------------------------------------------------|
            TSS |         12 |              0 |         -1 |         -1 |
            TTS |          7 |              0 |         -1 |         -1 |
------------------------------------------------------------------------|
            UTR | uniq. pred |    unique anno |      sens. |      spec. |
------------------------------------------------------------------------|
                |  true positive = 1 bound. exact, 1 bound. <= 20bp off |
 UTR exon level |          0 |              0 |       -nan |       -nan |
------------------------------------------------------------------------|
 UTR base level |          0 |              0 |       -nan |       -nan |
------------------------------------------------------------------------/



augustus --species=Kronos_manual --AUGUSTUS_CONFIG_PATH=/global/scratch/users/skyungyong/Kronos/5.Annotations/Braker/config first.gb.test | tee firsttest.out

---------------------------------------------\
                 | sensitivity | specificity |
---------------------------------------------|
nucleotide level |       0.851 |       0.935 |
---------------------------------------------/

----------------------------------------------------------------------------------------------------------\
           |  #pred |  #anno |      |    FP = false pos. |    FN = false neg. |             |             |
           | total/ | total/ |   TP |--------------------|--------------------| sensitivity | specificity |
           | unique | unique |      | part | ovlp | wrng | part | ovlp | wrng |             |             |
----------------------------------------------------------------------------------------------------------|
           |        |        |      |                287 |                156 |             |             |
exon level |    471 |    340 |  184 | ------------------ | ------------------ |       0.541 |       0.391 |
           |    471 |    340 |      |  158 |   33 |   96 |  120 |    3 |   33 |             |             |
----------------------------------------------------------------------------------------------------------/

----------------------------------------------------------------------------\
transcript | #pred | #anno |   TP |   FP |   FN | sensitivity | specificity |
----------------------------------------------------------------------------|
gene level |   161 |   126 |   34 |  127 |   92 |        0.27 |       0.211 |
----------------------------------------------------------------------------/

------------------------------------------------------------------------\
            UTR | total pred | CDS bnd. corr. |   meanDiff | medianDiff |
------------------------------------------------------------------------|
            TSS |         16 |              0 |         -1 |         -1 |
            TTS |          6 |              0 |         -1 |         -1 |
------------------------------------------------------------------------|
            UTR | uniq. pred |    unique anno |      sens. |      spec. |
------------------------------------------------------------------------|
                |  true positive = 1 bound. exact, 1 bound. <= 20bp off |
 UTR exon level |          0 |              0 |       -nan |       -nan |
------------------------------------------------------------------------|
 UTR base level |          0 |              0 |       -nan |       -nan |
------------------------------------------------------------------------/

So the NLR-specific parameters are much better!

optimize_augustus.pl --cpus=40 --kfold=40 --species=Wheat_NLR first.gb.train
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

------------------------------------------------------------------------\
            UTR | total pred | CDS bnd. corr. |   meanDiff | medianDiff |
------------------------------------------------------------------------|
            TSS |          4 |              0 |         -1 |         -1 |
            TTS |          7 |              0 |         -1 |         -1 |
------------------------------------------------------------------------|
            UTR | uniq. pred |    unique anno |      sens. |      spec. |
------------------------------------------------------------------------|
                |  true positive = 1 bound. exact, 1 bound. <= 20bp off |
 UTR exon level |          0 |              0 |       -nan |       -nan |
------------------------------------------------------------------------|
 UTR base level |          0 |              0 |       -nan |       -nan |
------------------------------------------------------------------------/
 
Here, our aim is to define highly variable NLR group within Kronos. Other wheat species may have divergent NLRs, whose the close homologs may be missing in Kronos. These genes may not be correctly predicted. This is OK. These genes will be filtered out anyways, as they will not offer any evolutionary information in the Shannon Entropy analyses of Kronos NLRs. 


