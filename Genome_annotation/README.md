
# Protein-coding Gene Preidction: v1.0 annotation

## Data availability 
Protein coding genes can be downloaded from Zenodo.  
• `v1.0 annotation`: [https://zenodo.org/records/10215402](https://zenodo.org/records/10215402)  
• `v2.0 annotation`: [https://zenodo.org/records/14189805](https://zenodo.org/records/14189805)  
• `v2.1 annotation [Final]`: [https://zenodo.org/records/15539216](https://zenodo.org/records/14189805) 

## Software version
```
trim_galore v0.6.6
cutadapt v3.7
sratoolkit v3.1.1
hisat v2.2.1
samtools v1.9
stringtie v2.1.7
trinity v2.15.1
pasa v2.5.3
transdecoder v5.7.1
braker v3.0.6
funannotate v1.8.15
segmasker v1.0.0
augustus v3.3.2
snap v2013_11_29
ginger v1.0.1
miniprot v0.12
evidencemodeler v2.1.0
blast v2.15.0
```

---

## Annotation v1.0
### 1. Paired-end Short-read Transcriptome Data Processing

This step processes publicly available RNA-seq datasets for Kronos. Reads were downloaded, adapter-trimmed, aligned to the Kronos genome, and assembled both genome-guided and de novo to support gene structure prediction.

**📥 Inputs**  
• `v1_rnaseq.list`: List of NCBI SRA accessions  
• `Kronos.collapsed.chromosomes.fa`: Kronos reference genome  
• `Kronos.collapsed.chromosomes.masked.fa`: Kronos reference genome (masked)  

**📥 Outputs**  
• `all.merged.sorted.bam`: Merged and sorted RNA-seq alignments by mapping  
• `transcripts.fasta`: Trinity-assembled transcripts (de novo + genome-guided)   
• `stringtie.gtf`: Genome-guided transcript models from StringTie  
• `sample_mydb_pasa.sqlite.assemblies.fasta`: PASA-refined transcript structures  

---

**Download RNA-seq datasets from NCBI**  
```
while read -r accession; do 
    sratoolkit.3.1.1-centos_linux64/bin/prefetch ${accession}
    sratoolkit.3.1.1-centos_linux64/bin/fasterq-dump -O . -e ${Numthreads} ${accession}
done < v1_rnaseq.list
```

---

**Adapter trimming and quality filtering**  
```
ls *.fastq | cut -d "_" -f 1 | sort -u | while read accession; do 
    trim_galore --paired -j 8 -a "${accession}_1.fastq" "${accession}_2.fastq"
done
```
---
**Genome-guided mapping and transcript assembly**  
• Build genome index
```
hisat2-build -p 20 Kronos.collapsed.chromosomes.fa Kronos #v1.0 genome was used 
```
• Map RNA-seq reads
```
for lib in $(ls *_val_1.fq); do
  prefix=$(echo $lib | cut -d "_" -f 1)
  read_1=$lib
  read_2=$(echo $lib | sed 's/1_val_1/2_val_2/')
  hisat2 -p 56 -x Kronos -1 $read_1 -2 $read_2 --dta -S $prefix.mapped.sam
  echo 'done' > $prefix.done
done
```
• Sort and filter alignments
```
for sam in *.mapped.sam; do
  bam="${sam%.sam}.bam"
  #get primary alignments 
  samtools view -@ 56 -q 20 -h -b -F 260 "$sam" | samtools sort -@ 56 -o "$bam" 
  samtools index "$bam"
done
```
• Merge and assemble with StringTie
```
samtools merge -@ 56 -h SRX10965366.mapped.bam -o all.merged.bam *.mapped.bam #make sure to heve an header
samtools sort -@ 56 all.merged.bam > all.merged.sorted.bam

stringtie -o stringtie.gtf -p 56 --conservative all.merged.sorted.bam
```
---
**De novo and genome-guided transcriptome assembly with Trinity**  
Note: Due to memory limitations (>1.6 TB input), Trinity was run on individually normalized libraries.  

• Normalize each library
```
#run for for each pair 
singularity run trinity.sif Trinity --verbose --max_memory 90G --just_normalize_reads \
                                    --seqType fq --CPU 40 --left $left --right $right --output trinity_$prefix
```
• Create sample list of normalized reads
```
#list all normalized reads
ls -d trinity_* | while read folder; do
    prefix=$(echo $folder | cut -d "_" -f 2)
    left=$(ls $(pwd)/$folder\/insilico_read_normalization/*_1_val_1*.fq)
    right=$(ls $(pwd)\/$folder\/insilico_read_normalization/*_2_val_2*.fq)
    echo $prefix $prefix $left $right
done > sample.list
```
• Run Trinity (de novo and genome-guided)
```
#merge all outputs into transcripts.fasta
#run trinity de novo
singularity run trinity.sif Trinity --verbose --seqType fq --max_memory 1500G --CPU 56 --samples_file sample.list

#run trinity genome-guided
singularity run trinity.sif Trinity --verbose --max_memory 250G --CPU 56 --genome_guided_max_intron 10000 --genome_guided_bam all.merged.sorted.bam
```
---
**Transcript refinement with PASA**  
```
singularity exec pasapipeline.v2.5.3.simg /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \
            -c /usr/local/src/PASApipeline/sample_data/sqlite.confs/alignAssembly.config -r -C -R \
            --CPU 56 --ALIGNERS gmap --TRANSDECODER -g Kronos.collapsed.chromosomes.fa \
            -t transcripts.fasta \ #this includes all trinity-assembled transcripts (genome-guided and de novo)
            --trans_gtf stringtie.gtf #stringtie-based assembly outputs
```
---
### 2. Gene Prediction with BRAKER  
BRAKER was used to generate gene models using both RNA-seq alignment evidence and protein homology. Protein sequences from the Poales clade (TAXID: 38820) were downloaded from UniProt.  

**📥 Inputs**  
• `all.merged.sorted.bam`: Filtered transcriptome alignments (HISAT2 + SAMtools)  
• `uniprotkb_38820.fasta`: 2.85 million Poales proteins  

**📥 Outputs**  
• `braker.gtf`: Gene models predicted by BRAKER  

---

**Run BRAKER**  
```
singularity exec -B $PWD braker3.sif braker.pl --verbosity=3 \
    --genome=Kronos.collapsed.chromosomes.masked.fa \
    --bam=all.merged.sorted.bam \
    --prot_seq=uniprotkb_38820.fasta \
    --species=Kronos --threads 48 --gff3 \
    --workingdir=$wd/braker \
    --AUGUSTUS_CONFIG_PATH=$wd/config
```

### 3. Funannotate
```
inputs:
transcripts.fasta: de novo and genome-guided transcript assemblies from trinity
all.merged.sorted.bam: filtered transcritpome alignments produced using hisat and samtools
stringtie.gtf: transcript assemblies from stringtie
braker.gtf: braker gene models
braker.aa: protein sequences of braker gene models

outputs:
Triticum_kronos.filtered.gff3: funannotate gene models
```

Funannotate internally trains augustus and snap. However, we manually trained them and provided those parameters to funannotate. Gene models from BRAKER were searched against the IWGSC reference annotation (v1.0) and transcript assemblies from Trinity (protein sequences from transdecoder). Then, genes were selected if they had start and stop codons, having the same length with hits, sequence ideneity ≥ 99.5% and minimum protein lengths of 350. 6,000 genes were randomly selected for Augustus and SNAP, respectively. For Augustus, 5,700 were used as trainning set and 300 for testing. The worflow we followed is identical to the one described below in the GINGER section. 
```
blastp -query braker.aa -db trinity -max_target_seqs -max_hsps 1 -num_threads 56 -evalue 1e-10 -outfmt "6 std qlen slen" -out braker_vs_trinity.blast.out
blastp -query braker.aa -db iwgsc -max_target_seqs -max_hsps 1 -num_threads 56 -evalue 1e-10 -outfmt "6 std qlen slen" -out braker_vs_iwgsc.blast.out
cat  braker_vs_trinity.blast.out braker_vs_iwgsc.blast.out > braker.blast.out

#select genes
python select_genes_for_training.py 6000 
```

Funannotate was run as below. We made it to pick up the pre-trained SNAP parameters. 
```
funannotate predict \
-i Kronos.collapsed.chromosomes.masked.fa \
-o Funannotate \
-s "Triticum kronos" \
--transcript_evidence transcripts.fasta \
--repeats2evm \
--cpus 56 \
--ploidy 2 \ #2 was used as Kronos is homozygous and can be collapsed 
--rna_bam all.merged.sorted.bam \
--stringtie stringtie.gtf \
--augustus_species Kronos_maker \
--AUGUSTUS_CONFIG_PATH /global/scratch/users/skyungyong/Software/anaconda3/envs/funannotate/config/ \
--organism other \
--EVM_HOME /global/scratch/users/skyungyong/Software/anaconda3/envs/funannotate/opt/evidencemodeler-1.1.1/ \
--GENEMARK_PATH /global/scratch/users/skyungyong/Software/gmes_linux_64_4 \
```

Funannotate produced a lot of gene models than expected (137,303!). Gene models with high portion of low complexity regions were removed. This was not perfect filtering, but as different evidence would be combined at a later step, we did not dig in for further filtering. 
```
segmasker -in Triticum_kronos.proteins.fa -out Triticum_kronos.proteins.segmakser.out
python filter_genes_funannotate.py
```

### 4. GINGER
```
inputs:
transcripts.fasta: de novo and genome-guided transcript assemblies from trinity
all.merged.sorted.bam: filtered transcritpome alignments produced using hisat and samtools
left.norm.norm.fq: pre-normalized rnaseq data
right.norm.norm.fq: pre-normalized rnaseq data
braker.gtf: braker gene models
braker.aa: protein sequences of braker gene models

#protein evidence for SPALN (tritaest)
Triticum_aestivum.IWGSC.pep.all.fa
Triticum_turgidum.Svevo.v1.pep.all.fa
Triticum_dicoccoides.WEWSeq_v.1.0.pep.all.fa
Triticum_spelta.PGSBv2.0.pep.all.fa
Triticum_urartu.IGDB.pep.all.fa

outputs:
ginger_phase2.gff: ginger gene models
```

GINGER uses Nextflow to streamline annotations. We made a few changes here in the workflow
```
denovo.nf: transcript assembles with oases/velvet were not performed, as this required 17,000,000 Gb memory. 
abinitio.nf: Augustus and SNAP were trained manually
```

Augustus and SNAP were trainned similarly. Gene models from BRAKER were searched against the IWGSC reference annotation (v1.0) and transcript assemblies from Trinity. Then, genes were selected if they had start and stop codons, having the same length with hits, sequence ideneity ≥ 99.5% and minimum protein lengths of 350. This time, 8,500 genes were randomly selected for Augustus and SNAP, respectively. For Augustus, 5,700 were used as trainning set and 300 for testing. 
```
blastp -query braker.aa -db trinity -max_target_seqs -max_hsps 1 -num_threads 56 -evalue 1e-10 -outfmt "6 std qlen slen" -out braker_vs_trinity.blast.out
blastp -query braker.aa -db iwgsc -max_target_seqs -max_hsps 1 -num_threads 56 -evalue 1e-10 -outfmt "6 std qlen slen" -out braker_vs_iwgsc.blast.out
cat  braker_vs_trinity.blast.out braker_vs_iwgsc.blast.out > braker.blast.out

#select genes
python select_genes_for_training.py 8500 
```

AUGUSTUS was trained as below.
```
gff2gbSmallDNA.pl augustus.gff3 Kronos.collapsed.chromosomes.masked.fa 2000 genes.gb
randomSplit.pl genes.gb 400 # 400 test set
new_species.pl --species=Kronos_manual
etraining -species=Kronos_manual genes.gb.train
optimize_augustus.pl --species=Kronos_manual --cpus=48 --UTR=off genes.gb.train
```

SNAP was trained as below.
```
gff3_to_zff.pl genome.dna snap.gff3 > genome.ann # genome.dna = Kronos.collapsed.chromosomes.masked.fa
fathom -validate genome.ann genome.dna 
fathom -categorize 1000 genome.ann genome.dna 
fathom -export 100 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl Kronos . > Kronos_manual.hmm
```

GINGER was run roughly as below. Note that we had to make changes in the next flow modules for the modifications mentioned earlier in this section.
```
nextflow -C nextflow.config run mapping.nf
nextflow -C nextflow.config run denovo.nf
nextflow -C nextflow.config run abinitio.nf
nextflow -C nextflow.config run homology.nf
phase0.sh nextflow.config
phase1.sh nextflow.config > phase1.log
phase2.sh 50
summary.sh nextflow.config
```


### 5. Miniprot
```
inputs:
uniprotkb_38820.fasta: 2,850,097 protein sequences from Poales (TAXID: 38820) downloaded from UniProt

outputs:
miniprot.gff3: protein evidence alignments
```

We used miniprot to align the protein sequences from UniProt. 
```
miniprot -t 56 --gff --outc=0.95 -N 0 Kronos.collapsed.chromosomes.fa uniprotkb_38820.fasta > miniprot.gff3
```



### 6. EvidenceModler
```
inputs:
braker.gff: braker gene models
Triticum_kronos.filtered.gff3: funannotate gene models
ginger_phase2.gff: ginger gene models
miniprot.gff3: protein evidence alignments
sample_mydb_pasa.sqlite.pasa_assemblies.gff3: pasa transcript assemblies
sample_mydb_pasa.sqlite.assemblies.fasta.transdecoder.genome.gff3: pasa transcript assembles translated by transdecoder

outputs:
Kronos.EVM.gff3: evidencemodeler gene models
```

Finally, all evidence was combined by EvidenceModeler. There were some modifications within the inputs.
```
abinitio.gff3: this combines braker.gff, Triticum_kronos.filtered.gff3, ginger_phase2.gff and sample_mydb_pasa.sqlite.assemblies.fasta.transdecoder.genome.gff
homology.gff3: this includes miniprot.gff3 and an intermediate output from ginger (spaln alignments) 
transcripts.gff3: this includes sample_mydb_pasa.sqlite.pasa_assemblies.gff3
```

EvidenceModeler was run as below with the specified weights. 
```
singularity exec EVidenceModeler.v2.1.0.simg EVidenceModeler --sample_id Kronos --genome genome.fa --weights weights.txt --gene_predictions abinitio.gff3 --protein_alignments homology.gff3 --transcript_alignments transcripts.gff3 --repeats repeat.gff3 --CPU 56 -S --segmentSize 100000 --overlapSize 10000

#weight.txt
ABINITIO_PREDICTION     funannotate   3       #funannotat gene models
ABINITIO_PREDICTION     braker  3             #braker gene models
ABINITIO_PREDICTION     ginger  3             #ginger gene models
ABINITIO_PREDICTION     gingers  3.5          #single-exon genes from GINGER
PROTEIN homology        2                     #ginger's SPALN
PROTEIN                  miniprot       1     #protein mapping with miniprot
OTHER_PREDICTION        transdecoder    2.5   #pasa transcript assembles translated by transdecoder
TRANSCRIPT               pasa  8              #pasa transcript assemblies
```

### 7. PASA
```
Inputs:
Kronos.EVM.gff3: evidencemodeler gene models

Outputs:
Kronos.EVM.pasa.gff3: evidencemodeler gene models updated by pasa
```

PASA was run one more time to update UTRs and isoforms. 
```
singularity exec pasapipeline.v2.5.3.simg /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl -C -R -c alignAssembly.config -g Kronos.collapsed.chromosomes.masked.v1.1.fa -t transcripts.fasta --trans_gtf stringtie.gtf --TRANSDECODER --ALT_SPLICE --ALIGNERS gmap
singularity exec pasapipeline.v2.5.3.simg /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl -A -L -c compare.config  -g Kronos.collapsed.chromosomes.masked.v1.1.fa -t transcripts.fasta --annots Kronos.EVM.gff3
```

### 8. Gene Model Selection
```
Inputs:
Kronos.EVM.pasa.gff3: evidencemodeler gene models updated by pasa
pasa.transdecoder.pep.complete.fa: pasa transcript assemblies translated with transdecoder. These all have start and stop codons. 

protein databases:
Aegilops_tauschii.Aet_v4.0.pep.all.fa
Hordeum_vulgare.MorexV3_pseudomolecules_assembly.pep.all.fa
Triticum_aestivum.IWGSC.pep.all.fa
Triticum_urartu.IGDB.pep.all.fa
Avena_sativa_ot3098.Oat_OT3098_v2.pep.all.fa                    
Hordeum_vulgare_goldenpromise.GPv1.pep.all.fa
Triticum_dicoccoides.WEWSeq_v.1.0.pep.all.fa  all.prot.evidence.fa
Avena_sativa_sang.Asativa_sang.v1.1.pep.all.fa
Lolium_perenne.MPB_Lper_Kyuss_1697.pep.all.fa
Triticum_spelta.PGSBv2.0.pep.all.fa
Brachypodium_distachyon.Brachypodium_distachyon_v3.0.pep.all.fa
Secale_cereale.Rye_Lo7_2018_v1p1p1.pep.all.fa
Triticum_turgidum.Svevo.v1.pep.all.fa

Outputs:
Kronos.EVM.pasa.gff3: evidencemodeler gene models updated by pasa
```
Lastly, high-confidence genes were selected based on the evidence. 

```
#search against the protein databases
blastp -query Kronos.EVM.pasa..pep.fa -num_threads 56 -evalue 1e-10 -max_target_seqs 10 -max_hsps 1 -outfmt "6 std qlen slen" -out Kronos.v1.0.against.all.all_prot.max10 -db all_prot 
blastp -query Kronos.EVM.pasa..pep.fa -num_threads 56 -evalue 1e-10 -max_target_seqs 10 -max_hsps 1 -outfmt "6 std qlen slen" -out Kronos.v1.0.against.all.pasa_orf.max10 -db pasa_orf

#create version 1 annotations
#the input gff has to be correctly sorted first
python generate_v1.0_annot.py
```

----

# Protein-coding Gene Preidction: v2.0 annotation

```
sratoolkit v3.1.1
minimap v2.28-r1209
samtools v1.20
stringtie v2.1.7
isoquant v3.5.2
```

The second version of genome annotation integrates publicly available long-read sequencing data for Triticum.

### 1. Long-read Transcriptome Data Download

We donwloaded the long-read RNA-seq data from the NCBI. The list can be found in **v2_rnaseq.list**. 
```
while read -r accession; do 
    sratoolkit.3.1.1-centos_linux64/bin/prefetch ${accession}
    sratoolkit.3.1.1-centos_linux64/bin/fasterq-dump -O . -e ${Numthreads} ${accession}
done < v2_rnaseq.list
```

### 2. Alignment

These long-reads were aligned to the genome.
```
for fq in *.fq; do
  prefix=$(echo $fq | cut -d "." -f 1)
  minimap2 -I 12G -t 56 -x splice:hq -a -o "${prefix}.sam" Kronos.collapsed.chromosomes.masked.v1.1.broken.fa ${fq}
  samtools view -@56 -h -b "${prefix}.sam" | samtools sort -@ 56 > "${prefix}.bam"
done

#merge all bam files
samtools merge -@56 -h ERR11193282.bam all.long-read.merged.bam *.bam
samtools index -@56 all.long-read.merged.bam

#separate for each chromosome
for chromosome in \
  1A 1A_296625417 1B 1B_362283996 2A 2A_389606086 2B 2B_416081101 \
  3A 3A_354343362 3B 3B_427883679 4A 4A_376933649 4B 4B_351648618 \
  5A 5A_305547233 5B 5B_360298581 6A 6A_294206980 6B 6B_365632995 \
  7A 7A_370147894 7B 7B_378890030 Un; do
  samtools view -@56 -h -b all.long-read.merged.bam $chromosome > all.long-read.merged.${chromosome}.bam
done
```

### 3. Assembly

Transcripts were then assemblied, using stringtie and isoquant. The short-read transcriptome alignments obtained during the first version of annotation was included to enhance splicing site detection.
```
#de novo assembly using stringtie
stringtie -p 4 -v -o Kronos.${chromosome}.stringtie.denovo.gtf --mix all.short-read.merged.${chromosome}.bam all.long-read.merged.${chromosome}.bam

#reference-annotation guided using stringtie
stringtie -p -4 -o Kronos.${chromosome}.stringtie.guided.gtf -p 4 -G Kronos.v1.0.all.gff3 --mix all.short-read.merged.${chromosome}.bam all.long-read.merged.${chromosome}.bam

#assemblies using isoquant
isoquant.py --threads 56 --reference Kronos.collapsed.chromosomes.masked.v1.1.broken.fa --illumina_bam all.short-read.merged.bam --output Isoquant_Kronos --data_type pacbio_ccs --bam $bam
```

### 4. Updating Annotations



# Protein-coding Gene Preidction: v2.1 annotation

In this version, existing genes in the v2.0 annotation that overlap with low-confidence NLRs were discarded. High and medium-confidence NLRs were added instead. 
```
#input:
v2.0 annotations
nlr annotations
nlr confidence (high/medium/low)

#generate v2.1 annotation
generate_v2.1_annot.py

#sort gff
grep -v '^#' Kronos.v2.1.initial.gff3 | \
awk 'BEGIN {OFS="\t"} { 
    if ($3 == "gene") type = 1;
    else if ($3 == "mRNA") type = 2;
    else if ($3 == "exon") type = 3;
    else if ($3 == "CDS") type = 4;
    else type = 5;
    print $0, type 
}' | sort -k1,1 -k4,4n -k10,10n | cut -f1-9 > Kronos.v2.1.gff3  # Sort by chromosome, start, then type order

```


# None-coding RNA Preidction: 

```
cmscan v1.1.5
```


## RFAM search


```
#download the databases: rfam v15.0
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
cmpress Rfam.cm

#preparezvalues for each chromosome
cat zvalue.list
1A 1200.879162
1B 1417.658372
2A 1591.626378
2B 1657.049066
3A 1518.249256
3B 1728.269574
4A 1535.707034
4B 1399.361912
5A 1440.551718
5B 1462.269252
6A 1248.596746
6B 1467.16489
7A 1506.943932
7B 1532.01879
Un 421.039088

#for each chromosome
while read -r chromosome zvalue; do
  cmscan -Z ${zvalue} --cut_ga --rfam --nohmmonly --tblout ${chromosome}.Rfam.tblout --fmt 2 --cpu 56 --clanin Rfam.clanin Rfam.cm Kronos.v1.1.${chromosome}.fa
done < seqLengths.list

./tRNAscan-SE_installed/bin/tRNAscan-SE -E -o tRNAscan-SE.out -f tRNAscan-SE.ss -s tRNAscan-SE.iso -m tRNAscan-SE.stats -c ./tRNAscan-SE_installed/bin/tRNAscan-SE.conf ../Final/Kronos.collapsed.chromosomes.v1.1.fa
