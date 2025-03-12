
# Protein-coding Gene Preidction: v1.0 annotation

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
```


The first version of genome annotation largely focused on the integration of short-read sequencing data produced for Kronos and the consensus of multiple gene prediction software

### 1. Paired-end Short-read Transcriptome Data Processing
```
inputs:
Publicly available RNA-seq data for Kronos

outputs:
all.merged.sorted.bam: filtered transcriptome alignments
transcripts.fasta: de novo and genome-guided transcript assemblies from trinity
stringtie.gtf: transcript assemblies from stringtie
sample_mydb_pasa.sqlite.assemblies.fasta: transcript assemblies from pasa
```

We donwloaded the paired-end RNA-seq data from the NCBI. The list can be found in **v1_rnaseq.list**. 
```
while read -r accession; do 
    sratoolkit.3.1.1-centos_linux64/bin/prefetch ${accession}
    sratoolkit.3.1.1-centos_linux64/bin/fasterq-dump -O . -e ${Numthreads} ${accession}
done < v1_rnaseq.list
```
Remove adapters and low-quality reads from the datasets, using trim_galore and cutadapt. This generated about 1.6 Tb of fastq files.
```
ls *.fastq | cut -d "_" -f 1 | sort -u | while read accession; do 
    trim_galore --paired -j 8 -a "${accession}_1.fastq" "${accession}_2.fastq"
done
```

The RNA-seq data were mapped to the genome by hisat, processed by samtools and assembled by stringtie. 
```
# index the genome (v1.0)
hisat2-build -p 20 Kronos.collapsed.chromosomes.fa Kronos

# align the reads
for lib in $(ls *_val_1.fq); do
  prefix=$(echo $lib | cut -d "_" -f 1)
  read_1=$lib
  read_2=$(echo $lib | sed 's/1_val_1/2_val_2/')
  hisat2 -p 56 -x Kronos -1 $read_1 -2 $read_2 --dta -S $prefix.mapped.sam
  echo 'done' > $prefix.done
done

# filter and sort with samtools
for sam in *.mapped.sam; do
  bam="${sam%.sam}.bam"
  samtools view -@ 56 -q 20 -h -b -F 260 "$sam" | samtools sort -@ 56 -o "$bam"
  samtools index "$bam"
done

#merge all bamfiles
samtools merge -@ 56 -h SRX10965366.mapped.bam -o all.merged.bam *.mapped.bam
samtools sort -@ 56 all.merged.bam > all.merged.sorted.bam

#assemble with stringtie
stringtie -o stringtie.gtf -p 56 --conservative all.merged.sorted.bam
```


De novo assembly was done using Trinity v2.15.1. We initially tried running Trinity on the 1.6 Tb of paired-end fastq files all at once. After two weeks, Trinity was still stuck at the insilico normalization step with about 35-45% progress. We, therefore, had to take some other ways around. Each pair will be normalized first, and then Trinity was run. This took a few days, producing transcripts of ~1 Gb. 
```
#for each pair 
singularity run trinity.sif Trinity --verbose --max_memory 90G --just_normalize_reads --seqType fq --CPU 40 --left $left --right $right --output trinity_$prefix

#list all normalized reads
ls -d trinity_* | while read folder; do
    prefix=$(echo $folder | cut -d "_" -f 2)
    left=$(ls $(pwd)/$folder\/insilico_read_normalization/*_1_val_1*.fq)
    right=$(ls $(pwd)\/$folder\/insilico_read_normalization/*_2_val_2*.fq)
    echo $prefix $prefix $left $right
done > sample.list

#run trinity de novo
singularity run trinity.sif Trinity --verbose --seqType fq --max_memory 1500G --CPU 56 --samples_file sample.list

#run trinity genome-guided
singularity run trinity.sif Trinity --verbose --max_memory 250G --CPU 56 --genome_guided_max_intron 10000 --genome_guided_bam all.merged.sorted.bam

#combine the oututs to transcripts.fa
#process assemblies
TransDecoder.LongOrfs -t transcripts.fa
TransDecoder.Predict -t TransDecoder.LongOrfs -t transcripts.fa
```

Finally, these assembled transcripts was processed by PASA, and orfs were predicted by transdecoder.
```
singularity exec pasapipeline.v2.5.3.simg /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl -c /usr/local/src/PASApipeline/sample_data/sqlite.confs/alignAssembly.config -r -C -R --CPU 56 --ALIGNERS gmap --TRANSDECODER -ALT_SPLICE -g Kronos.collapsed.chromosomes.masked.fa -t transcripts.fasta --trans_gtf stringtie.gtf

#process assemblies
TransDecoder.LongOrfs -t sample_mydb_pasa.sqlite.assemblies.fasta
TransDecoder.Predict -t TransDecoder.LongOrfs -t sample_mydb_pasa.sqlite.assemblies.fasta
```


### 2. BRAKER
```
inputs:
all.merged.sorted.bam: filtered transcritpome alignments produced using hisat and samtools
uniprotkb_38820.fasta: 2,850,097 protein sequences from Poales (TAXID: 38820) downloaded from UniProt

outputs:
braker.gtf: braker gene models
braker.aa: protein sequences of braker gene models
```

BRAKER was run as below.
```
singularity exec -B $PWD braker3.sif braker.pl --verbosity=3 \
    --genome=Kronos.collapsed.chromosomes.masked.fa \
    --bam=all.merged.sorted.bam \
    --prot_seq=uniprotkb_38820.fasta \
    --species=Kronos --threads 48 --gff3 \
    --workingdir=$wd/braker \
    --AUGUSTUS_CONFIG_PATH=$wd/config

#decorate utr
#utr is deleted during evidencemodler step, so this is not needed.
python stringtie2utr.py -g braker.gtf -s GeneMark-ETP/rnaseq/stringtie/transcripts_all.merged.sorted.gff -o braker.utr.gtf
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


blastp -query Kronos.v1.1.coordinate_fixed.pasa_updated.manual_fix.all_final.pep.fa -num_threads 56 -evalue 1e-10 -max_target_seqs 10 -max_hsps 1 -outfmt "6 std qlen slen" -out Kronos.v1.1.against.all.all_prot.max10 -db database/all_prot #pasa.orf #all_prot



(base) [skyungyong@ln002 MAKER]$ cd transcripts/
(base) [skyungyong@ln002 transcripts]$ ls
all.est.fasta  sample_mydb_pasa.sqlite.assemblies.fasta  stringtie.fasta
(base) [skyungyong@ln002 transcripts]$ cd ..
(base) [skyungyong@ln002 MAKER]$ ls
Runs  abinitio  proteins  transcripts
(base) [skyungyong@ln002 MAKER]$ rm -r transcripts/
(base) [skyungyong@ln002 MAKER]$ cd proteins/
(base) [skyungyong@ln002 proteins]$ ls
Aegilops_tauschii.Aet_v4.0.pep.all.fa                            Hordeum_vulgare.MorexV3_pseudomolecules_assembly.pep.all.fa  Triticum_aestivum.IWGSC.pep.all.fa            Triticum_urartu.IGDB.pep.all.fa
Avena_sativa_ot3098.Oat_OT3098_v2.pep.all.fa                     Hordeum_vulgare_goldenpromise.GPv1.pep.all.fa                Triticum_dicoccoides.WEWSeq_v.1.0.pep.all.fa  all.prot.evidence.fa
Avena_sativa_sang.Asativa_sang.v1.1.pep.all.fa                   Lolium_perenne.MPB_Lper_Kyuss_1697.pep.all.fa                Triticum_spelta.PGSBv2.0.pep.all.fa           get_pasa_orfs.py
Brachypodium_distachyon.Brachypodium_distachyon_v3.0.pep.all.fa  Secale_cereale.Rye_Lo7_2018_v1p1p1.pep.all.fa                Triticum_turgidum.Svevo.v1.pep.all.fa         pasa.transdecoder.pep.complete.fa
