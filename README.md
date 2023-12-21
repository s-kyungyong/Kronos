# Kronos Genomics

This repository describes the computational pipelines we used to build and analyze our Kronos genome (*Triticum durum* cv Kronos). Kronos is an allotetraploid (AABB).

## Genome assembly

### Quality control

#### HiFi reads

We have the following HiFi reads as bam files from Revio. In this step, the bamfiles are converted into a fastq file, and the reads are filtered. 

```
ls -lha m84066_2305*
-rwxr-xr-x 1 skyungyong ucb  45G Jun 12 20:07 m84066_230503_201048_s2.hifi_reads.default.bam
-rwxr-xr-x 1 skyungyong ucb  69M Jun 12 20:07 m84066_230503_201048_s2.hifi_reads.default.bam.pbi
-rwxr-xr-x 1 skyungyong ucb 436G Jun 12 21:55 m84066_230508_174311_s3.hifi_reads.default.bam
-rwxr-xr-x 1 skyungyong ucb  87M Jun 12 21:56 m84066_230508_174311_s3.hifi_reads.default.bam.pbi
-rwxr-xr-x 1 skyungyong ucb 422G Jun 12 23:52 m84066_230508_181330_s4.hifi_reads.default.bam
-rwxr-xr-x 1 skyungyong ucb  83M Jun 12 23:52 m84066_230508_181330_s4.hifi_reads.default.bam.pbi
-rwxr-xr-x 1 skyungyong ucb 312G Jun 13 01:14 m84066_230512_175630_s1.hifi_reads.default.bam
-rwxr-xr-x 1 skyungyong ucb  66M Jun 13 01:15 m84066_230512_175630_s1.hifi_reads.default.bam.pbi
-rwxr-xr-x 1 skyungyong ucb 346G Jun 13 02:54 m84066_230517_182502_s1.hifi_reads.default.bam
-rwxr-xr-x 1 skyungyong ucb  67M Jun 13 02:54 m84066_230517_182502_s1.hifi_reads.default.bam.pbi
-rwxr-xr-x 1 skyungyong ucb 341G Jun 13 04:38 m84066_230523_181907_s2.hifi_reads.default.bam
-rwxr-xr-x 1 skyungyong ucb  69M Jun 13 04:38 m84066_230523_181907_s2.hifi_reads.default.bam.pbi
```

Convert the bam files into a single fasta file with bam2fastq.

```
bam2fastq --version
bam2fastq 3.0.0 (commit v3.0.0)

Using:
  pbbam     : 2.3.0 (commit v2.3.0-3-g502c299)
  pbcopper  : 2.2.0 (commit v2.2.0-9-gde47bd7)
  boost     : 1.77
  htslib    : 1.15
  zlib      : 1.2.11
```
```
reads=$(ls *.default.bam)
bam2fastq -o Kronos.HiFi -j 52 $reads
```

Filter the reads with hifiadapterfilt. This step is not necessary. Often, the HiFi reads come out pretty clean. However, we occasionally had cases in other species where not-so-clean HiFi reads resulted in adapters contained inside the contigs. We thus included this step for the quality control. 
```
/HiFiAdapterFilt/hifiadapterfilt.sh -p Kronos -t 54
```

As can be seen from the statistics below, possible contaminants are extremely rare. 

```
cat Kronos.HiFi.stats
Started on Wed Jun 14 21:12:00 PDT 2023
For the Kronos.HiFi dataset:
Removing reads containing adapters a minimum of 44 bp in length and 97% match.

Number of ccs reads: 39390905
Number of adapter contaminated ccs reads: 13 (3.30025e-05% of total)
Number of ccs reads retained: 39390892 (100% of total)
```

#### Hi-C data

We have the following paired-end libraries from the Hi-C protocol. 

```
ls -lha  KVK-KRONOS-*/*.fastq.gz
-rwxr-xr-x 1 skyungyong ucb 49G Jul 26 11:06 KVK-KRONOS-874809/kvk-kronos-874809_S3HiC_R1.fastq.gz
-rw-r--r-- 1 skyungyong ucb 50G Jul 26 10:46 KVK-KRONOS-874809/kvk-kronos-874809_S3HiC_R2.fastq.gz
-rwxr-xr-x 1 skyungyong ucb 55G Jul 26 11:21 KVK-KRONOS-HIC2-1066975/1066975_S1_R1_001.fastq.gz
-rwxr-xr-x 1 skyungyong ucb 56G Jul 26 11:37 KVK-KRONOS-HIC2-1066975/1066975_S1_R2_001.fastq.gz
-rwxr-xr-x 1 skyungyong ucb 60G Jul 26 12:27 KVK-KRONOS-HIC2-1067017/kvk-kronos-hic2-1067017_S3HiC_R1.fastq.gz
-rwxr-xr-x 1 skyungyong ucb 62G Jul 26 12:01 KVK-KRONOS-HIC2-1067017/kvk-kronos-hic2-1067017_S3HiC_R2.fastq.gz
```

We use fastp to filter the reads. 

```
fastp --version
fastp 0.23.2

fastqc -version
FastQC v0.12.1
```

We will allow auto-detection of the adapters (--detect_adapter_for_pe) and trim polyG trimming (-g). The reads that have less than 20 quality score and are shorter than 50 nucleotides will be removed.

```
for folder in $(ls -d  KVK-KRONOS-*)
  do

  date
  echo STARTED: $folder
  cd $folder

  read1=$(ls *_R1*.fastq.gz)
  read2=$(ls *_R2*.fastq.gz)

  out1=$(ls *_R1*.fastq.gz | sed 's/fastq.gz/trimmed.fq.gz/g')
  out2=$(ls *_R2*.fastq.gz | sed 's/fastq.gz/trimmed.fq.gz/g')

  fastp --thread 54 -g --detect_adapter_for_pe -q 20 -l 50 --in1 $read1 --in2 $read2 --out1 $out1 --out2 $out2 -h bssh1.html &> bssh1.log
  fastqc -t 54 $out1 $out2

  date
  echo DONE: $folder
  cd ..

done
```

----------

## Pre-assembly assessment

Kronos is an allotetraploid species (AABB). In the field, Kronos rarely out-crosses. We therefore think that the heterozygosity should be extremely low or residual, and we aim to generate collapsed haplotypes (AB). An available Durum wheat genome ([Svevo](https://www.nature.com/articles/s41588-019-0381-3)) is 10.45G in size. We also roughly estimate that each haplotype will be 5-6G in size, totalling up to 10-12G for the collapsed haplotypes. We will evaluate our assumptions with GenomeScope. 

```
jellyfish --version
jellyfish 2.2.10

genomescope2 --version
GenomeScope 2.0
```

```
jellyfish count -C -m 21 -s 50000000000 -t 20 Kronos.HiFi.fastq -o kmer_counts.jf
jellyfish histo -h 5000000 -t 20 kmer_counts.jf > reads.histo
genomescope2 -p 4 -i reads.histo -o genomescope --verbose  
```

We can also perform similar analysis for Svevo. We first downloaded the paired-end short reads from the NCBI, filtered them and evaluated k-mer. 
```
fasterq-dump-orig.2.11.2 --version
fasterq-dump-orig.2.11.2 : 2.11.2

trim_galore --version
version 0.6.6

cutadapt --version
3.7
```

```
cat Svevo_SRA.list | while read accession; do fasterq-dump-orig.2.11.2 -e 40 -t ./tmp $accession; done
ls *.fastq | cut -d "_" -f 1 | sort -u |  while read accession; do trim_galore --illumina -j 8 --paired $accession\_1.fastq $accession\_2.fastq ; done

jellyfish count -C -m 21 -s 50000000000 -t 56 *.fq -o svevo.kmer_counts.jf
jellyfish histo -h 5000000 -t 56 svevo.kmer_counts.jf > svevo.reads.histo
genomescope2 -p 4 -i svevo.reads.histo -o svevo.genomescope --verbose
```

This is the GenomeScope statistics.  

|    | Kronos             ||              Svevo ||
|----|---------|-----------|---------|-----------|
|    | min | max | min | max |
| Homozygous (aaaa) | 88.853%  | 91.6029% |  89.1867%     |     91.4669% |
| Heterozygous (not aaaa)  | 8.39712% | 11.147% | 8.5331%    |       10.8133% |
| aaab | 0%  | 0.632431% | 0%      |          0.520259%|
| aabb | 2.87586% | 3.34577%  | 2.64974%    |      3.05399%|
| aabc | 1.65953%  | 2.9789% | 1.96852%     |     3.04119%|
| abcd | 3.86173% | 4.1899%  | 3.91484%     |     4.19786%|
| Genome Haploid Length | 2,605,665,036 bp | 2,612,433,029 bp  |  2,749,147,211 bp | 2,753,880,954 bp |
| Genome Repeat Length | 2,357,672,871 bp | 2,363,796,726 bp  |  2,493,441,344 bp | 2,497,734,789 bp |
|Genome Unique Length       |   247,992,165 bp  |  248,636,303 bp | 255,705,867 bp   | 256,146,166 bp |
|Model Fit                 |    28.7973%       |   89.812% | 26.1724%    |      82.8801% |
|Read Error Rate            |   0.0707465%     |   0.0707465% | 0.190909%    |     0.190909% |

The statistics looks quite similar!


We can compare our Kronos statistics to the GenomeScope result for the hexaploid wheat in [this paper](https://www.nature.com/articles/s41467-020-14998-3). The data is in Fig. S21. Here, the genome size is esimated as haplotype size x ploidy. Although this isn't technically correct, it seems to give the right estimate. See the discussion [here](https://github.com/schatzlab/genomescope/issues/107).

|    | Kronos | Svevo| Bread wheat |
|----|---------|-----------| -----------|
| ploidy (p) | 4 | 4 | 6 | 
| Haplotype size (Gb) | 2.612  | 2.754 | 2.354 |
| Genome size (Gb)  | 10.45 | 11.02 | 14.1 |
| unique (%) | 9.52 | 9.3 | 8 |
| heterozygosity (%) | 9.5 | 9.4 |  10.1 |
| err (%) | 0.0707 | 0.191 | 0.506 |
| dup | 0.654 | 2.67 | 0.836 | 

The statistics looks fairly similar. I believe that just like other reference wheat genomes, we can generate collapsed haplotypes (AB) for our Kronos genome. 

## Genome assembly and assessment

Genome assembly is done with hifiasm. Because the residual heterozygosity is low, and we aim to generate collapsed haplotypes (AB), we will only use the HiFi reads at the assembly stage.

```
hifiasm --version
0.19.5-r587
```

```
hifi=Kronos.HiFi.filt.fastq.gz
hifiasm -l0 -t 54 -o l0-hic $hifi
```

This is the very last line of the log. It took about 61 hours with 54 CPUs and a peak memory of 615 GB. 
```
[M::main] Real time: 218150.098 sec; CPU: 9384831.006 sec; Peak RSS: 615.502 GB
```

We can quickly look at the assembly statistics with QUAST and compleasm. 

```
awk '/^S/{print ">"$2"\n"$3}' l0.bp.p_ctg.gfa | fold > Kronos.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' l0.bp.a_ctg.gfa | fold > Kronos.a_ctg.fa
quast.py --fast -t 20 -o quast Kronos.p_ctg.fa Kronos.a_ctg.fa

VERSION=0.2.2
export SINGULARITY_CACHEDIR=/global/scratch/users/skyungyong/temp
singularity exec docker://huangnengcsu/compleasm:v${VERSION} compleasm run -l poales_odb10 -t 20 -o p_busco -a Kronos.p_ctg.fa
 exec docker://huangnengcsu/compleasm:v${VERSION} compleasm run -l poales_odb10 -t 20 -o a_busco -a Kronos.a_ctg.fa
```

Here is the statistics. 

|    | p_ctg | a_ctg |
|----|---------|-----------|
| # contigs  | 3317  | 2530 |
| Length (Gb)  | 10.55 | 0.074 |
| Largest contig (Mb) | 346.20  | 0.976 |
| N50 (Mb) | 40.96 | 0.0306  |
| L50| 50 | 6  | 849 |
||            ||
| Complete | 5.43% | 0.71% |
| Duplicated | 94.40% | 0.16% |
| Fragmented | 0.10% | 0.25% |
| Missing | 0.06% | 98.88% |

The associate contigs (a_ctg) include a lot of fragments that are potentially not useful. Most of these are likeley plasmids or repeats, which will be later discarded. Some might have been separated based on residual hetrozygosity. For now, let's combine the primary and associate contigs into a single file. 

```
cat Kronos.p_ctg.fa Kronos.a_ctg.fa > Kronos.draft.fa
```

Some may be interested in haplotypes, so we can also generate haplotype-resolved assemblies with the HiFi and HiC data together. However, as we focus on the collapsed genome, I will record the pipeline [here](https://github.com/s-kyungyong/Kronos/blob/main/Haplotypes/README.md).




## Scaffolding and assessment

Now, we will use our Hi-C data to scaffold the contigs. We will follow [this Omni-C protocol](https://omni-c.readthedocs.io/en/latest/index.html) for mapping and use yahs for scaffolding. 

```
samtools --version
samtools 1.15.1
Using htslib 1.16

bwa
Version: 0.7.17-r1188

pairtools --version
pairtools, version 1.0.2

yahs --version
1.2a.2
```

```
samtools faidx Kronos.draft.fa
bwa index Kronos.draft.fa

#align Hi-C read pairs
bwa mem -o aligned.sam -5SP -T0 -t52 /Kronos.draft.fa <(zcat 0.HiC/KVK-*/*R1*trimmed.fq.gz) <(zcat 0.HiC/KVK-*/*R2*trimmed.fq.gz)

#process the alignments
samtools view -@56 -h aligned.sam  \
pairtools parse --min-mapq 30 --walks-policy 5unique \
--max-inter-align-gap 30 --nproc-in 56 --nproc-out 56 --chroms-path Kronos.draft.fa | \
pairtools sort --tmpdir=./tmp --nproc 56 | pairtools dedup --nproc-in 56 \
--nproc-out 56 --mark-dups --output-stats stats.txt | pairtools split --nproc-in 56 \
--nproc-out 56 --output-pairs mapped.pairs --output-sam - |samtools view -bS -@56 | \
samtools sort -@56 -o mapped.PT.bam ; samtools index mapped.PT.bam

yahs -o YaHS -e GATC,GANTC,CTNAG,TTAA Kronos.draft.fa mapped.PT.bam
```

We can quickly check the length of the scaffolds. We expect 14 largest scaffolds (7 chromosomes x AB). 
```
python -c "from Bio import SeqIO; print('\n'.join([f'{record.id} {round(len(record.seq)/1000000, 3)} Mb' for record in SeqIO.parse('YaHS_scaffolds_final.fa', 'fasta')]))" | sort -r -nk 2 | head -n 20

scaffold_1 864.152 Mb
scaffold_2 828.542 Mb
scaffold_3 795.82 Mb
scaffold_4 767.866 Mb
scaffold_5 766.027 Mb
scaffold_6 759.128 Mb
scaffold_7 753.477 Mb
scaffold_8 733.6 Mb
scaffold_9 731.153 Mb
scaffold_10 720.28 Mb
scaffold_11 708.843 Mb
scaffold_12 699.697 Mb
scaffold_13 624.303 Mb
scaffold_14 600.444 Mb
scaffold_15 3.636 Mb
scaffold_16 3.606 Mb
scaffold_17 2.687 Mb
scaffold_18 2.63 Mb
scaffold_19 2.481 Mb
scaffold_20 2.432 Mb
```
It looks like the largest 14 scaffolds may be the chromosomes!

Let's generate a contact map and visualize through JuiceBox.
```
/global/scratch/users/skyungyong/Software/yahs/juicer pre -a -o out_JBAT YaHS.bin YaHS_scaffolds_final.agp Kronos.draft.fa.fai >out_JBAT.log 2>&1
java -jar juicer_tools.1.9.9_jcuda.0.8.jar pre out_JBAT.txt out_JBAT.hic <(cat out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')
```
The two outputs, out_JBAT.hic and out_JBAT.assembly, can be loaded into [Juicebox](https://github.com/aidenlab/Juicebox/wiki/Download). Set the scale as below:
```
grep 'scale factor' out_JBAT.log
[I::main_pre] scale factor: 8
```

For visualization, we will use the following setups:
```
Show: Log(Observed+1)
Normalization: Balanced
Color Range: 0-13-19 (for the entire assembly)
```

```
juicer pre YaHS.bin YaHS_scaffolds_final.agp Kronos.draft.fa.fai | sort -k2,2d -k6,6d -T ./ --parallel=30 -S256G | awk 'NF' > alignments_sorted.txt
python -c "from Bio import SeqIO; print('\n'.join([f'{record.id}\t{len(record.seq)}' for record in SeqIO.parse('YaHS_scaffolds_final.fa', 'fasta')]))"  > scaf.length
java -jar juicer_tools.1.9.9_jcuda.0.8.jar pre alignments_sorted.txt out.hic scaf.length
```
We did not observe any abnormal features from the contact map and are happy to move on!

## Scaffolding renaming

We will remove chloroplast and mitocondiral genomes into separate files and reassign the scaffold names. For the plasmids, we will download these two accessions from the NCBI. The wheat reference genome can be downloaded from [EnsemblPlants](https://plants.ensembl.org/Triticum_aestivum/Info/Index).
```
Triticum aestivum chloroplast, complete genome: NC_002762.1
Triticum aestivum cultivar Chinese Yumai mitochondrion, complete genome: NC_036024.1
Triticum_aestivum.IWGSC.dna.toplevel.fa
```

As the scaffolds are too large, minimap won't be able to run. Let's break the assemblies first. 
```
python break_fa.py YaHS_scaffolds_final.fa
python break_fa.py Triticum_aestivum.IWGSC.dna.toplevel.fa
```

We can then run minimap v2.24-r1122. This step took > 72 hours. It is recommended to split the query and submit multiple jobs.
```
cat NC_002762.1.fasta NC_036024.1.fasta > Triticum_aestivum.plasmids.fa
minimap2 -x asm5 -t 52 ../Triticum_aestivum.IWGSC.dna.toplevel.broken.fa YaHS_scaffolds_final.broken.fa > minimap.ref.paf
sort -k1,1 -k3,3n minimap.ref.paf > minimap.ref.sorted.paf
minimap2 -x asm5 -t 52 ../Triticum_aestivum.plasmids.fa YaHS_scaffolds_final.fa > minimap.plasmid.paf
sort -k1,1 -k3,3n minimap.plasmid.paf > minimap.plasmid.sorted.paf
```

Run the following script to reassign the scaffold names and separate plasmid DNAs. This will compare 14 Kronos scaffolds and 21 Wheat chromosomes and transfer the chromosome IDs. 
```
python process_scaffolds.py minimap.plasmid.sorted.paf minimap.ref.sorted.paf YaHS_scaffolds_final.fa Kronos.collapsed
```

| Chromosomes  | 1A | 1B | 2A | 2B | 3A | 3B | 4A | 4B | 5A | 5B | 6A | 6B | 7A | 7B | Un | 
|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|
| scaffold ID  | 14  | 11 | 3 | 2 | 6 | 1 | 4 | 12 | 10 | 9 | 13 | 8 | 7 | 5 | - | 
| size  | 600443981| 708842986| 795820389| 828541533| 759128228| 864152387| 767865717| 699696956| 720280059| 731153026| 624303373| 733599645| 753476766| 766026795| 211250944 | 
| N's | 4400| 16000| 7400| 18400| 3800| 20400| 13400| 18600| 4400| 21400| 5200| 20400| 5200| 20000| 731600 | 
| unambiguous base pairs | 600439581| 708826986| 795812989| 828523133| 759124428| 864131987| 767852317| 699678356| 720275659| 731131626| 624298173| 733579245| 753471566| 766006795| 210519344 |


51,529,289 base pairs where separated as chloroplast, and 8,923,050 as mitochondria. 


minimap2 -x asm5 -t 20 ../../Triticum_aestivum.plasmids.fa Kronos.draft.fa > min
imap.plasmid.paf
sort -k1,1 -k3,3n minimap.plasmid.paf > minimap.plasmid.sorted.paf
python ../process_scaffolds.py ../minimap.plasmid.sorted.paf None Kronos.draft.fa Kronos.contigs
quast -t 20 --fast Kronos.contigs.genomic.fa


## Repeat masking

We will use [HiTE](https://github.com/CSU-KangHu/HiTE) for repeat masking. This step took slightly more than two weeks. 
```
export _CACHEDIR=/global/scratch/users/skyungyong/temp
singularity pull HiTE.sif docker://kanghu/hite:3.0.0 nb
WD=$(pwd)
singularity run -B ${host_path}:${container_path} --pwd /HiTE  HiTE.sif python main.py --genome $WD/Kronos.collapsed.chromosomes.fa --thread 56 --outdir $WD/Kronos_output --recover 1 --plant 1 --classified 1 --domain 1
```

Then, we can mask the genome with the constructed library. To save time, we seperated each chromosome and masked it individually across multiple nodes with ReapeatMasker v4.1.5. 
```
RepeatMasker -xsmall -e ncbi -pa 56 -q -no_is -norna -nolow -div 40 -gff -lib /confident_TE.cons.fa.classified -dir $prefix\/ -cutoff 225 $prefix\/$prefix\.fa
```
For the unplaced scaffold (Un), RepeatMasker somehow got stuck in the processrepeat script. It turned out the software gets infinetly stuck in a while loop (around line 404). We didn not spend a lot of time troubleshooting this. Rather, we removed the problematic sequence (N_8755) from confident_TE.cons.fa.classified and re-run RepeatMasker. 

Then all masked fasta files are merged into one
```
cat */*.fa.masked > Kronos.collapsed.chromosomes.masked.fa
cat */*.fa.out > Kronos.collapsed.chromosomes.masked.rm.out
cat */*.out.gff > Kronos.collapsed.chromosomes.masked.out.gff
bedtools maskfasta  -fi Kronos.collapsed.chromosomes.masked.fa -bed Kronos.collapsed.chromosomes.masked.out.gff -fo Kronos.collapsed.chromosomes.hard-masked.fa
```


### Telomere search


According to [Telobase](http://cfb.ceitec.muni.cz/telobase/), the telomere sequences for *Triticum* is TTTAGGG. We can double check before we use this sequence.
```
tidk search ../../2.Scaffold/Kronos.collapsed.chromosomes.fa -s TTTAGGG -o tidk.search -d search
tidk plot --tsv search/tidk.search_telomeric_repeat_windows.tsv
```
It looks like 8 scaffolds may have telomeres at the both ends. 6 may have telomeres only in one end. Some telomeres might have not been scaffolded properly due to complexity. 


# Soft-mask the genome
singularity run -B ${host_path}:${container_path} --pwd /HiTE HiTE.sif RepeatMasker -xsmall -lib Haplotype-1/confident_TE.cons.fa.classified -dir Haplotype-1 -pa 20 SH1353.haplotype-1.fa

## Genome annotation: RNA-seq

We donwloaded the paired-end RNA-seq data from the NCBI. The list can be found in [SRA.list](https://github.com/s-kyungyong/Kronos/blob/main/RNAseq/SRA.list). Let's first remove adapters and low-quality reads from the libraries. Trim_galore and cutadapt versions were v0.6.6 and v3.7. 
```
ls *.fastq | cut -d "_" -f 1 | sort -u | while read accession; do trim_galore -a -j 8 --paired $accession\_1.fastq $accession\_2.fastq ; done
```
This generated about 1.6 Tb of fastq files in total

### De novo assembly

De novo assembly is done with Trinity v2.15.1. We initially tried running Trinity on the 1.6 Tb of paired-end fastq files all at once. After two weeks, Trinity was still stuck at the insilico normalization step with about 35-45% progress. We, therefore, had to take some other ways around.

Let's normalize each paired-end library first. We submitted ~50 jobs across different nodes to speed this step. 
```
read=/global/scratch/users/skyungyong/Kronos/4.RNAseq
prefix=$1
left=$read/$prefix\_1_val_1.fq
right=$read/$prefix\_2_val_2.fq

singularity run -B $PWD /global/scratch/users/skyungyong/Software/trinity.sif Trinity --verbose --max_memory 90G --just_normalize_reads --seqType fq --CPU 40 --left $left --right $right --output trinity_$prefix
```

cat trinity_S*/insilico_read_normalization/*1_val_1*.fq > right.norm.all.fq
cat trinity_S*/insilico_read_normalization/*2_val_2*.fq > left.norm.all.fq

Then, create a file that describes the samples and run Trinity.
```
ls -d trinity_* | while read folder; do prefix=$(echo $folder | cut -d "_" -f 2); left=$(ls $(pwd)/$folder\/insilico_read_normalization/*_1_val_1*.fq); right=$(ls $(pwd)\/$folder\/insilico_read_normalization/*_2_val_2*.fq); echo $prefix $prefix $left $right; done > sample.list

singularity run -B $PWD /global/scratch/users/skyungyong/Software/trinity.sif Trinity  --verbose --seqType fq --max_memory 1500G --CPU 56 --samples_file $PWD/sample.list
```
After re-normalization, Trinity produced ~56 Gb of the paired-end reads. It took a few days to run the software, the resulting transcripts were ~1 Gb. 

Let's process the transcripts with TransDecoder v5.7.1.

```
TransDecoder.LongOrfs -t trinity_out_dir.Trinity.fasta

wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmsearch --cpu 56 -E 1e-10 --domtblout pfam.domtblout Pfam-A.hmm longest_orfs.pep

# diamond version is 2.0.15
diamond makedb --in uniprotkb_taxonomy_id_38820_2023_12_08.fasta --db uniprot_38820_diamond # see below to check where these sequences came from
diamond blastp --threads 56 --evalue 1e-5 --db uniprot_38820_diamond.dmnd --max-target-seqs 1 --out diamond.outfmt6 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --query longest_orfs.pep

TransDecoder.Predict --retain_pfam_hits pfam.domtblout --retain_blastp_hits diamond.outfmt6 -t trinity_out_dir.Trinity.fasta
```

Not all transcripts would be in good quality, and we would need to process them later on. 
    

### Mapping

The RNA-seq data can be also mapped to the genome and processed. We are using hisat v2.2.1 to map the paired end libraries. The reads will be aligned as below. 
```
# index the genome 
hisat2-build -p 20 ../2.Scaffold/Kronos.collapsed.chromosomes.fa Kronos

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
```

We will also merge the bam files into one and use it later for the genome annotation with BRAKER. 
```
samtools merge -@ 56 -h SRX10965366.mapped.bam -o all.merged.bam *.mapped.bam
samtools sort -@ 56 all.merged.bam > all.merged.sorted.bam
```

We can then use stringtie v2.2.1 to process this. 
```
stringtie -o stringtie.gtf -p 56 --conservative all.merged.sorted.bam
```

### Miniprot

/global/scratch/users/skyungyong/Software/miniprot/miniprot -t 56 --gff --outc=0.95 -N 0 ../3.Repeat/Kronos_output_latest/Kronos.collapsed.chromosomes.fa ../5.Annotations/Braker/uniprotkb_taxonomy_id_38820_2023_12_08.fasta

### Braker

To run BRAKER, we need three inputs 1) a soft-masked genome (Kronos.collapsed.chromosomes.masked.fa), 2) mapped RNA-seq data (all.merged.sorted.bam), and 3) protein datasets. We downloaded 2,850,097 sequences from UniProt by searching for the following: 

```
(taxonomy_id:38820)
```

Then, Run BRAKER.
```
singularity build braker3.sif docker://teambraker/braker3:latest
singularity exec -B $PWD:$PWD braker3.sif cp -r /usr/share/augustus/config .

singularity exec braker3.sif braker.pl --verbosity=3 \
    --genome=Kronos.collapsed.chromosomes.masked.fa \
    --bam=all.merged.sorted.bam \
    --prot_seq=protkb_taxonomy_id_38820_2023_12_08.fasta \
    --species=Kronos_collapsed --threads 48 --gff3 \
    --nocleanup \
    --workingdir=./braker \
    --AUGUSTUS_CONFIG_PATH=./config
```
Then, UTR was added as below.

### Ginger

conda create ginger && conda activate ginger
conda install mamba python=3.9
mamba install salmon gffread nextflow=21.04 pasa bwa bowtie2 transdecoder hisat2 samtools stringtie velvet oases trinity gmap cd-hit seqkit spaln=3.0.0 augustus snap 
perl-carp perl-pathtools perl-data-dumper perl-db_file perl-findbin perl-uri perl-exporter perl-dbi perl-parallel-forkmanager perl-getopt-long matplotlib

VELVET
make 'LONGSEQUENCES=1' 'MAXKMERLENGTH=57' 'CATEGORIES=2'

git clone https://github.com/i10labtitech/GINGER.git && cd GINGER
mkdir util/mapping util/evaluation
no velvet -> 16,777,216 Giga bytes
make


wge https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/triticum_aestivum/pep/Triticum_aestivum.IWGSC.pep.all.fa.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/triticum_turgidum/pep/Triticum_turgidum.Svevo.v1.pep.all.fa.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/triticum_dicoccoides/pep/Triticum_dicoccoides.WEWSeq_v.1.0.pep.all.fa.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/triticum_spelta/pep/Triticum_spelta.PGSBv2.0.pep.all.fa.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/triticum_urartu/pep/Triticum_urartu.IGDB.pep.all.fa.gz
gunzip * 


samtools view -h -b -f 3 all.merged.sorted.bam > 

### Maker

We will first train SNAP and AUGUSTUS. AUGUSTUS is retrained to get different parameters, hoping that it will capture gene structures missed in BRAKER's run. Let's first get reliable gene models. We will consider gene models reliable if BRAKER's and Trinity's models agree with 100% coverage and identity. 

diamond makedb --in ../Trinity/trinity_out_dir.Trinity.fasta.transdecoder.pep --db trinity.transdecoder
diamond blastp --threads 56 --evalue 1e-10 --db trinity.transdecoder --max-target-seqs 1 --out braker.against.trinity.transdecoder.diamond.out --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --query ../Braker/braker_rerun/braker.aa

less braker.against.trinity.transdecoder.diamond.out | awk '$13 == $14 && $7 == 1 && $9 == 1 && $8 == $10 && $3 > 98 && $4 > 200 {print}' > train.initial.list
less train.initial.list | wc -l
8662

We will randomly split these into two sets and use to train SNAP and AUGUSTUS, respectively. Some of these transcripts come from the same gene, so we will filter those out. 
python select_genemodels.py train.initial.list ../Braker/braker_rerun/braker.gff3
augustus.train.gff has 3734 genes
snap.train.gff has 3734 genes



#SNAP
#Download gff3_to_zff.pl from here: https://biowize.wordpress.com/2012/06/01/training-the-snap-ab-initio-gene-predictor/
perl gff3_to_zff.pl < snap.train.gff3 > genome.ann
#The order of the scaffolds needs to be the same with the genome.ann file
grep '>' genome.ann | cut -d ">" -f 2 |  while read line; do awk -v seq=$line -v RS=">" '$1 == seq {print R
S $0; exit}' ../../../3.Repeat/Kronos_output_latest/Kronos.collapsed.chromosomes.fa; done > genome.dna

fathom -gene-stats genome.ann genome.dna > gene-stats.log
fathom -validate genome.ann genome.dna > gene.validate
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.*
fathom -validate export.ann export.dna
forge export.ann export.dna

hmm-assembler.pl Kronos . > Sohab.hmm



gff2gbSmallDNA.pl ../augustus.train.gff ../../Braker/Kronos.collapsed.chromosomes.masked.fa 2000 genes
randomSplit.pl genes.gb 200
new_species.pl --species=Kronos_re
etraining --species=Kronos_re genes.gb.train
augustus --species=Kronos_re genes.gb.test | tee first-test.out

*******      Evaluation of gene prediction     *******

---------------------------------------------\
                 | sensitivity | specificity |
---------------------------------------------|
nucleotide level |       0.893 |       0.757 |
---------------------------------------------/

----------------------------------------------------------------------------------------------------------\
           |  #pred |  #anno |      |    FP = false pos. |    FN = false neg. |             |             |
           | total/ | total/ |   TP |--------------------|--------------------| sensitivity | specificity |
           | unique | unique |      | part | ovlp | wrng | part | ovlp | wrng |             |             |
----------------------------------------------------------------------------------------------------------|
           |        |        |      |                414 |                310 |             |             |
exon level |   1458 |   1354 | 1044 | ------------------ | ------------------ |       0.771 |       0.716 |
           |   1458 |   1354 |      |  127 |    8 |  279 |  127 |    7 |  176 |             |             |
----------------------------------------------------------------------------------------------------------/

----------------------------------------------------------------------------\
transcript | #pred | #anno |   TP |   FP |   FN | sensitivity | specificity |
----------------------------------------------------------------------------|
gene level |   315 |   200 |   84 |  231 |  116 |        0.42 |       0.267 |
----------------------------------------------------------------------------/



optimize_augustus.pl --species=Sohab --cpus=8 --UTR=off genes.gb.train
augustus --species=Sohab genes.gb.test | tee second-test.out

MAKER will be run with trained SNAP and AUGUSTUS as well as 



wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/triticum_aestivum/pep/Triticum_aestivum.IWGSC.pep.all.fa.gz

### Evience modeler 



### NLR annotation

