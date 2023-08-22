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

Filter the reads with hifiadapterfilt. This step is not necessary. Often, the HiFi reads come out pretty clean. However, we occasionally had cases where not-so-clean HiFi reads resulted in adapters contained inside the contigs. We thus included this step for the quality control. 
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

Kronos is an allotetraploid species (AABB). In the field, Kronos rarely out-crosses. We therefore think that the heterozygosity should be extremely low or residual, and we aim to generate collapsed haplotypes (AB). An available Durum wheat genome ((Svevo)[https://www.nature.com/articles/s41588-019-0381-3]) is 10.45G in size. We also roughly estimate that each haplotype will be 5-6G in size, totalling up to 10-12G for the collapsed haplotypes. We will evaluate our assumptions with GenomeScope. 


```


```


## Genome assembly
