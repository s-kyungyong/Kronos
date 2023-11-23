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
singularity exec docker://huangnengcsu/compleasm:v${VERSION} compleasm run -l poales_odb10 -t 20 -o a_busco -a Kronos.a_ctg.fa
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
python process_scaffolds.py
```

| Chromosomes  | 1A | 1B | 2A | 2B | 3A | 3B | 4A | 4B | 5A | 5B | 6A | 6B | 7A | 7B | Un | 
|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|
| scaffold ID  | 14  | 11 | 3 | 2 | 6 | 1 | 4 | 12 | 10 | 9 | 13 | 8 | 7 | 5 | - | 
| size  | 600443981| 708842986| 795820389| 828541533| 759128228| 864152387| 767865717| 699696956| 720280059| 731153026| 624303373| 733599645| 753476766| 766026795| 211250944 | 
| N's | 4400| 16000| 7400| 18400| 3800| 20400| 13400| 18600| 4400| 21400| 5200| 20400| 5200| 20000| 731600 | 
| unambiguous base pairs | 600439581| 708826986| 795812989| 828523133| 759124428| 864131987| 767852317| 699678356| 720275659| 731131626| 624298173| 733579245| 753471566| 766006795| 210519344 |


51,529,289 base pairs where separated as chloroplast, and 8,923,050 as mitochondria. 



## Repeat masking

We will use [HiTE](https://github.com/CSU-KangHu/HiTE) for repeat masking. 

```
export SINGULARITY_CACHEDIR=/global/scratch/users/skyungyong/temp
singularity pull HiTE.sif docker://kanghu/hite:3.0.0
WD=$(pwd)
singularity run -B ${host_path}:${container_path} --pwd /HiTE  HiTE.sif python main.py --genome $WD/Kronos.collapsed.chromosomes.fa --thread 56 --outdir $WD/Kronos_output --recover 1 --annotate 1 --plant 1 --classified 1 --domain 1

# soft-mask genomes
singularity run -B ${host_path}:${container_path} --pwd /HiTE HiTE.sif RepeatMasker -e ncbi -pa 40 -q -no_is -norna -nolow -div 40 -gff -lib confident_TE.cons.fa.classified -cutoff 225 ${your_genome_path} && mv ${your_genome_path}.out ${HiTE_out_dir}/HiTE.out && mv ${your_genome_path}.tbl ${HiTE_out_dir}/HiTE.tbl && mv ${your_genome_path}.out.gff ${HiTE_out_dir}/HiTE.gff
```

According to [Telobase](http://cfb.ceitec.muni.cz/telobase/), the telomere sequences for *Triticum* is TTTAGGG. We can double check before we use this sequence.
```
tidk search ../../2.Scaffold/Kronos.collapsed.chromosomes.fa -s TTTAGGG -o tidk.search -d search
tidk plot --tsv search/tidk.search_telomeric_repeat_windows.tsv
```
It looks like 8 scaffolds may have telomeres at the both ends. 6 may have telomeres only in one end. Some telomeres might have not been scaffolded properly due to complexity. 


# Soft-mask the genome
singularity run -B ${host_path}:${container_path} --pwd /HiTE HiTE.sif RepeatMasker -xsmall -lib Haplotype-1/confident_TE.cons.fa.classified -dir Haplotype-1 -pa 20 SH1353.haplotype-1.fa

## RNA-seq

ls *.fastq | cut -d "_" -f 1 | sort -u | while read accession; do trim_galore -a -j 8 --paired $accession\_1.fastq $accession\_2.fastq ; done


We are using hisat v2.2.1 to map the paired end libraries. The reads will be aligned as below. 
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
samtools -@56 -q 30 -h -b
```
However, some of the sequencing data are big. For instance, SRX10965365, SRX10965366, and SRX10965367 are 460G, 340G and 510G in size, respectively. We will map individual or some combined paired-end libraries, and then merge them into a single alignment file later to speed up this process.



