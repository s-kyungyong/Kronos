# Genome assembly

## Data availability 
Sequencing data were deposted in the NCBI under the BioProject assession, PRJNA1213727. The following runs include:  
• `HiFi reads`: SRR32063042  
• `PacBio sequencing data (bam)`: SRR32063043  
• `Hi-C data`: SRR32063044  

The Kronos reference genome can be assessed through the NCBI and Zenodo.   
• `v1.0 genome`: https://zenodo.org/records/10215402  
• `v1.1 genome [Final]`: https://zenodo.org/records/11106422  
• `v1.1 genome [no Un]`: JBLYQA000000000  
• `Haplotype resolved assembly (contigs)`:   

Note that the genome acessible through the NCBI does not have **Un** sequences. Furthermore, the following genome regions were hard-masked due to some similarity to mitochondria.
```
2A	795820389	239269473..239297248
2B	828541533	290778354..290821184
2B	828541533	290821385..290871421
3B	864152387	598683503..598723947
3B	864152387	598724148..598754566
4A	767865717	738117070..738151875
5B	731153026	206904209..206929881
6B	733599645	394727123..394769363
6B	733599645	394769564..394869826
7B	766026795	742909075..742941679
```

## Software version

```
GenomeScope v2.0
Jellyfish v2.2.10
trim_galore v0.6.6
cutadapt v3.7
bam2fastq v3.0.0
hifiasm v0.19.5-r587
fastp v0.23.2
bwa v0.7.17-r1188
pairtools v1.0.2
samtools v1.15.1
yahs v1.2a.2
juicer v1.9.9
Juicebox v2.20.00
minimap v2.24-r1122
BLAST v2.15.0
SyRI v1.7.0
plotsr v1.1.0
```
---
### 1. Quality Control
This step removes adapters from sequencing data. 

**📥 Inputs**  
• `*.hifi_reads.default.bam`: HiFi data in bam format
• `*.fastq.gz`: Hi-C data in fastq format  

**📥 Outputs**  
• `Kronos.HiFi.filt.fastq.gz`: HiFi data in fastq format
• `*.filtered.fastq.gz`: trimmed filtered Hi-C data in fastq format  

⚙️ **Convert to FastQ**  
About 50X HiFi reads were obtained from Revio. 
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
```
reads=$(ls *.default.bam)
bam2fastq -o Kronos.HiFi -j 52 $reads
```

⚙️ **Remove Adapters (optional)**  
Here, all filtered reads were combined.
```
HiFiAdapterFilt/hifiadapterfilt.sh -p Kronos -t 54

#check the statistics
cat Kronos.HiFi.stats
Started on Wed Jun 14 21:12:00 PDT 2023
For the Kronos.HiFi dataset:
Removing reads containing adapters a minimum of 44 bp in length and 97% match.

Number of ccs reads: 39390905
Number of adapter contaminated ccs reads: 13 (3.30025e-05% of total)
Number of ccs reads retained: 39390892 (100% of total)
```

⚙️ **Trim Hi-C datasets**  
```
ls -lha  KVK-KRONOS-*/*.fastq.gz
-rwxr-xr-x 1 skyungyong ucb 49G Jul 26 11:06 KVK-KRONOS-874809/kvk-kronos-874809_S3HiC_R1.fastq.gz
-rw-r--r-- 1 skyungyong ucb 50G Jul 26 10:46 KVK-KRONOS-874809/kvk-kronos-874809_S3HiC_R2.fastq.gz
-rwxr-xr-x 1 skyungyong ucb 55G Jul 26 11:21 KVK-KRONOS-HIC2-1066975/1066975_S1_R1_001.fastq.gz
-rwxr-xr-x 1 skyungyong ucb 56G Jul 26 11:37 KVK-KRONOS-HIC2-1066975/1066975_S1_R2_001.fastq.gz
-rwxr-xr-x 1 skyungyong ucb 60G Jul 26 12:27 KVK-KRONOS-HIC2-1067017/kvk-kronos-hic2-1067017_S3HiC_R1.fastq.gz
-rwxr-xr-x 1 skyungyong ucb 62G Jul 26 12:01 KVK-KRONOS-HIC2-1067017/kvk-kronos-hic2-1067017_S3HiC_R2.fastq.gz
```
```
#out can be formatted as *.filtered.fastq.gz
fastp --thread 54 -g --detect_adapter_for_pe -q 20 -l 50  \
      --in1 $read1 --in2 $read2 --out1 $out1 --out2 $out2 \
      -h bssh1.html &> bssh1.log
fastqc -t 54 $out1 $out2
```

---
### 2. Genome Assessment
GenomeScope was used to assess genomic characteristics of Kronos.

**📥 Inputs**  
• `Kronos.HiFi.filt.fastq.gz`: HiFi data in bam format


⚙️ **Run GenomeScope**  
```
#count k-mer
jellyfish count -C -m 21 -s 50000000000 -t 20 Kronos.HiFi.filt.fastq -o kmer_counts.jf
jellyfish histo -h 5000000 -t 20 kmer_counts.jf > reads.histo

genomescope2 -p 4 -i reads.histo -o genomescope --verbose  
```
⚙️ **Comparison to Svevo**    
The same analysis was performed for Svevo for comparison.
```
#download sequencing data
while read -r accession; do 
    sratoolkit.3.1.1-centos_linux64/bin/prefetch ${accession}
    sratoolkit.3.1.1-centos_linux64/bin/fasterq-dump -O . -e ${Numthreads} ${accession}
done < Svevo_SRA.list

#trim
ls *.fastq | cut -d "_" -f 1 | sort -u |  while read accession; do
    trim_galore --illumina -j 8 --paired $accession\_1.fastq $accession\_2.fastq
done

#genomescope
jellyfish count -C -m 21 -s 50000000000 -t 56 *.fq -o svevo.kmer_counts.jf
jellyfish histo -h 5000000 -t 56 svevo.kmer_counts.jf > svevo.reads.histo
genomescope2 -p 4 -i svevo.reads.histo -o svevo.genomescope --verbose
```

⚙️ **Statistics (Kronos vs Svevo)**  
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

⚙️ **Statistics (Kronos vs Chinese Spring)**  
Note: The statistics for bread weak comes Fig. S21 from [this paper](https://www.nature.com/articles/s41467-020-14998-3). Here, the genome size was esimated as haplotype size x ploidy. Although this isn't technically correct, it seems to give the right estimate. See the discussion [here](https://github.com/schatzlab/genomescope/issues/107).  
|    | Kronos | Svevo| Bread wheat |
|----|---------|-----------| -----------|
| ploidy (p) | 4 | 4 | 6 | 
| Haplotype size (Gb) | 2.612  | 2.754 | 2.354 |
| Genome size (Gb)  | 10.45 | 11.02 | 14.1 |
| unique (%) | 9.52 | 9.3 | 8 |
| heterozygosity (%) | 9.5 | 9.4 |  10.1 |
| err (%) | 0.0707 | 0.191 | 0.506 |
| dup | 0.654 | 2.67 | 0.836 | 

---

### 3. Genome Assembly
Both haplotype-collapsed and haplotype-resolved assemblies were generated. 

**📥 Inputs**  
• `Kronos.HiFi.filt.fastq.gz`: HiFi data in fastq format
• `*.fastq.gz`: Hi-C data in fastq format (for haplotype-resolved assembly)

**📥 Outputs**  
• `*Kronos.draft.fa`: Haplotype-collapsed assembly
• `Kronos.haplotype_resolved.HAP1.contigs.v1.0.fa`: Haplotype-resolved assembly (Haplotype 1)
• `Kronos.haplotype_resolved.HAP2.contigs.v1.0.fa`: Haplotype-resolved assembly (Haplotype 2)

⚙️ **Create Haplotype-collapsed Assembly**  
As Kronos' heterozygosity is low, and we aimed to generate collapsed haplotypes (AB). This took 61 hours and 615 Gb of a peak memory.
```
hifiasm -l0 -t 54 -o l0 Kronos.HiFi.filt.fastq.gz

#convert gfa to fasta files
awk '/^S/{print ">"$2"\n"$3}' l0.bp.p_ctg.gfa | fold > Kronos.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' l0.bp.a_ctg.gfa | fold > Kronos.a_ctg.fa

#concatnate primary and associate contigs
cat Kronos.p_ctg.fa Kronos.a_ctg.fa > Kronos.draft.fa
```

⚙️ **Create Haplotype-resolved Assembly**  
```
hifiasm -l0 -t 54 -o l0-hic --s-base -1 -h1 $hic_pair1 --h2 $hic_pair2 Kronos.HiFi.filt.fastq.gz

#convert gfa to fasta files
awk '/^S/{print ">"$2"\n"$3}' l0.bp.p_ctg.gfa | fold > Kronos.haplotype_resolved.HAP1.contigs.v1.0.fa
awk '/^S/{print ">"$2"\n"$3}' l0.bp.a_ctg.gfa | fold > Kronos.haplotype_resolved.HAP2.contigs.v1.0.fa
```

### 4. Scaffolding 
Haplotype-collapsed assembly was scaffolded with Hi-C datasets. 

**📥 Inputs**  
• `*Kronos.draft.fa`: Haplotype-collapsed assembly
• `*.fastq.gz`: Hi-C data in fastq format (for haplotype-resolved assembly)

**📥 Outputs**  
• `*Kronos.draft.fa`: Haplotype-collapsed assembly
• `Kronos.haplotype_resolved.HAP1.contigs.v1.0.fa`: Haplotype-resolved assembly (Haplotype 1)
• `Kronos.haplotype_resolved.HAP2.contigs.v1.0.fa`: Haplotype-resolved assembly (Haplotype 2)

⚙️ **Indexing and Alignment** 
```
samtools faidx Kronos.draft.fa
bwa index Kronos.draft.fa

#align Hi-C read pairs
bwa mem -o aligned.sam -5SP -T0 -t52 Kronos.draft.fa <(zcat 0.HiC/KVK-*/*R1*trimmed.fq.gz) <(zcat 0.HiC/KVK-*/*R2*trimmed.fq.gz)
```
⚙️ **Omni-C pipeline** 
```
#instruction: https://omni-c.readthedocs.io/en/latest/index.html
samtools view -@56 -h aligned.sam  \
pairtools parse --min-mapq 30 --walks-policy 5unique \
--max-inter-align-gap 30 --nproc-in 56 --nproc-out 56 --chroms-path Kronos.draft.fa | \
pairtools sort --tmpdir=./tmp --nproc 56 | pairtools dedup --nproc-in 56 \
--nproc-out 56 --mark-dups --output-stats stats.txt | pairtools split --nproc-in 56 \
--nproc-out 56 --output-pairs mapped.pairs --output-sam - |samtools view -bS -@56 | \
samtools sort -@56 -o mapped.PT.bam 
samtools index mapped.PT.bam
```
⚙️ **Scaffolding** 
```
yahs -o YaHS -e GATC,GANTC,CTNAG,TTAA Kronos.draft.fa mapped.PT.bam
```

Let's check the scaffold lengths. We expect 14 largest scaffolds to correspond to 7 chromosomes from A and B subgenomes. 
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
It looks like the largest 14 scaffolds are the chromosomes! The other anchored sequences are all smaller than 4 Mb. 

Let's generate a contact map and visualize through JuiceBox. The two outputs, out_JBAT.hic and out_JBAT.assembly, can be loaded into [Juicebox](https://github.com/aidenlab/Juicebox/wiki/Download). Set the scale as below. We did not observe any significant abnormality, so we did not manually change anything.
```
./yahs/juicer pre -a -o out_JBAT YaHS.bin YaHS_scaffolds_final.agp Kronos.draft.fa.fai > out_JBAT.log 2>&1
java -jar juicer_tools.1.9.9_jcuda.0.8.jar pre out_JBAT.txt out_JBAT.hic <(cat out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')

grep 'scale factor' out_JBAT.log
[I::main_pre] scale factor: 8
```

### 3C. Final Processing

Remove chloroplast and mitocondiral genomes into separate files and reassign the scaffold names. For the plasmids, we will download these two accessions from the NCBI below. The wheat reference genome can be downloaded from [EnsemblPlants](https://plants.ensembl.org/Triticum_aestivum/Info/Index).
```
Triticum aestivum chloroplast, complete genome: NC_002762.1
Triticum aestivum cultivar Chinese Yumai mitochondrion, complete genome: NC_036024.1
Triticum_aestivum.IWGSC.dna.toplevel.fa
```

As the scaffolds are too large, minimap runs slow. Let's break the assemblies first. 
```
python break_fa.py YaHS_scaffolds_final.fa
python break_fa.py Triticum_aestivum.IWGSC.dna.toplevel.fa
```

We can then run minimap v2.24-r1122. This step took > 72 hours. It is recommended to split the query and submit multiple jobs.
```
cat NC_002762.1.fasta NC_036024.1.fasta > Triticum_aestivum.plasmids.fa
minimap2 -x asm5 -t 52 Triticum_aestivum.IWGSC.dna.toplevel.broken.fa YaHS_scaffolds_final.broken.fa > minimap.ref.paf
sort -k1,1 -k3,3n minimap.ref.paf > minimap.ref.sorted.paf
minimap2 -x asm5 -t 52 ../Triticum_aestivum.plasmids.fa YaHS_scaffolds_final.fa > minimap.plasmid.paf
sort -k1,1 -k3,3n minimap.plasmid.paf > minimap.plasmid.sorted.paf
```

Run the following script to reassign the scaffold names and separate plasmid DNAs. This will compare 14 Kronos scaffolds and 21 Wheat chromosomes and transfer the chromosome IDs. 
```
python process_scaffolds.py minimap.plasmid.sorted.paf minimap.ref.sorted.paf YaHS_scaffolds_final.fa Kronos.collapsed
```

This step creates the Kronos reference v1.0, with the following statistics. 

| Chromosomes  | 1A | 1B | 2A | 2B | 3A | 3B | 4A | 4B | 5A | 5B | 6A | 6B | 7A | 7B | Un | 
|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|
| scaffold ID  | 14  | 11 | 3 | 2 | 6 | 1 | 4 | 12 | 10 | 9 | 13 | 8 | 7 | 5 | - | 
| size  | 600443981| 708842986| 795820389| 828541533| 759128228| 864152387| 767865717| 699696956| 720280059| 731153026| 624303373| 733599645| 753476766| 766026795| 211250944 | 
| N's | 4400| 16000| 7400| 18400| 3800| 20400| 13400| 18600| 4400| 21400| 5200| 20400| 5200| 20000| 731600 | 
| unambiguous base pairs | 600439581| 708826986| 795812989| 828523133| 759124428| 864131987| 767852317| 699678356| 720275659| 731131626| 624298173| 733579245| 753471566| 766006795| 210519344 |

That's our genome v1.0! In the version 1.1, chromosomes 1B, 2A, 2B, 3A, 3B, 5A, 6A and 6B are flipped to make their orientations consistent with the Chinese Spring genome. 


##4. Synteny analyses

The Kronos reference genome was compared to Svevo and Chinese Spring (IWGSC v1.0), which were downloaded from Plant Ensembl. For global synteny, the similarity search was performed using minimap.
```
#run minimap
minimap2 --eqx -c -f 0.05 -K4g -t 30 -x asm5 Kronos.collapsed.chromosomes.masked.v1.1.fa Triticum_turgidum.Svevo.v1.dna.toplevel.fa -o Kronos_vs_Svevo.eqx_asm5.paf
minimap2 --eqx -c -f 0.05 -K4g -t 30 -x asm5 Kronos.collapsed.chromosomes.masked.v1.1.fa Triticum_aestivum.IWGSC.dna.toplevel.fa -o Kronos_vs_CS.eqx_asm5.paf

#filter
awk '$1 != "Un" && $1 !~ /[0-9]+D/' Kronos_vs_Svevo.eqx_asm5.paf > Kronos_vs_Svevo.eqx_asm5.filtered.paf
awk '$1 != "Un" && $1 !~ /[0-9]+D/' Kronos_vs_CS.eqx_asm5.paf > Kronos_vs_CS.eqx_asm5.filtered.paf

#plot
python synteny_minimap_synteny.py Kronos_vs_Svevo.eqx_asm5.filtered.paf 1A 1A
```

Structural variants were detecte using syri.
```
#analyze genomic variants
syri -r Kronos.collapsed.chromosomes.masked.v1.1.fa -q Triticum_turgidum.Svevo.v1.dna.toplevel.fa -c Kronos_vs_Svevo.eqx_asm5.filtered.paf -F P -k --prefix Kronos_vs_Svevo.
syri -r Kronos.collapsed.chromosomes.masked.v1.1.fa -q Triticum_aestivum.IWGSC.dna.toplevel.fa -c Kronos_vs_CS.eqx_asm5.filtered.paf -F P -k --prefix Kronos_vs_CS.

#visualize as plots
#genomeList files simply have location of each genome and its name as two columns
plotsr --sr Kronos_vs_Svevo.syri.out --genomes genomeList.1.txt -o Kronos_vs_Svevo.pdf
plotsr --sr Kronos_vs_CS.syri.out --genomes genomeList.2.txt -o Kronos_vs_CS.pdf
```

For local synteny, BLAST v2.15.0 was used.
```
makeblastdb -in db.fa -out db -dbtype 'nucl'
blastn -query query.fa -db db -outfmt "6 std qlen slen" -out query_vs_hit.blast.out -evalue 1e-6
python synteny_blast_synteny.py query_vs_hit.blast.out 1A_start 1A_start --alignment_length 10000 --hstart 0 --hend 20000000
```

for chromosome in 1A 1B 2A 2B 3A 3B 4A 4B 5A 5B 6A 6B 7A 7B; do
  fasta=GCF_002162155.2_WEW_v2.1_genomic.renamed.fa
  length=$(awk -v chr=$chromosome '$1 == chr {print $2}' ${fasta}.fai)
  tail_start=$((length - 20000000 + 1))

  samtools faidx $fasta ${chromosome}:1-20000000 > Zavitan_${chromosome}.head.fa
  samtools faidx $fasta ${chromosome}:${tail_start}-${length} > Zavitan_${chromosome}.tail.fa
done

for fa in Svevo*fa; do makeblastdb -in ${fa} -out $(echo $fa | cut -d "." -f -2) -dbtype 'nucl'; done

for fa in Kronos*.fa; do 
  # Extract base name
  base=$(basename $fa .fa)
  region=$(echo $base | cut -d "_" -f 2)

  # Svevo comparison
  svevo_db=$(echo $fa | sed 's/Kronos/Svevo/')
  blastn -query ${fa} -db ${svevo_db} -num_threads 48 \
    -outfmt "6 std qlen slen" \
    -out Kronos_vs_Svevo.${region}.blastn.out -evalue 1e-4

  # Zavitan comparison
  zavitan_db=$(echo $fa | sed 's/Kronos/Zavitan/')
  blastn -query ${fa} -db ${zavitan_db} -num_threads 48 \
    -outfmt "6 std qlen slen" \
    -out Kronos_vs_Zavitan.${region}.blastn.out -evalue 1e-4
done


for fa in Kronos*.fa; do 
  base=$(basename $fa | sed "s/.fa//g")
  region=$(echo $base | cut -d "_" -f 2)
  minimap2 -x asm5 -t 40 -o Kronos_vs_Svevo.${region}.minimap.out ../0.Syri/Triticum_turgidum.Svevo.v1.dna.toplevel.fa ${fa}
  minimap2 -x asm5 -t 40 -o Kronos_vs_Zavitan.${region}.minimap.out ../0.Syri/GCF_002162155.2_WEW_v2.1_genomic.renamed.fa ${fa}
done

