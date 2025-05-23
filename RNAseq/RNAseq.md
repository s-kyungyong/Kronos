- RNAseq filtering
```
# Construct the fastp command
            cmd = [
                fastp_path,
                "--in1", in1,
                "--out1", out1,
                "--in2", in2,
                "--out2", out2,
                "-q", "20",
                "--detect_adapter_for_pe",
                "-w", "16"
            ]
```
- Parameters for STAR
1. Indexing
```
--runThreadN 56
--runMode genomeGenerate
--genomeDir /~
--genomeFastaFiles /~
--sjdbGTFfile /annotation.gtf or gff3
--sjdOverhang Decide later when we get the dataset and look at the read length
--sjdbGTFtagExonParentTranscript Parent
```
2. runnning:
```
--runThreadN 56
--genomeDir /path to genomeDir (where the indexed genomea is stored)
--readFilesIn
--outFilterMultimapNmax 5 (the output of multimappers (reads mapping to multiple loci) is controlled by this command. Default N=10. If a read maps to <= N loci, do not discard it)
--outFileNamePrefix /path to output dir prefix (follows the SRR* run names)
--outSAMprimaryFlag AllBestScore (change the default behavior to keeping every same-score alignments as the primary alignments for downstream analysis)
```
Multimappers: 
a. The number of loci Nmap a read maps to is given by NH:i:Nmap field. A value of 1 corresponds to unique mappers, while values >1 correspond to multi-mappers. The mapping quality MAPQ (column 5) is 255 for unique mapping reads, and int(-10*log10(1-1/Nmap)) for multi-mapping reads. (might use the i=1 to find unique mapper)
b. Might need to consider the definition of "unique" in the manual. The unique doesn't necessarily mean unique primary alignment over the whole read or the unique mapping alignment but not the multi-mapping alignment. 

3. Counting the alignments using featureCounts with feature of count multi-mapping reads and multi-overlapping reads
```
-primary (if specified, only primary alignments will be counted. Primary and secondary are identified using bit 0x100 in the Flag field of SAM/BAM files)
```
SAMtools:
256 for filtering out secondary alignment
4 for unmapped (if specify 4, it will filter out all unmapped alignments)

- Annotation gene rescue:
1. Script for appending missing gene modules

# Procedures:
- Creating reference database, so-called indexing, for different read length data (150bp/100bp):
```
mkdir 150bp && STAR --runThreadN 56 --runMode genomeGenerate --genomeDir 150bp --genomeFastaFiles /global/scratch/projects/vector_kvklab/KS-Kronos_remapping/Kronos.collapsed.chromosomes.masked.v1.1.fa --sjdbOverhang 149 --sjdbGTFfile ../../KS-Makido/Processing/65.Final_Final/Kronos.v2.0.gff3 --sjdbGTFtagExonParentTranscript Parent

mkdir 100bp && STAR --runThreadN 56 --runMode genomeGenerate --genomeDir 100bp --genomeFastaFiles /global/scratch/projects/vector_kvklab/KS-Kronos_remapping/Kronos.collapsed.chromosomes.masked.v1.1.fa --sjdbOverhang 99 --sjdbGTFfile ../../KS-Makido/Processing/65.Final_Final/Kronos.v2.0.gff3 --sjdbGTFtagExonParentTranscript Parent
```
150 bp database: `/global/scratch/projects/vector_kvklab/KS-Kronos_remapping/RNAseqDB/150bp/` and 100 bp database: 

```
for pair1 in $(ls ../../wheat_RNAseq/VanGessel_all_100bp/*.filtered.fastq); do
  prefix=$(basename ${pair1} | cut -d "." -f 1)
  STAR --runThreadN 56 --genomeDir /global/scratch/projects/vector_kvklab/wheat_RNAseq/100bp/  --readFilesIn ${pair1} --outFilterMultimapNmax 5 --outFilterMismatchNoverReadLmax 0.04 --outFileNamePrefix ${prefix}. --outSAMtype BAM SortedByCoordinate --outSAMprimaryFlag AllBestScore
  done
```

Using samtools v1.20, extract primary alignments. We define unique alignments as an alignment without two or more equivalent best score. 
```
samtools flagstats SRR12969761.Aligned.sortedByCoord.out.UniqPrimary.bam
47715226 + 0 in total (QC-passed reads + QC-failed reads)
47715226 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
47715226 + 0 mapped (100.00% : N/A)
47715226 + 0 primary mapped (100.00% : N/A)
47715226 + 0 paired in sequencing
23857613 + 0 read1
23857613 + 0 read2
47715226 + 0 properly paired (100.00% : N/A)
47715226 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
(snp) [skyungyong@n0166 ChenY_all_150bp]$ samtools flagstats SRR12969761.Aligned.sortedByCoord.out.bam
52029741 + 0 in total (QC-passed reads + QC-failed reads)
51742335 + 0 primary
287406 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
52029741 + 0 mapped (100.00% : N/A)
51742335 + 0 primary mapped (100.00% : N/A)
51742335 + 0 paired in sequencing
25871178 + 0 read1
25871157 + 0 read2
51742312 + 0 properly paired (100.00% : N/A)
51742312 + 0 with itself and mate mapped
23 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```


for bam in $(ls *.UniqPrimqry.bam); do
        prefix=$(echo $bam | cut -d "." -f 1 )
        output=${prefix}.featurecounts
        featureCounts -a /global/scratch/projects/vector_kvklab/KS-Kronos_remapping/Kronos.v2.0.gtf -o ${output} -Q 20 --primary -p -C -T 28 --largestOverlap --fracOverlap 0.1 ${bam}

done





