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
```
