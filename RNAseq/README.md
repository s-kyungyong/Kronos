
For IsoSeq: 

We downloaded available long-read RNA-seq datasets as forms of CCS or HiFi reads with focus on Sequel I and II. There is no datasets directly from Kronos, and these are from *Triticum*. The accessions are available in SRA.list

```
Download fastq files.
```
cat Run.list| while read accession; do
  prefetch "${accession}"
  fasterq-dump -e 56 "${accession}"
done
```

The native Iso-Seq pipeline relies on BAM files, and we only have fastq files. Most of the runs don't seem to include the bam files as well, so we will simply trim adapters and remove poly-A tails if possible and then align those reads with minimap. This is not the most sophisticated way, but no other chocies. 

python filter.py

for fq in *.filtered.fq; do
  prefix=$(echo $fq | cut -d "." -f 1)
  minimap2 -I 12G -t 56 -x splice:hq -a -o "${prefix}.sam" /global/scratch/projects/vector_kvklab/KS-Kronos_remapping/Reference/Kronos.collapsed.chromosomes.masked.v1.1.broken.fa Stringtie/LongReads/"${prefix}.filtered.fq"
  samtools view -@56 -h -b Stringtie/LongReads"${prefix}.sam" | samtools sort -@ 56 > Stringtie/"${prefix}.bam"
done
  samtools merge -@56 all.bam *.bam



For short reads:

hisat2-build -p 56 Kronos.collapsed.chromosomes.masked.v1.1.broken.fa Kronos.broken
