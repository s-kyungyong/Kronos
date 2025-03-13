



Quantification

```
#prepare database
grep ">"Kronos.collapsed.chromosomes.masked.v1.1.fa | cut -d " " -f 1 | cut -d ">" -f 2 > decoys.txt
gffread -x Kronos.v2.1.transcripts.fa -g Kronos.collapsed.chromosomes.masked.v1.1.fa Kronos.v2.1.gff3
cat Kronos.v2.1.transcripts.fa Kronos.collapsed.chromosomes.masked.v1.1.fa > Kronos.gentrome.fa

#index gentrome
salmon index -t Kronos.gentrome.fa -d decoys.txt -p 30 -i salmon_index
```

```
#map and quantify:
#for paired-end data:

indir=$1

for fq1 in ${indir}/*_1.filtered.fastq; do
  prefix=$(basename "$fq1" | cut -d "_" -f 1) 
  fq2=$(echo $fq | sed 's/_1.filtered/_2.filtered/g')

  salmon quant -l A -1 "$fq1" -2 "$fq2" -p 40 \
    -g /global/scratch/projects/vector_kvklab/KS-Kronos_remapping/RNAseqDB-Salmon/Kronos.v2.1.gtf \
    -i /global/scratch/projects/vector_kvklab/KS-Kronos_remapping/RNAseqDB-Salmon/salmon_index/ \
    -o "${prefix}" --validateMappings
done

#for single_end
indir=$1

for fq1 in ${indir}/*.filtered.fastq; do
  prefix=$(basename "$fq1" | cut -d "." -f 1) 

  salmon quant -l A -r "$fq1" -p 40 \
    -g /global/scratch/projects/vector_kvklab/KS-Kronos_remapping/RNAseqDB-Salmon/Kronos.v2.1.gtf \
    -i /global/scratch/projects/vector_kvklab/KS-Kronos_remapping/RNAseqDB-Salmon/salmon_index/ \
    -o "${prefix}" --validateMappings
done
```







minimap2 -x asm5 -t 20 ../../Triticum_aestivum.plasmids.fa Kronos.draft.fa > min
imap.plasmid.paf
sort -k1,1 -k3,3n minimap.plasmid.paf > minimap.plasmid.sorted.paf
python ../process_scaffolds.py ../minimap.plasmid.sorted.paf None Kronos.draft.fa Kronos.contigs
quast -t 20 --fast Kronos.contigs.genomic.fa


Repeat annotation: temp
singularity exec -B $(pwd):$(pwd) /global/scratch/users/skyungyong/Software/EDTA.sif EDTA.pl --genome /global/scratch/users/skyungyong/Kronos/5.Annotations/Final/Kronos.collapsed.chromosomes.masked.v1.1.fa --species others --step all --sensitive 1 --anno 1 --evaluate 1 --threads 56 --cds Kronos.v2.0.cds.fa --rmlib

./tRNAscan-SE_installed/bin/tRNAscan-SE -E -o tRNAscan-SE.out -f tRNAscan-SE.ss -s tRNAscan-SE.iso -m tRNAscan-SE.stats -c ./tRNAscan-SE_installed/bin/tRNAscan-SE.conf ../Final/Kronos.collapsed.chromosomes.v1.1.fa

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
for lib in $(ls SRX10965365*_val_1.fq); do
  prefix=$(echo $lib | cut -d "_" -f 1)
  read_1=$lib
  read_2=$(echo $lib | sed 's/1_val_1/2_val_2/')
  hisat2 -p 56 -x /global/scratch/projects/vector_kvklab/KS-Kronos_remapping/Reference/Kronos.broken -1 $read_1 -2 $read_2 --no-discordant --no-mixed --dta -S /global/scratch/projects/vector_kvklab/KS-IsoSeq-HiFi/Stringtie/ShortReads/${prefix}.sam
done


stringtie -p 56 -v --mix LongReads/all.sorted.bam ShortReads/all.merged.bam


python merge_gtf.py
gffread -w PostProcessing/Kronos.transcripts.fa -g ../Kronos.collapsed.chromosomes.masked.fa Kronos.all.gtf
TransDecoder.LongOrfs -t Kronos.transcripts.fa
/usr/local/bin/diamond makedb --threads 40 --db all.evidence.fa --in all.evidence.fa #v2.1.9.163

/usr/local/bin/diamond blastp --threads 40 --db ../../database/all.evidence.fa --evalue 1e-10 --max-hsps 1 --max-target-seqs 1 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --out longest.against.db.dmnd.out --query longest_orfs.pep



gffread -w Kronos.stringtie.transcripts.fa -g /global/scratch/users/skyungyong/Kronos/5.Annotations/Final/Final_Final_for_release/Kronos.collapsed.chromosomes.masked.v1.1.fa Kronos.stringtie.all.gtf

__
mamba crate -n isoquant python=3.8 isoquant
activate isoquant
python filter.py

fq=$(ls *.filtered.fq)
isoquant.py --threads 56 --reference /global/scratch/projects/vector_kvklab/KS-Kronos_remapping/Reference/Kronos.collapsed.chromosomes.masked.v1.1.broken.fa --genedb /global/scratch/projects/vector_kvklab/KS-IsoSeq-HiFi/Stringtie/Guided/Kronos.v1.0.all.recoordinated.gtf --illumina_bam /global/scratch/projects/vector_kvklab/KS-IsoSeq-HiFi/Stringtie/ShortReads/all.merged.bam --output Isoquant_Kronos --data_type pacbio_ccs --fastq $fq


# now testing stringtie denovo only
gffcompare -r Kronos.v1.0.all.gff3 -o Kronos isoquant.all.split.repositioned.gtf  stringtie.denovo.all.split.renamed.repositioned.gtf stringtie.guided.all.split.renamed.repositioned.gtf
/global/scratch/users/skyungyong/Software/TransDecoder-TransDecoder-v5.7.1/util/gtf_genome_to_cdna_fasta.pl Kronos.combined.gtf
 /global/scratch/users/skyungyong/Kronos/5.Annotations/Final/Final_Final_for_release/Kronos.collapsed.chromosomes.masked.v1.1.fa > transcripts.fa

/global/scratch/users/skyungyong/Software/TransDecoder-TransDecoder-v5.7.1/util/gtf_to_alignment_gff3.pl Kronos.combined.gtf > transcripts.gff3

/global/scratch/users/skyungyong/Software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t transcripts.fa

# extract complete
python extract_complete.py longest_orfs.pep > longest_orfs.complete.pep

#plant ensemble only
diamond makedb --threads 40 --db all.prot.evidence.fa --in all.prot.evidence.fa
diamond blastp --threads 40 --db db/all.prot.evidence.fa --evalue 1e-10 --max-hsps 1 --max-target-seqs 1 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --out longest.against.db.dmnd.out --query longest_orfs.complete.pep
awk '($8 - $7 + 1)/($13 - 1) > 0.97 && ($10 - $9 + 1)/($14) > 0.97 {print}' longest.against.db.dmnd.out > longest.against.db.dmnd.out.filtered
/global/scratch/users/skyungyong/Software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict --retain_blastp_hits transcripts.fa.transdecoder_dir/longest.against.db.dmnd.out.filtered -t transcripts.fa

/global/scratch/users/skyungyong/Software/TransDecoder-TransDecoder-v5.7.1/util/cdna_alignment_orf_to_genome_orf.pl transcripts.fa.transdecoder.gff3 transcripts.gff3 transcripts.fa > transcripts.fa.transdecoder.genome.gff3




gffcompare -r ../Kronos.v1.0.all.gff3 ../Stringtie_denovo_transdecoder/transcripts.fa.transdecoder.genome.complete_only.gff3
gffread -F -M -d gffcmp.duplicates -K ../Kronos.v1.0.all.gff3 ../Stringtie_denovo_transdecoder/transcripts.fa.transdecoder.genome.complete_only.gff3 > gffread.out

  172203 reference transcripts loaded.
  25 duplicate reference transcripts discarded.
  227208 query transfrags loaded.

These genes need to be fixed
['TrturKRN1A01G039440', 'TrturKRN1A01G039430']
['TrturKRN1A01G039460', 'TrturKRN1A01G039430']
['TrturKRN1A01G039480', 'TrturKRN1A01G039430']
['TrturKRN1B01G019200', 'TrturKRN1B01G019210']
['TrturKRN1B01G028300', 'TrturKRN1B01G028290']
['TrturKRN1B01G032370', 'TrturKRN1B01G032360']
['TrturKRN1B01G057670', 'TrturKRN1B01G057680']
['TrturKRN1B01G069790', 'TrturKRN1B01G069780']
['TrturKRN2B01G070470', 'TrturKRN2B01G070460']
['TrturKRN3B01G013460', 'TrturKRN3B01G013440', 'TrturKRN3B01G013430']
['TrturKRN3B01G013480', 'TrturKRN3B01G013430']
['TrturKRN4A01G043880', 'TrturKRN4A01G043890']
['TrturKRN6A01G013540', 'TrturKRN6A01G013530']
['TrturKRN7A01G004660', 'TrturKRN7A01G004650']
['TrturKRN7A01G073810', 'TrturKRN7A01G073800']

1) Identical genes: "=" or "c" : 20179
2) Nothing predicted in V2: 44962
3) Unique complete genes in V2: 4643
4) Low -> High: 





 gffread -g /global/scratch/users/skyungyong/Kronos/5.Annotations/Final/Final_Final_for_release/Kronos.collapsed.chromosomes.masked.v1.1.fa -x Kronos.cds.fa -w Kronos.exon.fa -y Kronos.pep.fa ../transcripts.fa.transdecoder.genome.gff3

less ../Kronos.tracking.coding.only | awk '{print $2 "\t" $4}' | sort | uniq -c > Kronos.tracking.summary
less Kronos.tracking.summary | awk '{print $2}' | sort | uniq -c > Kronos.tracking.summary.gene.counts


1) Identical genes: no changes in V1 and V2
high: 19034
low:  11902

1) Version 2 unique transcripts
in gffcompare, all transcripts of a gene are 'u' (unknown) as there are no overlapping genes. One of these transripts have CDS_type=complete and protein=complete
python select_novel_genes.py > add_new_genes.list
== 4596 new genes were added

2) Improved in Version 2
classified as 'low quality' in V1, but in V2, they are complete genes
== 6779 genes improved


repeat annotation with EDTA

wget https://trep-db.uzh.ch/downloads/trep-db_complete_Rel-19.fasta.gz

For the annotation of noncoding RNAs, tRNAscan-SE software52 was used to predict the tRNAs with eukaryotic parameters. miRNAs, rRNAs, and snRNAs were detected using Infernal cmscan53 to search the Rfam database54. The rRNAs and the corresponding subunits were annotated with RNAmmer v1.255.
