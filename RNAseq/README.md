
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



gffcompare -r Kronos.v1.0.all.recoordinated.gtf -o Kronos isoquant.all.split.gtf  stringtie.denovo.all.split.renamed.gtf stringtie.guided.all.split.renamed.gtf
/global/scratch/users/skyungyong/Software/TransDecoder-TransDecoder-v5.7.1/util/gtf_genome_to_cdna_fasta.pl Kronos.combined.repositioned.gtf /global/scratch/users/skyungyong/Kronos/5.Annotations/Final/Final_Final_for_release/Kronos.collapsed.chromosomes.masked.v1.1.fa > transcripts.fa

/global/scratch/users/skyungyong/Software/TransDecoder-TransDecoder-v5.7.1/util/gtf_to_alignment_gff3.pl Kronos.combined.repositioned.gtf > transcripts.gff3

/global/scratch/users/skyungyong/Software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t transcripts.fa

# extract complete
python extract_complete.py longest_orfs.pep > longest_orfs.complete.pep

#plant ensemble only
diamond makedb --threads 40 --db all.prot.evidence.fa --in all.prot.evidence.fa
diamond blastp --threads 40 --db db/all.prot.evidence.fa --evalue 1e-10 --max-hsps 1 --max-target-seqs 1 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --out longest.against.db.dmnd.out --query longest_orfs.complete.pep
awk '($8 - $7 + 1)/($13 - 1) > 0.97 && ($10 - $9 + 1)/($14) > 0.97 {print}' longest.against.db.dmnd.out > longest.against.db.dmnd.out.filtered
/global/scratch/users/skyungyong/Software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict --retain_blastp_hits transcripts.fa.transdecoder_dir/longest.against.db.dmnd.out.filtered -t transcripts.fa

/global/scratch/users/skyungyong/Software/TransDecoder-TransDecoder-v5.7.1/util/cdna_alignment_orf_to_genome_orf.pl transcripts.fa.transdecoder.gff3 transcripts.gff3 transcripts.fa > transcripts.fa.transdecoder.genome.gff3


Now refine the annotation with maker
from uniprot: Poaceae AND (taxonomy_id:4479)


#do it for three 
python /global/scratch/users/skyungyong/Software/BRAKER/scripts/stringtie2fa.py -f stringtie.denovo.all.split.renamed.gtf -g /global/scratch/users/skyungyong/Kronos/5.Annotations/Final/Final_Final_for_release/Kronos.collapsed.chromosomes.masked.v1.1.fa -o stringtie.denovo.all.split.renamed.fa

/global/scratch/users/skyungyong/Software/gmst.pl --strand both stringtie.denovo.all.split.renamed.fa.mrna --output stringtie.denovo.all.split.renamed.fa.mrna.gmst.out --format GFF

python /global/scratch/users/skyungyong/Software/BRAKER/scripts/gmst2globalCoords.py -t stringtie.denovo.all.split.renamed.repositioned.gtf -p stringtie.denovo.all.split.renamed.repositioned.no_Ns.fa.mrna.gmst.out -o stringtie.denovo.gmst.global.gtf -g /global/scratch/users/skyungyong/Kronos/5.Annotations/Final/Final_Final_for_release/Kronos.collapsed.chromosomes.masked.v1.1.fa
