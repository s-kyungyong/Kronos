
## Protein-coding Gene Preidction: v1.0 annotation

The first version of genome annotation largely focused on the integration of short-read sequencing data produced for Kronos. 

### Paired-end Short-read Transcriptome Data Processing

We donwloaded the paired-end RNA-seq data from the NCBI. The list can be found in [v1_rnaseq.list](https://github.com/s-kyungyong/Kronos/blob/main/RNAseq/SRA.list](https://github.com/s-kyungyong/Kronos/blob/main/Genome_annotation/v1_rnaseq.list). 
```
while read -r accession; do 
    sratoolkit.3.1.1-centos_linux64/bin/prefetch ${accession}
    sratoolkit.3.1.1-centos_linux64/bin/fasterq-dump -O . -e ${Numthreads} ${accession}
done < v1_rnaseq.list
```
Let's first remove adapters and low-quality reads from the libraries. Trim_galore and cutadapt versions were v0.6.6 and v3.7. This generated about 1.6 Tb of fastq files in total
```
ls *.fastq | cut -d "_" -f 1 | sort -u | while read accession; do trim_galore -a -j 8 --paired $accession\_1.fastq $accession\_2.fastq ; done
```

The RNA-seq data can be mapped to the genome and processed. Reads will be mapped to the genome with hisat, filtered with samtools and assembled with stringtie. 
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
```


De novo assembly is done with Trinity v2.15.1. We initially tried running Trinity on the 1.6 Tb of paired-end fastq files all at once. After two weeks, Trinity was still stuck at the insilico normalization step with about 35-45% progress. We, therefore, had to take some other ways around. Each pair will be normalized first, and then Trinity is run. This took a few days, producing transcripts of ~1 Gb. 
```
#for each pair 
singularity run -B $PWD trinity.sif Trinity --verbose --max_memory 90G --just_normalize_reads --seqType fq --CPU 40 --left $left --right $right --output trinity_$prefix

#list all normalized reads
ls -d trinity_* | while read folder; do
    prefix=$(echo $folder | cut -d "_" -f 2)
    left=$(ls $(pwd)/$folder\/insilico_read_normalization/*_1_val_1*.fq)
    right=$(ls $(pwd)\/$folder\/insilico_read_normalization/*_2_val_2*.fq)
    echo $prefix $prefix $left $right
done > sample.list

#run trinity
singularity run -B $PWD trinity.sif Trinity  --verbose --seqType fq --max_memory 1500G --CPU 56 --samples_file $PWD/sample.list
```


### BRAKER

To run BRAKER, we need transcriptome and protein evidence. **all.merged.sorted.bam** produced in the pipeline above will be the transcriptome evidence. For protein evidence, 2,850,097 sequences that belong to Poales (TAXID: 38820) were downloaded from UniProt. Then BRAKER can be run:

```
singularity exec -B $PWD braker3.sif braker.pl --verbosity=3 \
    --genome=Kronos.collapsed.chromosomes.masked.fa \
    --bam=all.merged.sorted.bam \
    --prot_seq=uniprotkb_taxonomy_id_38820_2023_12_08.fasta \
    --species=Kronos --threads 48 --gff3 \
    --workingdir=$wd/braker \
    --AUGUSTUS_CONFIG_PATH=$wd/config
```

### Funannotate
```
funannotate predict \
-i /global/scratch/users/skyungyong/Kronos/3.Repeat/Kronos_output_latest/RepeatMasking/Kronos.collapsed.chromosomes.masked.fa \
-o /global/scratch/users/skyungyong/Kronos/5.Annotations/Funannotate \
-s "Triticum kronos" \
--transcript_evidence /global/scratch/users/skyungyong/Kronos/5.Annotations/PASA/transcripts.fasta \
--repeats2evm \
--cpus 56 \
--ploidy 2 \
--rna_bam /global/scratch/users/skyungyong/Kronos/5.Annotations/Ginger/alignments/all.merged.ginger.sorted.bam \
--stringtie /global/scratch/users/skyungyong/Kronos/5.Annotations/Stringtie/stringtie.gtf \
--augustus_species Kronos_maker \
--AUGUSTUS_CONFIG_PATH /global/scratch/users/skyungyong/Software/anaconda3/envs/funannotate/config/ \
--organism other \
--EVM_HOME /global/scratch/users/skyungyong/Software/anaconda3/envs/funannotate/opt/evidencemodeler-1.1.1/ \
--GENEMARK_PATH /global/scratch/users/skyungyong/Software/gmes_linux_64_4 \
```

### Miniprot

```
/global/scratch/users/skyungyong/Software/miniprot/miniprot -t 56 --gff --outc=0.95 -N 0 ../3.Repeat/Kronos_output_latest/Kronos.collapsed.chromosomes.fa ../5.Annotations/Braker/uniprotkb_taxonomy_id_38820_2023_12_08.fasta
```



### EvidenceModler

singularity exec -B /global/scratch/users/skyungyong/Kronos/ /global/scratch/users/skyungyong/Kronos/5.Annotations/EVM/EVidenceModeler.v2.1.0.simg EVidenceModeler --sample_id Kronos --genome $(pwd)/genome.fa --weights /global/scratch/users/skyungyong/Kronos/5.Annotations/EVM/V5/EVM_outputs/weights.txt --gene_predictions $(pwd)/abinitio.gff3 --protein_alignments $(pwd)/homology.gff3 --transcript_alignments $(pwd)/transcripts.gff3 --repeats $(pwd)/repeat.gff3 --CPU 56 -S --segmentSize 100000 --overlapSize 10000


ABINITIO_PREDICTION     funannotate   3
ABINITIO_PREDICTION     braker  3
ABINITIO_PREDICTION     ginger  3
ABINITIO_PREDICTION     gingers  3.5
PROTEIN homology        2
PROTEIN                  miniprot       1
OTHER_PREDICTION        transdecoder    2.5
TRANSCRIPT               pasa  8

