# NLR manual curation

TrturKRN2B02G099990 is a low-confidence NLR

## Data Availability
Curated NLRs were deposited in [Zenodo](https://zenodo.org/records/15539721).

## Software Versions
```
orfipy v0.0.4
hmmsearch v3.4
nlr-annotator v2.1b
maker v3.01.03
seqkit v
star v2.7.11b
deeptools v3.5.5
samtools v
interproscan v5.68-100.0
apllo v2.0.6
agat v0.8.0
```

---

### 1. Putative NLR Loci Detection
The Kronos genome is large, and for targeted manual curation, NLR loci need to be first extracted from the genome. NLRs were loosely defined as NB-ARC domain-containing genes or proteins. Genomic regions containing NB-ARC domains were identified and extracted with 15,000 flanking sequences from both ends. You may choose to increase the flanking size, as a small number of genes could not be fully contained in this region.

**📥 Inputs**  
• `Kronos.collapsed.chromosomes.masked.v1.1.fa`: Kronos genome  

**📥 Outputs**  
• `Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa`: NLR loci  


⚙️ **identify open reading frames**  
```
orfipy --procs 56 --bed orfs.bed --pep orfs.aa.fa --max 45000 --ignore-case --partial-3 --partial-5 Kronos.collapsed.chromosomes.masked.v1.1.fa
```
⚙️ **Search for NB-ARC domains**  
```
hmmsearch --cpu 56 --domE 1e-4 -E 1e-4 --domtblout orfs.against.NBARC.out PF00931.hmm orfs.aa.fa
```
⚙️ **Extract NLR loci**  
We used 15,000 flanking regions, but some NLRs had really long introns could could not be contained in extracted sequences. It may better to increase this threshold. 
```
python crop_genome.py --hmm orfs.aa.fa.against.NBARC.out --genome Kronos.collapsed.chromosomes.masked.v1.1.fa
```

⚙️ **NLR annotator**  
We later learned that some divergent NB-ARC domains cannot be properly detected by the HMM approach and additionally incoporated NLR-Annotator.
```
java -jar NLR-Annotator-v2.1b.jar -t 40 -x ./NLR-Annotator/src/mot.txt -y ./NLR-Annotator/src/store.txt -i Kronos.collapsed.chromosomes.masked.v1.1.fa -o NLRannotator.whole-genome.out -g NLRannotator.whole-genome.gff3

#crop genome to NB-ARC-containing loci
python crop_genome.py --nlrannot NLRannotator.whole-genome.gff3 --genome Kronos.collapsed.chromosomes.masked.v1.1.fa
```

Alternatively, the two outputs can be combined.
```
python crop_genome.py --hmm orfs.aa.fa.against.NBARC.out --nlrannot NLRannotator.whole-genome.gff3 --genome Kronos.collapsed.chromosomes.masked.v1.1.fa
```
---
### 2. Gene model prediction
Initial gene models were predicted with MAKER. 

**📥 Inputs**  
• `Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa`: NLR loci    
• `proteins.fa`: NLRs from 18 Poaceae species [(source)](https://zenodo.org/records/13627395) and 415 reference NLRs from [RefPlantNLR](https://zenodo.org/records/3936022).  
• `est.fa`: Stringtie transcripts assememblies (v1.0 annotation (short reads) and v2.0 annotation (short/long-reads). See [this record](https://github.com/s-kyungyong/Kronos/tree/main/Genome_annotation).  
• `Kronos_collapsed` & `Kronos.hmm`: The two ab initio parameters obtained in the v1 annotation: Augustus from BRAKER and SANP from GINGER. See [this folder](https://github.com/s-kyungyong/Kronos/tree/main/Genome_annotation/abinitio_parameters))  

**📥 Outputs**  
• `Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.maker_out.gff3`: Initial NLR prediction results  

⚙️ **Run MAKER**  
Default control files from MAKER will be used. In **maker_opts.ctl**, some parameters were modified as below. 
```
est=est.fa                        #est evidence
protein=proteins.fa               #protein evidence
snaphmm=Kronos.hmm                #paramters trained as part of ginger, v1.0 annotation
augustus_species=Kronos_collapsed #parameters trained as part of braker, v1.0 annotation
```

Due to some incompatibility between OpenMPI in our computer cluster and MAKER, each contig was separated and MAKER was run in parallel. 
```
#separate genome
mkdir split_genome
seqkit split -i -O split_genome Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa

#run maker for each directory
maker -RM_off -genome ${dir}.fa

#collect all outputs
for dir in $(ls -d *); do 
    gff3_merge -d ${dir}/${dir}.maker.output/${dir}_master_datastore_index.log -o ${dir}/${dir}.gff3
done

#combine into one file
cat */*.gff3 > Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.maker_out.gff3
```

---
### 3. Additional gene models 
Other than MAKER annotations, intermediate and final annotation files produced during the version 1 and 2 annotations were also included.  

**📥 Inputs**   
• `Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa`: NLR loci    
• `Kronos.v1.0.all.gff3`: version 1.0 annotation  
• `Kronos.v2.0.all.gff3`: version 2.0 annotation  
• `v1_abinitio.gff3`: annotations from BRAKER, Ginger and Funannotate produced as part of v 1.0 annotation. Re-coordinated file is available in this folder.   

**📥 Outputs**    
• `all_models.recoordinated.gff3`: Additional NLR gene models  


⚙️ **Re-coordinate to match NLR loci**  
```
grep ">" Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa | cut -d ">" -f 2 | sort -u > coordinates.list 
python recoordinate_gff3.py Kronos.v1.0.all.gff3
python recoordinate_gff3.py Kronos.v2.0.gff3
python recoordinate_gff3.py v1_abinitio.gff3
cat Kronos.v1.0.all.recoordinated.gff3 Kronos.v2.0.recoordinated.gff3 v1_abinitio.recoordinated.gff3  > all_models.recoordinated.gff3
```

---
### 4. Transcriptome evidence

**📥 Inputs**   
• `Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa`: NLR loci    

**📥 Outputs**  
• `NLR.merged.bam`: transcriptome alignments  
• `NLR.merged.bigwig`: transcriptome coverage  

⚙️ **Run STAR**  
```
#index
STAR --runMode genomeGenerate \
     --genomeFastaFiles Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa \
     --genomeDir GenomeDir

#align transcriptome data from Kronos to the putative NLR loci
#Nnte that if coverage is high, the alignment may not be visuzlied in Apollo genome browser.
#use publicly available RNAseq data (reads.list)

while read -r read1 read2; do
    # Extract the prefix from the first read filename
    prefix=$(basename "$read1" | cut -d "_" -f 1)

    # Run STAR with the specified parameters
    STAR --runThreadN 56 \
        --genomeDir GenomeDir \
        --readFilesIn "$read1" "$read2" \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterMultimapNmax 3 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --outFileNamePrefix "${prefix}." \
        --readFilesCommand zcat
done < reads.list

#filter alignments
for bam in *.bam; do
        samtools view -F 260 -q 20 -@ 56 -b ${bam} > ${bam}.filtered.bam
done

#merge and sort
samtools merge -@ 56 NLR.merged.bam *.filtered.bam
samtools index -@ 56 NLR.merged.bam

#collapse the bam file to a bigwig file.
bamCoverage -b NLR.merged.bam -o NLR.merged.bigwig --binSize 1 --normalizeUsing None
```

---
### 5. Domain annotation
Lastly, domain prediction on 6-frame translated genomic sequenes can be loaded.

**📥 Inputs**   
• `Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa`: NLR loci     

**📥 Outputs**    
• `Kronos.aa.fa.gff3`: Domain prediction on translated genomic sequences  

⚙️ **Predict domains across genomic sequences**  
```
#identify open reading frames (ORFs)
orfipy --procs 56 --bed orfs.bed --pep orfs.aa.fa --max 45000 --ignore-case --partial-3 --partial-5 Kronos.collapsed.chromosomes.masked.v1.1.fa

#predict pfam domains
interproscan.sh \
    -i orfs.aa.fa \
    --appl PFAM-37.0 \
    --disable-precalc \
    --output-dir ./output \
    --tempdir ./temp \
    --cpu 40 \
```

---
### 6. Manual curation
For manual curation, all these datasets need to be loaded into the Apollo Genome Browser. We conducted annotations for each chromosome. It took about 3 weeks to annotate all loci. After curation, there was a lot of downstream QC and manual examination over multiple iterations.

**📥 Inputs**   
• `Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa`: NLR loci     
• `Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.maker_out.gff3`: Initial NLR prediction results  
• `all_models.recoordinated.gff3`: Additional NLR gene models  
• `NLR.merged.bam`: transcriptome alignments  
• `NLR.merged.bigwig`: transcriptome coverage  


⚙️ **Load into genome browser**  
```
#create reference sequence
perl ./Apollo/bin/prepare-refseqs.pl --fasta Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa --out .

#load annotations
#first, extract genes in that chromosome
cat all_models.recoordinated.gff3 Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.maker_out.gff3 | \
    awk '$1 == ${chromosome} && ($3 == "gene" || $3 == "mRNA" || $3 == "exon" || $3 == "CDS") {print}' > all.genes.gff3

#get rid of models with in-frame stop codon using agat and then load
agat_sp_flag_premature_stop_codons.pl --gff all.genes.gff3 --fa Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa --out all.genes.fixed.gff3
perl ./Apollo/bin/flatfile-to-json.pl --trackLabel models --type mRNA --className mRNA --out . --gff all.genes.fixed.gff3

#change and load interproscan coordinates
python recoordinate_ipr.py orfs.aa.fa.gff3 all.nlr_loci.fa ${chromosome} > ${chromosome}.ipr.gff3
perl ./Apollo/bin/flatfile-to-json.pl --trackLabel iprscan --type match:iprscan --out . --gff ${chromosome}.ipr.gff3

#load transcriptome data
perl ./Apollo/bin/add-bam-track.pl --bam_url NLR.merged.bam --label bam --in ./trackList.json
perl ./Apollo/bin/add-bw-track.pl --bw_url NLR.merged.bigwig --label bigwig --in ./trackList.json
```

