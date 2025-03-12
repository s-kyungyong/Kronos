# *Triticum turgidum* var 'Kronos'


## Preprint

This repository describes the computational pipelines we used to build and analyze our Kronos genome (*Triticum durum* cv Kronos). Please, refer to [out preprint: ]() for more information. 

---
## Databases

Our genome and annotations have been incoporated into multiple databases, which allow visualization of genes and mutations and sequence search. Please, use the databases listed below:

[GrainGenes](https://wheat.pw.usda.gov/GG3/genome_browser)  

[Plant Ensembl]()  

[Dubcovsky Lab](https://dubcovskylab.ucdavis.edu/)  


---
## Data Sharing through Zenodo

Most of the datasets can be accessed through [Zenodo: Chromosome-level genome assembly of Triticum turgidum var 'Kronos'](https://zenodo.org/records/10215402). This repository includes **4** versions, each of which hosts different datasets. Please, use the information below to nevigate the datasets. 

| version                                   | Contents | Comments |
|-------------------------------------------|---------|-----------|
| [1](https://zenodo.org/records/10215402)  | Genome assembly v1.0  |  |
| [2](https://zenodo.org/records/11106422)  | **Genome assembly v1.1** and Annotation v1.0 | Final genome assembly  |
| [3](https://zenodo.org/records/14189805)  | Genome annotation v2.0  |  |
| [4](https://zenodo.org/records/14853918)  | Exome-capture sequencing data remapped for 1,440 Kronos EMS mutants v1.0 |


---
## Final Versions
For clarity, here are the final versions of our datasets

```
Genome assembly:   v1.1   Zenodo:11106422
Genome annoation:  v2.1   Zenodo:         #This is nearly identicial to v2.0 but includes manually curated NLRs. Some v2.0 annotaions were updated.
Repeat annotation: v1.0
Exome capture:     v1.1   Zenodo:         #This is identical to v1.0, but the variant effect prediction was performed on the v2.1 annotation instead of v2.0
promoter capture:  v1.0   Zenodo:
NLR annotation:    v2.1   Zenodo:         #This is included in the v2.1 annotation

```


## Data Availability

### Sequencing Data
Sequencing data were deposted in the NCBI under the BioProject assession, PRJNA1213727. The following runs include: 
```
SRR32063042: HiFi reads
SRR32063043: PacBio sequencing data
SRR32063044: Hi-C sequencing data
```

### Genome Assemblies
The Kronos reference genome can be assessed through the NCBI and Zenodo. To learn more about how these genomes were generated, please refer to [Genome_assembly](https://github.com/s-kyungyong/Kronos/tree/main/Genome_assembly).
```
https://zenodo.org/records/11106422: the Kronos reference genome v1.1** This is the genome we used for our analysis and other databases are hosting. 
https://zenodo.org/records/10215402: the Kronos reference genome v1.0
: the Kronos reference genome v1.1
```

### Genome Annotations 
The Kronos reference genome annotations can be assessed through the NCBI and Zenodo. Our approaches to generate genome annotations v1.0 and v2.0 can be assessed through [Genome_annotations](https://github.com/s-kyungyong/Kronos/tree/main/Genome_annotation).
```
: the Kronos reference genome v2.1** This is the genome annotation we used for our analysis and other databases are hosting. 
https://zenodo.org/records/11106422: the Kronos reference genome v2.0
https://zenodo.org/records/14853918: the Kronos reference genome annotation v1.0
: the Kronos reference genome annotation v2.0
```

Manually curated NLRs are included in the annotation v2.1, updating some existing annotations of v2.0, as recorded in [NLR_analyses](https://github.com/s-kyungyong/Kronos/tree/main/NLR_anlyses). The sequences and labels for NLRs can be accessed through fasta files uploaded in Zenodo. 
```
: the Kronos reference genome v2.1** This is the genome annotation we used for our analysis and other databases are hosting.

```




### Ginger


AUGUSTUS training
compare PASA and BRAKER annotation. Extract the gene models, if they are perfect matches
cp ../../../Braker/braker_rerun/braker.gff3 .
cp ../../../Braker/braker_rerun/braker.aa .
makeblastdb -in trinity_out_dir.Trinity.fasta.transdecoder.pep -out trinity -dbtype 'prot'
blastp -num_threads 40 -evalue 1e-10 -max_target_seqs 1 -max_hsps 1 -outfmt "6 std qlen slen" -out braker_vs_trinity.blast.out -query braker.aa -db trinity

makeblastdb -in ../../../Triticum_aestivum.IWGSC.pep.all.fa -out IWGSC -dbtype 'prot'
blastp -num_threads 40 -evalue 1e-10 -max_target_seqs 1 -max_hsps 1 -outfmt "6 std qlen slen" -out braker_vs_IWGSC.blast.out -query braker.aa -db IWGSC

python select_genes.py
python select_genes.py
Total number of selected gene models: 23785
augustus.gff3 has 8500 genes
snap.gff3 has 8500 genes


gff2gbSmallDNA.pl ../augustus.gff3 ../../../../../3.Repeat/Kronos_output_latest/RepeatMasking/Kronos.collapsed.chromosomes.masked.fa 2000 genes.gb
randomSplit.pl genes.gb 400l
new_species.pl --species=Kronos_manual --AUGUSTUS_CONFIG_PATH=/global/scratch/users/skyungyong/Kronos/5.Annotations/Braker/config
singularity exec -B /global/scratch/users/skyungyong/Kronos/ /global/scratch/users/skyungyong/Kronos/5.Annotations/Braker/braker3.sif etraining --AUGUSTUS_CONFIG_PATH=/global/scratch/users/skyungyong/Kronos/5.Annotations/Braker/config --species=Kronos_manual genes.gb.train
singularity exec -B /global/scratch/users/skyungyong/Kronos/ /global/scratch/users/skyungyong/Kronos/5.Annotations/Braker/braker3.sif augustus --AUGUSTUS_CONFIG_PATH=/global/scratch/users/skyungyong/Kronos/5.Annotations/Braker/config --species=Kronos_manual genes.gb.test | tee first-test.out


cp ../../../../../3.Repeat/Kronos_output_latest/RepeatMasking/Kronos.collapsed.chromosomes.masked.fa genome.dna
samtools view -h -b -f 3 all.merged.sorted.bam > 

g17351.t3
g22576.t1

### Evience modeler 
mkdir abinitio
singularity exec -B /global/scratch/users/skyungyong/Kronos/ EVidenceModeler.v2.1.0.simg perl /usr/local/bin/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl ../Braker/braker_rerun/augustus.hints.gff3 | awk '$1 != "" {$2="augustus1"}1' OFS="\t" | sed 's/model./aug1_model./g' | sed 's/;;/;/g' | sed 's/;.exon/.exon/g' > abinitio/augustus_1.gff3
python ginger_int_to_evmgff.py ../../Ginger/ginger_augustus.gff > temp.gff3 && | awk '$1 != "" {$2="augustus2"}1' OFS="\t" temp.gff3 > augustus_2.gff3
python ginger_int_to_evmgff.py ../../Ginger/ginger_snap.gff > snap.gff3
singularity exec -B /global/scratch/users/skyungyong/Kronos/ ../EVidenceModeler.v2.1.0.simg perl /usr/local/bin/EvmUtils/misc/GeneMarkHMM_GTF_to_EVM_GFF3.pl ../../Braker/braker_rerun/GeneMark-ETP/genemark.gtf > genemark.gff3
python fix_genemark.py genemark.gff3 > temp.gff3 && mv temp.gff3 genemark.gff3
singularity exec -B /global/scratch/users/skyungyong/Kronos/ ../EVidenceModeler.v2.1.0.simg perl /usr/local/bin/EvmUtils/misc/braker_GTF_to_EVM_GFF3.pl ../../Braker/braker_rerun/braker.gtf | awk '$1 != "" {$2="braker"}1' OFS="\t" > braker.gff3
python ginger2evm.py ../../Ginger/ginger_phase2.gff > ginger.gff3
cat abinitio/*gff3 > abinitio.gff3

mkdir homology
singularity exec -B /global/scratch/users/skyungyong/Kronos/ EVidenceModeler.v2.1.0.simg python3 /usr/local/bin/EvmUtils/misc/miniprot_GFF_2_EVM_GFF3.py miniprot.gff3 > miniprot.gff3
python ginger2homology.py ../../Ginger/ginger_homology.gff  > ginger_homology.gff3
cat homology/*.gff3 > homology.gff3

cp ../PASA/sample_mydb_pasa.sqlite.pasa_assemblies.gff3 .
awk '$1 != "" {$2="pasa"}1' OFS="\t" sample_mydb_pasa.sqlite.pasa_assemblies.gff3 > transcript_alignments.gff3

other evidence
 singularity exec -B  /global/scratch/users/skyungyong/Kronos/ pasapipeline.v2.5.3.simg /usr/local/src/PASApipeline/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta sample_mydb_pasa.sqlite.assemblies.fasta --pasa_transcripts_gff3 sample_mydb_pasa.sqlite.pasa_assemblies.gff3

 cp sample_mydb_pasa.sqlite.assemblies.fasta.transdecoder.genome.gff3 ../EVM/others/transdecoder.gff3


V2
#Only use ginger and braker final outputs
cat braker.gff3 ginger.gff3 > ../abinitio.gff3

cp ../V1/transcripts.gff3 .


mkdir others
(base) [skyungyong@n0151 V2]$ cp ../V1/others/transdecoder.gff3 others/
only get complete ones
python get_complete.py
cat transdecoder.complete.gff >> ../abinitio.gff3

mkdir homology && cd homology
ls
Aegilops_tauschii.Aet_v4.0.pep.all.fa                            Hordeum_vulgare_goldenpromise.GPv1.pep.all.fa  Triticum_spelta.PGSBv2.0.pep.all.fa
Avena_sativa_ot3098.Oat_OT3098_v2.pep.all.fa                     Lolium_perenne.MPB_Lper_Kyuss_1697.pep.all.fa  Triticum_turgidum.Svevo.v1.pep.all.fa
Avena_sativa_sang.Asativa_sang.v1.1.pep.all.fa                   Secale_cereale.Rye_Lo7_2018_v1p1p1.pep.all.fa  Triticum_urartu.IGDB.pep.all.fa
Brachypodium_distachyon.Brachypodium_distachyon_v3.0.pep.all.fa  Triticum_aestivum.IWGSC.pep.all.fa
Hordeum_vulgare.MorexV3_pseudomolecules_assembly.pep.all.fa      Triticum_dicoccoides.WEWSeq_v.1.0.pep.all.fa

/global/scratch/users/skyungyong/Software/miniprot/miniprot -t 56 --gff --outc=0.95 -N 0 /global/scratch/users/skyungyong/Kronos/3.Repeat/Kronos_output_latest/RepeatMasking/Kronos.collapsed.chromosomes.masked.fa protein.evidence.fasta > minimap.gff3

 

cat *.fa > protein.evidence.fasta


V3 - Bring the homology information from V1 (uniprot/miniprot + Ginger's alignment)
Weight similar to V1
bring augustus and genemark from braker with low weight.
cp ../V2/abinitio.gff3 .
cat ../V1/abinitio.gff3 | awk '$2 == "augustus1" || $2 == "GeneMark.hmm" {print}'  >> abinitio.gff3
cp ../V2/transcripts.gff3 .
cp ../V1/protein_alignments.gff3 homology.gff3

ABINITIO_PREDICTION     augustus1       1.5
ABINITIO_PREDICTION     GeneMark.hmm    1.5
ABINITIO_PREDICTION     braker  3.5
ABINITIO_PREDICTION     ginger  3.5
PROTEIN homology        2
PROTEIN                  miniprot       1
OTHER_PREDICTION        transdecoder    2
TRANSCRIPT               pasa   7

Overall, V3 looks better, but missing genes. Keep the same input but alter weight

V4
ABINITIO_PREDICTION     augustus1       1.3
ABINITIO_PREDICTION     GeneMark.hmm    1.0
ABINITIO_PREDICTION     braker  4
ABINITIO_PREDICTION     ginger  4
PROTEIN homology        2.3
PROTEIN                  miniprot       1.8
OTHER_PREDICTION        transdecoder    5
TRANSCRIPT               pasa   10

V5

separate ginger_multi ginger_single -> a bit higher for ginger_single
ABINITIO_PREDICTION     funannotate    3
ABINITIO_PREDICTION     braker  3
ABINITIO_PREDICTION     ginger  3
ABINITIO_PREDICTION     gingers  4
PROTEIN homology        1
PROTEIN                  miniprot       1
OTHER_PREDICTION        transdecoder    2
TRANSCRIPT               pasa   8



### NLR annotation

