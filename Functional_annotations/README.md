# Functional Annotations

## Data Availability
The functional annotations can be downloaded from Zenodo. 
```
https://zenodo.org/records/15539216: Functional annotations for annotation v2.1
```

## Methods

Functional annotations were performed using eggNOG-mapper v2.1.12 and InterProScan v5.68.100. 

## Software version
```
eggNOG-mapper v2.1.12
InterProScan v5.68.100
```
---


## eggNOG-mapper

All the parameters were set default, and the sequences were submitted to [the eggNOG-mapper server](http://eggnog-mapper.embl.de/).
```
Minimum hit e-value:  0.001
Minimum hit bit-score:  60
Percentage identity:  40
Minimum % of query coverage:  20
Minimum % of subject coverage:  20

Taxonomic Scope: Auto adjust per query
Orthology restrictions: Transfer annotations from any ortholog
Gene Ontology evidence: Transfer non-electronic annotations
PFAM refinement: Report PFAM domains from orthologs
SMART annotation: Skip SMART annotations
```


## InterProScan

Due to memory issues, each database was searched and then all outputs were concatnated. 
```
DATABASES="FunFam-4.3.0 SFLD-4 PANTHER-18.0 Gene3D-4.3.0 Hamap-2023_05 PRINTS-42.0 ProSiteProfiles-2023_05 Coils-2.2.1 SUPERFAMILY-1.75 SMART-9.0 CDD-3.20 PIRSR-2023_05 ProSitePatterns-2023_05 AntiFam-7.0 Pfam-37.0 MobiDBLite-2.0 PIRSF-3.10 NCBIfam-14.0"

for db in $DATABASES; do
    prefix=$(echo $db | cut -d "-" -f 1)
    mkdir -p $prefix  # Create directory if it doesn't exist

    /global/scratch/users/skyungyong/Software/interproscan-5.68-100.0/interproscan.sh \
        -i Kronos.v2.1.pep.fa \
        --disable-precalc \
        --output-dir $prefix \
        --cpu 16 \
        --goterms \
        --appl $db  # Run only the specific database
done
```

The commandline above enables the following analyses
```
FunFam (4.3.0) : Prediction of functional annotations for novel, uncharacterized sequences.
SFLD (4) : SFLD is a database of protein families based on hidden Markov models (HMMs).
PANTHER (18.0) : The PANTHER (Protein ANalysis THrough Evolutionary Relationships) Classification System is a unique resource that classifies genes by their functions, using published scientific experimental evidence and evolutionary relationships to predict function even in the absence of direct experimental evidence.
Gene3D (4.3.0) : Structural assignment for whole genes and genomes using the CATH domain structure database.
Hamap (2023_05) : High-quality Automated and Manual Annotation of Microbial Proteomes.
PRINTS (42.0) : A compendium of protein fingerprints - a fingerprint is a group of conserved motifs used to characterise a protein family.
ProSiteProfiles (2023_05) : PROSITE consists of documentation entries describing protein domains, families and functional sites as well as associated patterns and profiles to identify them.
Coils (2.2.1) : Prediction of coiled coil regions in proteins.
SUPERFAMILY (1.75) : SUPERFAMILY is a database of structural and functional annotations for all proteins and genomes.
SMART (9.0) : SMART allows the identification and analysis of domain architectures based on hidden Markov models (HMMs).
CDD (3.20) : CDD predicts protein domains and families based on a collection of well-annotated multiple sequence alignment models.
PIRSR (2023_05) : PIRSR is a database of protein families based on hidden Markov models (HMMs) and Site Rules.
ProSitePatterns (2023_05) : PROSITE consists of documentation entries describing protein domains, families and functional sites as well as associated patterns and profiles to identify them.
AntiFam (7.0) : AntiFam is a resource of profile-HMMs designed to identify spurious protein predictions.
Pfam (37.0) : A large collection of protein families, each represented by multiple sequence alignments and hidden Markov models (HMMs).
MobiDBLite (2.0) : Prediction of intrinsically disordered regions in proteins.
PIRSF (3.10) : The PIRSF concept is used as a guiding principle to provide comprehensive and non-overlapping clustering of UniProtKB sequences into a hierarchical order to reflect their evolutionary relationships.
NCBIfam (14.0) : NCBIfam is a collection of protein families based on Hidden Markov Models (HMMs).
```


## Annotation Transfer
We also added descriptions for each protein through annotation transfer. These annotations come from RefSeqs curated by the NCBI. For annotation transfer, query and hit coverage > 90%, sequence ideneity > 90% and E-value < 1e-10 were required. For sequences that fail to meet these criteria, annotations were lifted from eggNOG-mapper if present. Otherwise, sequences were annotated as hypothetical proteins. 

```
#download proteins
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/294/505/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/162/155/GCF_002162155.2_WEW_v2.1/GCF_002162155.2_WEW_v2.1_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/073/215/GCF_003073215.2_Tu2.1/GCF_003073215.2_Tu2.1_protein.faa.gz
gunzip *.gz 

#concatnate all proteins
cat *_protein.faa > NCBI_refseq.aa.fa 

#run dimaond
diamond makedb --in NCBI_refseq.aa.fa --db NCBI_refseq.aa
diamond blastp --masking 0 -d NCBI_refseq.aa --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --evalue 1-e04 --max-target-seqs 5 --query Kronos.v2.1.pep.fa --out Kronos.v2.1.against.NCBI.refseq.dmnd.out

#annottion transfer
python transfer_annotations.py
```
