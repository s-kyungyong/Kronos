# Functional Annotations

## Data Availability
The functional annotations can be downloaded from Zenodo.   
• `Functional annotations for annotation v2.1` [✨Final✨]: https://zenodo.org/records/15539216  

## Software version
```
eggNOG-mapper v2.1.12
InterProScan v5.68.100
```

---


## 1. eggNOG-mapper

📥 Inputs
• `Kronos.v2.1.pep.fa`: Kronos annotation v2.1  

📥 Outputs
• `Kronos.v2.1.pep.eggnog.tsv`: eggNOG-mapper annotation

⚙️ Run eggNOG-mapper
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


## 2. InterProScan

📥 Inputs
• `Kronos.v2.1.pep.fa`: Kronos annotation v2.1  

📥 Outputs  
• `Kronos.v2.1.pep.fa.InterPro.tsv`: InterProScan annotations  

⚙️ Run InterProScan

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

## 3. Annotation Transfer

📥 Inputs
• `Kronos.v2.1.pep.fa`: Kronos annotation v2.1  
• `RefSeqs.fa`: Selected RefSeq sequences from NCBI
• `Kronos.v2.1.pep.eggnog.tsv`: eggNOG-mapper annotation

📥 Outputs  
• `Kronos.v2.1.description`: Gene descriptions for Kronos annotation v2.1

⚙️ Transfer Gene Annotation
• Collect databases
```
#download proteins
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/294/505/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/162/155/GCF_002162155.2_WEW_v2.1/GCF_002162155.2_WEW_v2.1_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/073/215/GCF_003073215.2_Tu2.1/GCF_003073215.2_Tu2.1_protein.faa.gz
gunzip *.gz 

#concatnate all proteins
cat *_protein.faa > NCBI_refseq.aa.fa 
```

• Homology search
```
diamond makedb --in NCBI_refseq.aa.fa --db NCBI_refseq.aa
diamond blastp --masking 0 -d NCBI_refseq.aa --evalue 1-e04 --max-target-seqs 5 \
        --query Kronos.v2.1.pep.fa --out Kronos.v2.1.against.NCBI.refseq.dmnd.out \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen 
```

• Annotation transfer
For annotation transfer, query and hit coverage > 90%, sequence ideneity > 90% and E-value < 1e-10 were required. For sequences that fail to meet these criteria, annotations were lifted from eggNOG-mapper if present. Otherwise, sequences were annotated as hypothetical proteins. 
```
python transfer_annotations.py
```
