# Functional Annotations

## eggNOG-mapper

Functional annotations are from eggNOG-MAPPER v2.1.12. All the parameters were set default, and the sequences were submitted to [the eggNOG-mapper server](http://eggnog-mapper.embl.de/).
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

Domain annotations are from InterProScan. 
```
nohup docker run --rm \
    -v /Users/s.kyungyong/Desktop/Software/interproscan-5.68-100.0/data:/opt/interproscan/data \
    -v $PWD/output:/output \
    -v $PWD:/workspace \
    -v $PWD/temp:/temp \
    interpro/interproscan:5.68-100.0 \
    --input /workspace/Kronos.v2.1.pep.fa \
    --disable-precalc \
    --output-dir /output \
    --tempdir /temp \
    --cpu 48 \
    --appl FunFam-4.3.0,SFLD-4,Gene3D-4.3.0,SUPERFAMILY-1.75,SMART-9.0,CDD-3.20,PIRSR-2023_05,Pfam-37.0,PIRSF-3.10,NCBIfam-14.0 &
```
