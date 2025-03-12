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

