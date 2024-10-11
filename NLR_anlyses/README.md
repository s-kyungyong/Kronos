## NLR analyses

### NLR annotation

Locate the region of genomes in which NB-ARC domains are detected.
```
orfipy --procs 56 --bed orfs.bed --pep orfs.aa.fa --max 45000 --ignore-case --partial-3 --partial-5 Kronos.collapsed.chromosomes.masked.v1.1.fa
hmmsearch --cpu 56 --domE 1e-4 -E 1e-4 --domtblout orfs.against.NBARC.out PF00931.hmm orfs.aa.fa
python crop_genome.py orfs.against.NBARC.out Sohab.Hap1.chromosomes.masked.fa
```

Predict gene structures with MAKER. For protein evidence, NLR sequences for 18 Poaceae species were collected from (this repository)[https://zenodo.org/records/13627395].
```


```
