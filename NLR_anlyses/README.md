## Genome assembly and assessment

Let's use the same filtered input files and perform genome assembly with them.

```
orfipy --procs 56 --bed orfs.bed --pep orfs.aa.fa --max 45000 --ignore-case --partial-3 --partial-5 Kronos.collapsed.chromosomes.masked.v1.1.fa
hmmsearch --cpu 56 --domE 1e-4 -E 1e-4 --domtblout orfs.against.NBARC.out PF00931.hmm orfs.aa.fa
python crop_genome.py orfs.against.NBARC.out Sohab.Hap1.chromosomes.masked.fa
```
