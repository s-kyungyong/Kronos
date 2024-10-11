## NLR analyses

### NLR annotation

Locate the region of genomes in which NB-ARC domains are detected.
```
orfipy --procs 56 --bed orfs.bed --pep orfs.aa.fa --max 45000 --ignore-case --partial-3 --partial-5 Kronos.collapsed.chromosomes.masked.v1.1.fa
hmmsearch --cpu 56 --domE 1e-4 -E 1e-4 --domtblout orfs.against.NBARC.out PF00931.hmm orfs.aa.fa
python crop_genome.py orfs.against.NBARC.out Sohab.Hap1.chromosomes.masked.fa
```

Predict gene structures with MAKER. For protein evidence, NLR sequences for 18 Poaceae species were collected from [this repository](https://zenodo.org/records/13627395) and 415 reference NLRs from [RefPlantNLR](https://zenodo.org/records/3936022). For EST evidence, the transcripts assemembled by Stringtie with short-reads (v1) and long-reads (v2 annotations) were used. The two ab initio parameters obtained in the v1 annotation were used: augustus from braker and snap from ginger. 

```
est=est.fa #this is from the annotation step (de novo + mapping)
protein=proteins.fa #reference NLR sequences
snaphmm=Kronos.hmm
augustus_species=Kronos_collapsed
```

```
mkdir split_genome
cat Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa | grep ">" | cut -d ">" -f 2 | while read tig; do mkdir split_genome/${tig} && awk -v seq=$tig -v RS=">" '$1 ==
seq {print RS $0; exit}' Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa > split_genome/${tig}/${tig}.fa; done
```

for dir in $(ls -d *); do
  cd
  genome=${dir}.fa
  maker -RM_off -genome ${genome} &
  cd ..
done
