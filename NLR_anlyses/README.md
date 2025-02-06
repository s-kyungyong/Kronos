# NLR Curation and Analyses

## NLR Curation

### 1. NLR Loci Isolation
The Kronos genome is large, and for manual curation, NLR loci need to be first extracted from the genome. Locate the region of genomes in which NB-ARC domains are detected, with 15,000 flanking sequences from both ends. These steps will generate an output named *Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa*.
```
# Identify open reading frames (ORFs)
orfipy --procs 56 --bed orfs.bed --pep orfs.aa.fa --max 45000 --ignore-case --partial-3 --partial-5 Kronos.collapsed.chromosomes.masked.v1.1.fa
# Search for NB-ARC domains using HMMER
hmmsearch --cpu 56 --domE 1e-4 -E 1e-4 --domtblout orfs.against.NBARC.out PF00931.hmm orfs.aa.fa
# Crop genome to NB-ARC-containing loci
python crop_genome.py orfs.against.NBARC.out Sohab.Hap1.chromosomes.masked.fa
```

### 2. Initial Annotations
Gene models are needed to aid curation. Typically, if NLRs do not have any mutations that disrupt their gene structures (e.g. framshift mutations or mutations in splicing sites), evidence-based annotators and even ab initio annotators do good jobs predicting the correct structures. Gene structures will be predicted with MAKER v3.01.03. For protein evidence, NLR sequences for 18 Poaceae species were collected from [this repository](https://zenodo.org/records/13627395) and 415 reference NLRs from [RefPlantNLR](https://zenodo.org/records/3936022). For EST evidence, the transcripts assemembled by Stringtie with short-reads (v1) and long-reads (v2 annotations) were used. The two ab initio parameters obtained in the v1 annotation were used: Augustus from BRAKER and SANP from GINGER. 

Default control files from MAKER will be used. In **maker_opts.ctl**, these parameters were modified as below. 
```
est=est.fa                        #this is from the annotation step (de novo + mapping)
protein=proteins.fa               #reference NLR sequences
snaphmm=Kronos.hmm                #from ginger 
augustus_species=Kronos_collapsed #from braker
```

Due to some incompatibility between OpenMPI in our computer cluster and MAKER, each contig will be separated into a file and MAKER will run in parallel. 
```
#separate genome
mkdir split_genome
grep ">" Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa cut -d ">" -f 2 | \
while read tig; do
  mkdir split_genome/${tig}
  awk -v seq=$tig -v RS=">" '$1 == seq {print RS $0; exit}' \
  Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa > split_genome/${tig}/${tig}.fa
done
```
Run MAKER.
```
cd split_genome
for dir in $(ls -d *); do
  cd $dir
  maker -RM_off -genome ${dir}.fa #this basically runs: maker -RM_off -genome ${genome}
  cd ..
done
```

After MAKER completes, collect all the annotations as gff files. 
```
for dir in $(ls -d *); do 
    gff3_merge -d ${dir}/${dir}.maker.output/${dir}_master_datastore_index.log -o ${dir}/${dir}.gff3
done
```

### 3. Collecting Existing Annotations
Other than the MAKER annotations, intermediate annotation files produced during the first annotation (v1) and the two released annotation sets (v1 and v2) are also included.
```
ls
Kronos.v1.0.all.gff3  #version 1 annotation: available through Zenodo
Kronos.v2.0.gff3      #version 2 annotation: available through Zenodo
v1_abinitio.gff3      #includes annotations from BRAKER, Ginger and Funannotate used for EVM 
```

The coordinates in these GFF files need to be adjusted. 
```
less Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa | grep ">" | cut -d ">" -f 2 | sort -u > coordinates.list 
python recoordinate_gff3.py Kronos.v1.0.all.gff3
python recoordinate_gff3.py Kronos.v2.0.gff3
python recoordinate_gff3.py v1_abinitio.gff3
cat Kronos.v1.0.all.recoordinated.gff3 Kronos.v2.0.recoordinated.gff3 v1_abinitio.recoordinated.gff3  > all_models.recoordinated.gff3
```

### 4. Manual curation

perl Apollo/bin/prepare-refseqs.pl --fasta 1A.fa --out .
for feature in {augustus,snap,maker,KRNv1.0,KRNv2.0,v1Annot}; do perl Apollo/bin/flatfile-to-json.pl --trackLabel ${feature} --type mRNA --className mRNA --out . --gff 1A.gff3; done

