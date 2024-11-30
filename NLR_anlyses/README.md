## NLR Curation and Analyses

### NLR curation

To manually curate NLR gene models, we need to prepare some initial annotations. The goal here is to isolate genomic regions that contain the NB-ARC domain, create annotations with various tools and manually select/modify those. 

Locate the region of genomes in which NB-ARC domains are detected. We will use the genome v1.1. 
```
orfipy --procs 56 --bed orfs.bed --pep orfs.aa.fa --max 45000 --ignore-case --partial-3 --partial-5 Kronos.collapsed.chromosomes.masked.v1.1.fa
hmmsearch --cpu 56 --domE 1e-4 -E 1e-4 --domtblout orfs.against.NBARC.out PF00931.hmm orfs.aa.fa
python crop_genome.py orfs.against.NBARC.out Sohab.Hap1.chromosomes.masked.fa
```

Predict gene structures with MAKER v3.01.03. For protein evidence, NLR sequences for 18 Poaceae species were collected from [this repository](https://zenodo.org/records/13627395) and 415 reference NLRs from [RefPlantNLR](https://zenodo.org/records/3936022). For EST evidence, the transcripts assemembled by Stringtie with short-reads (v1) and long-reads (v2 annotations) were used. The two ab initio parameters obtained in the v1 annotation were used: augustus from braker and snap from ginger. 

We will use default control files from MAKER. In **maker_opts.ctl**, these parameters were modified as below. 
```
est=est.fa #this is from the annotation step (de novo + mapping)
protein=proteins.fa #reference NLR sequences
snaphmm=Kronos.hmm #from ginger 
augustus_species=Kronos_collapsed #from braker
```

Due to some incompatibility between OpenMPI in our computer cluster and MAKER, each contig will be separated into a file and MAKER will run in parallel. 
```
#separate genome
mkdir split_genome
cat Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa | grep ">" | cut -d ">" -f 2 | while read tig; do mkdir split_genome/${tig} && awk -v seq=$tig -v RS=">" '$1 == seq {print RS $0; exit}' Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa > split_genome/${tig}/${tig}.fa; done

#run MAKER
cd split_genome
for dir in $(ls -d *); do
  cd $dir
  genome=${dir}.fa
  sbatch ../run_maker.sh ${genome} #this basically runs: maker -RM_off -genome ${genome}
  cd ..
done
```

Once done, let's collect the outputs. We will do this for each chromosome. 
```
for dir in $(ls -d *); do 
    gff3_merge -d ${dir}/${dir}.maker.output/${dir}_master_datastore_index.log -o ${dir}/${dir}.gff3
done
```
```
java -jar /global/scratch/users/skyungyong/Software/NLR-Annotator/NLR-Annotator-v2.1b.jar -t 40  -x /global/scratch/users/skyungyong/Software/NLR-Annotator/src/mot.txt -y /global/scratch/users/skyungyong/Software/NLR-Annotator/src/store.txt -i  Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa -o NLRannotator.out -g NLRannotator.gff3
```

Finally, collect all annotations
ls
Kronos.v1.0.all.gff3  Kronos.v2.0.gff3  v1_abinitio.gff3 #braker/ginger/funannotate annotations used for EVM
less ../Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.fa | grep ">" | cut -d ">" -f 2 | sort -u > coordinates.list
