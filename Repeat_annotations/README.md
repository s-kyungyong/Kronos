# Repeat Annotation

### Repeat Identification with HiTE
Repeative elements annotated in the Kronos reference genome v1.0 and v1.1 are predicted by (HiTE)[https://github.com/CSU-KangHu/HiTE] v3.0.0.
```
singularity run HiTE.sif python main.py --genome Kronos.collapsed.chromosomes.fa \
     --thread 56 --outdir HiTE --recover 1 --annotate 1 \
     --plant 1 --classified 1 --domain 1 --recover 1 --debug 1
```

### Repeat Annotations with EDTA
The second version of repeat elements is produced by EDTA v2.2.2. Some input files need to be prepared. 

TREP 
```
