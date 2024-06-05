
For IsoSeq: 

We downloaded available long-read RNA-seq datasets with focus on Sequel I and II. If HiFi or CCS reads are uploaded into SRA, we directly used them. If subread bam files were uploaded, we processed them. 




We used data from two BioProject accessions: PRJEB15048 and PRJNA427246

wget https://downloads.pacbcloud.com/public/software/installers/smrtlink_7.0.1.66975.zip
unzip smrtlink_7.0.1.66975.zip 
bash smrtlink_7.0.1.66975.run --rootdir /global/scratch/users/skyungyong/Software/smrtlink
------
#### **A. Data download**
------


```
grep 'ERR' SRA.list | while read -r accession; do
  cd $accession
  wget -c -r -np -nH --cut-dirs=5 -R "index.html*" ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR159/"${accession}"/
  cd ..
```


/global/scratch/users/skyungyong/Software/smrtlink/smrtcmds/bin/bax2bam  --subread -o ERR1597909 *.bax.h5
/global/scratch/users/skyungyong/Software/smrtlink/smrtcmds/bin/ccs -j 56 --min-pass 1 ERR1597909.subreads.bam ERR1597909.ccs.bam
/global/scratch/users/skyungyong/Software/smrtlink/smrtcmds/bin/lima --dump-clips --peek-guess --isoseq -j 56 ERR1597909.ccs.bam barcodes.fa ERR1597909.ccs.lima.bam
/global/scratch/users/skyungyong/Software/smrtlink/smrtcmds/bin/isoseq3 refine --require-polya ERR1597909.ccs.lima.Clontech_5p--NEB_Clontech_3p.bam barcodes.fa ERR1597909.flnc.bam
/global/scratch/users/skyungyong/Software/smrtlink/smrtcmds/bin/isoseq3 cluster ERR1597909.flnc.bam unpolished.bam
/global/scratch/users/skyungyong/Software/smrtlink/smrtcmds/bin/isoseq3 polish unpolished.bam ERR1597909.subreads.bam polished.bam
