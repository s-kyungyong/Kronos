
For IsoSeq: We used data from two BioProject accessions: PRJEB15048 and PRJNA427246

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


/global/scratch/users/skyungyong/Software/smrtlink/smrtcmds/bin/bax2bam  --subread -o ERR1597909.subreads.bam *.bax.h5

