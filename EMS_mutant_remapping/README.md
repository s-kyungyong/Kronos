
SRA accessions for 1,479 experiemtns in [PRJNA258539](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=258539) were accessed, and the table in 'Capture_design.txt' was created. This file looks as below:

```
Mutant ID	SRA accession	Capture ID	Remap ID
Kronos0-1	SRX688079	KTC1	1
Kronos166	SRX688097	KTC1	1
Kronos199	SRX688114	KTC1	1
Kronos447	SRX688221	KTC1	1
Kronos499	SRX688241	KTC1	1
Kronos910	SRX688360	KTC1	1
```

The first and second columns contain the Kronos mutant line and its corresponding SRA accession. The third column (Capture ID) indicates a group of mutants which were processed for exome capture and sequenced in the same batch by [Krasileva et al., 2017](https://www.pnas.org/doi/10.1073/pnas.1619268114). We aimed to remap exome capture data for about 24 samples, by groupping three or so batches processed at similar times. This is indicated in the last column (Remap ID). Some batches did not have enough samples, and other samples needed to be brought to make the group around 24. 

Create folders
```
cat Capture_design.txt | awk '{print $4}' | grep -E '[0-9]+' | sort -u | while read int; do mkdir Remap-${int}; done
```

conda activate /global/scratch/users/skyungyong/Software/anaconda3/envs/alignment
module load sra-tools

