# Non-coding RNA annotation

## Data availability 

## Software version

```
FEELnc v0.2
Plant-LncRNA-pipeline v2
stringtie
agat
diamond
rfam v15.0
tRNAscan-SE
cmscan
```

---

## Long non-coding RNAs

### 1. Identify non-overlapping candidate transcripts  
Transcripts that are not overlapping with annotated coding genes will be selected as long non-coding RNA (lncRNA) candidates.   

**📥 Inputs**    
• `stringtie.gtf`: assembled transcripts produced during v1.0 annotations  
• `mikado.loci.gff3`: assembled transcripts during v2.0 annotations      
• `Kronos.collapsed.chromosomes.masked.v1.1.fa`: genome version v1.1  
• `Kronos.v2.1.gtf`: genome annotation v2.1  

**📥 Outputs**   
• `stringtie_candidate_lncRNA.gtf`: candidate lncRNAs from v1.0 annotations    
• `mikado_candidate_lncRNA.gtf`: candidate lncRNAs from v2.0 annotations    

```
#for mikado, extract ncRNAs. for stringtie, use all. 
awk '$3 == "ncRNA" {print $9}' mikado.loci.gff3 | cut -d ";" -f 1 | cut -d "=" -f 2 | sort -u > mikado.ncRNA.list
agat_sp_filter_feature_from_keep_list.pl --keep_list mikado.ncRNA.list --gff mikado.loci.gff3 -o mikado.loci.ncRNA.gff3
agat_convert_sp_gff2gtf.pl --gff mikado.loci.ncRNA.gff3 -o mikado.loci.ncRNA.gtf
```

```
#identify non-overlapping transcripts 
FEELnc_filter.pl -i mikado.loci.ncRNA.gtf -a Kronos.v2.1.gtf --monoex=-1 -s 200 -p 20 > mikado_candidate_lncRNA.gtf
FEELnc_filter.pl -i ../stringtie.flipped.gtf -a Kronos.v2.1.gtf --monoex=-1 -s 200 -p 20 > stringtie_candidate_lncRNA.gtf
```

```
#create gff files
cut -d ";" -f 2 stringtie_candidate_lncRNA.gtf |sed 's/ transcript_id //g' | sed 's/"//g' | sort -u > stringtie_candidate_lncRNA.txt
gffread -w stringtie_candidate_lncRNA.fa -g Kronos.collapsed.chromosomes.masked.v1.1.fa stringtie_candidate_lncRNA.gtf

cut -d ";" -f 2 mikado_candidate_lncRNA.gtf | sed 's/ transcript_id //g' | sed 's/"//g' | sort -u > mikado_candidate_lncRNA.txt
gffread -w mikado_candidate_lncRNA.fa -g Kronos.collapsed.chromosomes.masked.v1.1.fa mikado_candidate_lncRNA.gtf
```
---

### 2. Run Plant-LncRNA pipelines
We mostly followed existing pipelines available in [here](https://github.com/xuechantian/Plant-LncRNA-pipeline-v2) and [here](https://github.com/xuechantian/Plant-LncRNA-pipline).  

**📥 Inputs**  
• `stringtie_candidate_lncRNA.gtf`: candidate lncRNAs from v1.0 annotations  
• `mikado_candidate_lncRNA.gtf`: candidate lncRNAs from v2.0 annotations  
• `uniprotkb_taxonomy_id_4479.fasta`: uniprot protein sequences for poaceae 

**📥 Outputs**  
• `stringtie_ncRNA.gtf`: lncRNAs from v1.0 annotations  
• `mikado_lncRNA.gtf`: lncRNAs from v2.0 annotations  

```
#run the pipeline for both of the gtf files obtained from the previous step 
python Plant-LncRNA-pipeline-v2/Script/Feature_extraction.py -i mikado_candidate_lncRNA.fa -o mikado_PlantLncBoost_feature.csv
python Plant-LncRNA-pipeline-v2/Script/PlantLncBoost_prediction.py -i mikado_PlantLncBoost_feature.csv -m PlantLncBoost/Model/PlantLncBoost_model.cb -t 0.5 -o mikado_PlantLncBoost_prediction.csv
cpat.py -x ./Plant-LncRNA-pipeline-v2/Model/Plant_Hexamer.tsv -d ./Plant-LncRNA-pipeline-v2/Model/Plant.logit.RData -g mikado_candidate_lncRNA.fa -o mikado_CPAT_plant.output
diamond blastx --masking 0 -d uniprotkb_taxonomy_id_4479.fasta.dmnd -q mikado_candidate_lncRNA.fa -o mikado_uniprotoutput.txt
```

```
#in R

library(LncFinder)
library(seqinr)

mRNA <- seqinr::read.fasta(file ="./data/training/mRNA.fasta")
lncRNA <- seqinr::read.fasta(file ="./data/training/lncRNA.fasta")
frequencies <- make_frequencies(cds.seq = mRNA, lncRNA.seq = lncRNA, SS.features = FALSE, cds.format = "DNA", lnc.format = "DNA", check.cds = TRUE, ignore.illegal = TRUE)
plant = readRDS("./Model/Plant_model.rda")
Seqs <- seqinr::read.fasta(file ="mikado_candidate_lncRNA.fa")
Plant_results <- LncFinder::lnc_finder(Seqs, SS.features = FALSE, format = "DNA", frequencies.file = frequencies, svm.model = plant, parallel.cores = 2)
write.table(Plant_results, file ="mikado_plant-lncFinder.txt", sep ="\t",row.names =TRUE, col.names =TRUE,quote =FALSE)
```

```
#combine all results
Rscript prediction_insersection.sh mikado_candidate_lncRNA.txt mikado_PlantLncBoost_prediction.csv mikado_CPAT_plant.output plant-lncFinder.txt mikado_uniprotoutput.txt
grep -Fwf mikado_Final_lncRNA_results.txt mikado_candidate_lncRNA.gtf > mikado_lncRNA.gtf
```
---

### 3. Classification  

**📥 Inputs**  
• `stringtie_ncRNA.gtf`: lncRNAs from v1.0 annotations    
• `mikado_lncRNA.gtf`: lncRNAs from v2.0 annotations    
• `Kronos.v2.1.gtf`: genome annotation v2.1  

**📥 Outputs**  
• `lncRNA_classes.txt`: lncRNA classes  

```
#merge transcripts and classify
stringtie --merge mikado_final_lncRNA.gtf stringtie_final_lncRNA.gtf > mikado_stringtie_merged_final_lncRNA.gtf
FEELnc_classifier.pl -i mikado_stringtie_merged_final_lncRNA.gtf -a Glycine_max_longest.gtf > lncRNA_classes.txt
```
---


## Small non-coding RNAs

### 1. Rfam search  
**📥 Inputs**    
• `Kronos.collapsed.chromosomes.masked.v1.1.fa`: genome version v1.1    

**📥 Outputs**   
• `*.Rfam.tblout`: Rfam outputs  
 
```
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz #v15.0 
gunzip Rfam.cm.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
cmpress Rfam.cm
```

```
#for each chromosome
while read -r chromosome zvalue; do
  cmscan -Z ${zvalue} --cut_ga --rfam --nohmmonly --tblout ${chromosome}.Rfam.tblout --fmt 2 --cpu 56 --clanin Rfam.clanin Rfam.cm Kronos.v1.1.${chromosome}.fa
done < zvalue.list
```
---

### 2. tRNA search

**📥 Inputs**  
• `Kronos.collapsed.chromosomes.masked.v1.1.fa`: genome version v1.1  


**📥 Outputs**  
• `tRNAscan-SE.out`: tRNA scan outputs

```
tRNAscan-SE -E -o tRNAscan-SE.out -f tRNAscan-SE.ss -s tRNAscan-SE.iso -m tRNAscan-SE.stats -c tRNAscan-SE.conf Kronos.collapsed.chromosomes.v1.1.fa
```
---


1.6.1
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/fasta/triticum_aestivum/ncrna/Triticum_aestivum.IWGSC.ncrna.fa.gz
unitas.pl -species x -refseq Triticum_aestivum.IWGSC.ncrna.fa -input all.fq -threads 48


