#download the databases: rfam v15.0
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
cmpress Rfam.cm

#preparezvalues for each chromosome
cat zvalue.list
1A 1200.879162
1B 1417.658372
2A 1591.626378
2B 1657.049066
3A 1518.249256
3B 1728.269574
4A 1535.707034
4B 1399.361912
5A 1440.551718
5B 1462.269252
6A 1248.596746
6B 1467.16489
7A 1506.943932
7B 1532.01879
Un 421.039088

#for each chromosome
while read -r chromosome zvalue; do
  cmscan -Z ${zvalue} --cut_ga --rfam --nohmmonly --tblout ${chromosome}.Rfam.tblout --fmt 2 --cpu 56 --clanin Rfam.clanin Rfam.cm Kronos.v1.1.${chromosome}.fa
done < seqLengths.list

./tRNAscan-SE_installed/bin/tRNAscan-SE -E -o tRNAscan-SE.out -f tRNAscan-SE.ss -s tRNAscan-SE.iso -m tRNAscan-SE.stats -c ./tRNAscan-SE_installed/bin/tRNAscan-SE.conf ../Final/Kronos.collapsed.chromosomes.v1.1.fa

1.6.1
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/fasta/triticum_aestivum/ncrna/Triticum_aestivum.IWGSC.ncrna.fa.gz
unitas.pl -species x -refseq Triticum_aestivum.IWGSC.ncrna.fa -input all.fq -threads 48

pwd
/global/scratch/users/skyungyong/Kronos/5.Annotations/lncRNA

ls
mikado.loci.gff3  stringtie.flipped.gtf
awk '$3 == "ncRNA" {print $9}' mikado.loci.gff3 | cut -d ";" -f 1 | cut -d "=" -f 2 | sort -u > mikado.ncRNA.list
agat_sp_filter_feature_from_keep_list.pl --keep_list mikado.ncRNA.list --gff mikado.loci.gff3 -o mikado.loci.ncRNA.gff3
agat_convert_sp_gff2gtf.pl --gff mikado.loci.ncRNA.gff3 -o mikado.loci.ncRNA.gtf

stringtie --merge mikado.loci.nc.gtf stringtie.flipped.gtf > lncRNA/merged_transcripts.gtf

FEELnc_filter.pl -i ../mikado.loci.ncRNA.gtf -a /global/scratch/projects/vector_kvklab/KS-Kronos_Final_datasets/02.Annotations/00.Proteins/v2.1/Kronos.v2.1.gtf --monoex=-1 -s 200 -p 20 > mikado_candidate_lncRNA.gtf
FEELnc_filter.pl -i ../stringtie.flipped.gtf -a /global/scratch/projects/vector_kvklab/KS-Kronos_Final_datasets/02.Annotations/00.Proteins/v2.1/Kronos.v2.1.gtf --monoex=-1 -s 200 -p 20 > stringtie_candidate_lncRNA.gtf

cut -d ";" -f 2 stringtie_candidate_lncRNA.gtf |sed 's/ transcript_id //g' | sed 's/"//g' | sort -u > stringtie_candidate_lncRNA.txt
cut -d ";" -f 2 mikado_candidate_lncRNA.gtf |sed 's/ transcript_id //g' | sed 's/"//g' | sort -u > mikado_candidate_lncRNA.txt


gffread -w mikado_candidate_lncRNA.fa -g  /global/scratch/projects/vector_kvklab/KS-Kronos_Final_datasets/00.Genomes/99.Genomes/Kronos.collapsed.chromosomes.masked.v1.1.fa mikado_candidate_lncRNA.gtf
gffread -w stringtie_candidate_lncRNA.fa -g /global/scratch/projects/vector_kvklab/KS-Kronos_Final_datasets/00.Genomes/99.Genomes/Kronos.collapsed.chromosomes.masked.v1.1.fa stringtie_candidate_lncRNA.gtf

python Plant-LncRNA-pipeline-v2/Script/Feature_extraction.py -i mikado_candidate_lncRNA.fa -o mikado_PlantLncBoost_feature.csv
python Plant-LncRNA-pipeline-v2/Script/Feature_extraction.py -i stringtie_candidate_lncRNA.fa -o stringtie_PlantLncBoost_feature.csv

python Plant-LncRNA-pipeline-v2/Script/PlantLncBoost_prediction.py -i mikado_PlantLncBoost_feature.csv -m PlantLncBoost/Model/PlantLncBoost_model.cb -t 0.5 -o mikado_PlantLncBoost_prediction.csv

cpat.py -x ./Plant-LncRNA-pipeline-v2/Model/Plant_Hexamer.tsv -d ./Plant-LncRNA-pipeline-v2/Model/Plant.logit.RData -g mikado_candidate_lncRNA.fa -o mikado_CPAT_plant.output
cpat.py -x ./Plant-LncRNA-pipeline-v2/Model/Plant_Hexamer.tsv -d ./Plant-LncRNA-pipeline-v2/Model/Plant.logit.RData -g stringtie_candidate_lncRNA.fa -o stringtie_CPAT_plant.output

 diamond blastx --masking 0 -d /global/scratch/users/skyungyong/Kronos/10.Genome_comparisons/Annotation_filtering_evaluation/uniprotkb_taxonomy_id_4479_2025_02_18.fasta.dmnd -q mikado_candidate_lncRNA.fa -o mikado_uniprotoutput.txt
stringtie --merge mikado_final_lncRNA.gtf stringtie_final_lncRNA.gtf > mikado_stringtie_merged_final_lncRNA.gtf
