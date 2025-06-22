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
