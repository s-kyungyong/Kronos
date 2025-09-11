



#find the best hit and print CDS of the match, as well as query sequences
exonerate --model protein2genome --bestn 1 --showtargetgff yes \
    --ryo ">cds|t:%ti|q:%qi|%tab-%tae|strand:%tS|pid:%pi|score:%s\n%tcs\n>prot|q:%qi\n%qs\n" \
    TrturKRN2B02G077030.1.fa Svevo.2B.fa > 1.out



wget https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/linux64.table2asn.gz
agat_sp_fix_cds_phases.pl 
gt gff3validator your_annotations.gff3
./linux64.table2asn -M n -J -c w -euk -t template.sbt  -gaps-min 10 -l paired-ends -j "[organism=Triticum turgidum]" -i Kronos.v1.1.folded.masked.NoUn.fa -f Kronos.v2.1.clean.sorted.fixed.gff3 -o Kronos.v1.1.folded.masked.NoUn.sqn -Z -V b

