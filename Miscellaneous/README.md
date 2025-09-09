



#find the best hit and print CDS of the match, as well as query sequences
exonerate --model protein2genome --bestn 1 --showtargetgff yes \
    --ryo ">cds|t:%ti|q:%qi|%tab-%tae|strand:%tS|pid:%pi|score:%s\n%tcs\n>prot|q:%qi\n%qs\n" \
    TrturKRN2B02G077030.1.fa Svevo.2B.fa > 1.out
