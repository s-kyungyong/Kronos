

import sys

chromosome = sys.argv[1]
with open('interpro.gff3', 'w') as o:
  for line in open('../Kronos.collapsed.chromosomes.masked.v1.1.fa.NLR_loci.6frameipr.gff3', 'r'):
    if line.startswith(chromosome):
        fields = line.split('\t') 
        start = int( line.split()[3] )
        end   = int( line.split()[4] )

        if start < end:
            o.write(line)
        else:
            fields[3] = str( end ) 
            fields[4] = str( start ) 
            o.write("\t".join(fields) )
