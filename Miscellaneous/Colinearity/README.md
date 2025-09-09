# Colinearity analyses 

## Software version
```
mcscanx v1.0.0
jcvi v1.5.7
diamond v2.1.9
```

------
### 1. Colinearity analyses with MCScanX

**📥 Inputs**  
• `AET_High_confidence_gene_protein.longest.fasta`: longest protein sequences per gene for Ae. tauschii [source](http://aegilops.wheat.ucdavis.edu/ATGSP/annotation/)  
• `AET_High_confidence_gene.gff3`: gff file for Ae. tauschii  
• `Kronos.v2.1.pep.longest.fa`: longest protein sequences per gene for Kronos    
• `Kronos.v2.1.gff`: gff file for Kronos  

**📥 Outputs**    
• `Kronos_vs_AET.collinearity`: colinearity prediction from mcscanx  

⚙️ **Similarity search**  
```
diamond makedb --in AET_High_confidence_gene_protein.longest.fasta --db AET_High_confidence_gene_protein.longest
diamond blastp -q Kronos.v2.1.pep.longest.fa -d AET_High_confidence_gene_protein.longest \
               -o Kronos_vs_AET.blast -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
               --masking 0 --max-target-seqs 5 --evalue 1e-10 --threads 56
```

⚙️ **Reformat gff**  
```
python get_gff.py AET_High_confidence_gene.gff3 AET
python get_gff.py Kronos.v2.1.gff KR
cat KR.bed AET.bed > Kronos_vs_AET.bed
cat KR.gff AET.gff > Kronos_vs_AET.gff
```

⚙️ **MCScanX**  
```
MCScanX -a Kronos_vs_AET
```

⚙️ **Visualize a dot plot**  

```
grep -v ">" Kronos_vs_AET.collinearity | awk -F "\t" '{print $3 "\t" $2}' | grep 'KRN' > Kronos_vs_AET.anchors
python -m jcvi.graphics.dotplot --qbed Kronos.bed --sbed AET.bed --dpi 500 Kronos_vs_AET.anchors
```
