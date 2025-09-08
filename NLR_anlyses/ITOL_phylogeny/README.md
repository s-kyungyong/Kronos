# Phylogenetic tree


## Software Versions
```
hmmsearch v3.4
mafft v7.525
trimal v1.5.rev0
raxml-ng v1.2.2
```

---

### 1. Sequence alignment

**📥 Inputs**   
• `Kronos_and_known_NLRs.fa`: Kronos and reference NLR sequences (1136 sequences; 1089 Kronos sequences)  
• `Kronos_NBARC.hmm`: HMM profile constructed with reliable Kronos NLRs (1121 sequences in total; 1074 Kronos sequences)  
 
**📥 Outputs**    
• `Kronos_and_known_NLRs.hmmalign.msa.filtered.clean.fasta`: NB-ARC domain sequences of filtered NLRs  

```
#capture NB-ARC domains
hmmalign --trim -o Kronos_and_known_NLRs.hmmalign.sto Kronos_NBARC.hmm Kronos_and_known_NLRs.fa
esl-reformat fasta Kronos_and_known_NLRs.hmmalign.sto > Kronos_and_known_NLRs.hmmalign.fasta

#align NB-ARC domains
mafft --maxiterate 1000 --globalpair --thread 40 Kronos_and_known_NLRs.hmmalign.fasta > Kronos_and_known_NLRs.hmmalign.msa.fasta

#remove low coverage sequences. 
trimal -gt 0.3 -in Kronos_and_known_NLRs.hmmalign.msa.fasta -out Kronos_and_known_NLRs.hmmalign.msa.filtered.fasta
python remove_gappy_seqs.py Kronos_and_known_NLRs.hmmalign.msa.filtered.fasta Kronos_and_known_NLRs.hmmalign.msa.filtered.clean.fasta
15 gappy sequences removed
```
---
### 2. Phylogenetic tree

**📥 Inputs**    
• `Kronos_and_known_NLRs.hmmalign.msa.filtered.clean.fasta`: NB-ARC domain sequences of filtered NLRs  

**📥 Outputs**    
• `kronos_support.raxml.support`: ML tree  

```
# Step 1: ML tree search
raxml-ng \
  --msa Kronos.NLRs.reliable.hmm.msa.filtered.clean.fasta \
  --model LG+G8+F \
  --tree pars{10} \
  --search 50 \
  --threads 40 \
  --seed 12345 \
  --prefix kronos_ml

# Step 2: Bootstrap
raxml-ng --bootstrap \
  --msa Kronos.NLRs.reliable.hmm.msa.filtered.clean.fasta \
  --model LG+G8+F \
  --bs-trees 1000 \
  --threads 40 \
  --seed 12345 \
  --prefix kronos_bs

# Step 3: Compute support
raxml-ng --support \
  --tree kronos_ml.raxml.bestTree \
  --bs-trees kronos_bs.raxml.bootstraps \
  --prefix kronos_support
```
---
To visualize the phylogenetic tree, please use the files below:  

`Tree`  
• `kronos_support.raxml.support`: maximum-likelihood tree inferred from NB-ARC domains  

`Chromosome labels`  
• `00.chromosome_memberships_all`: 14 chromosomes are assigned with distinct colors   
• `00.chromosome_memberships_per_homologous_chromosome`✨: one color is assigned for each set of homologous chromosomes (e.g. 1A and 1B)  

`Branch colors`  
• `01.categorical_branch_colors_and_ring`: hvNLRs, NLR-IDs, hvNLR-IDs and known genes are colored in branches and mapped around an outer ring   
• `01.categorical_branch_colors_only`✨: hvNLRs, NLR-IDs, hvNLR-IDs and known genes are colored in bold branches  

`Known genes`  
• `02.known_NLRs`: cloned functional NLRs are indicated around an outer ring  

`Labels`  
• `03.HvNLRs_Kronos`: hvNLRs are indicated around an outer ring  
• `03.NLR-IDs_Kronos`: NLR-IDs are indicated around an outer ring  
  
`Domain architecture`   
• `04.domain_architecture`: multi-domain architectures will be visualized  

----

⚙️ **ITOL** 

Import the phylogenetic tree into ITOL and adjust parameters as below. Please note that we modified the tree in Illustrator, and therefore, main and supplemental figures are not fully reproducible with these parameters. The tree is unrooted. Make suere to root at the midpoint (see advanced option).

`Basic`
```
Mode: Circular

Rotation: 210
Arc: 350
Branch lengths: Use
Invert tree: No

Labels: Display (this shows dotted lines and labels. Dotted lines were removed in figures)
  Font: Arial
  Font style: 2px, black
  Position: Aligned
  Alignment: Left
  Rotation: On
  Shift: 0 px

Line stype: 2px, #7f7f7f, curved
Color gradient: off
Dashed line: 0.9 px, dotted lines
```

`Advanced`
```
Scaling factors: 1x, 1x
Inverted circle size: 0 px
Leaf sorting: Default
Invert sort order: No

Node IDs: Hide
Branch lengths: Hide
Bootstrap: Display
  Data source: bootstrap
  Display range: 70 to 100
  Legend: Off
  Symbol: o
  Fill color: #1b8d44
  Border color: NA
  Border width: 0 px
  Minimum size: 3px
  Maximum size: 3px
  Position on branch: 50%

Internal tree scale: Hide
Tree scale box: Display
  Label font size: 12px
  Label text: Tree scale:
  Line width: 2px
  Line color: black
  Fixed value: 0

Leaf node symbols: Hide
Internal node symbols: Hide
Auto collapse clades: NA
Delete branches: NA

Lable functions: NA
Auto assign taxonomy: NA
Midpoint root: ✨click once✨ 
```
