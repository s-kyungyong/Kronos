# Phylogenetic tree

To visualize the phylogenetic tree, please use the files below
• `kronos_support.raxml.support`: maximum-likelihood tree inferred from NB-ARC domains


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

Line stype: 2px, curved
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
Midpoint root: click once 
```
