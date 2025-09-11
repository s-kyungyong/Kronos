################################
####  PHAS LOCI BED COORD.  ####
################################
rule dcl_phasirna_subset:
    """
    Extract lines of PHAS loci significant at PhaseScore >= 40, subset loci of 21-nt and 24-nt phasiRNAs and classify groups of
    pre-/mid-/post-meiotic phasiRNAs
    """
    input:
        phas21= rules.dcl_phasirna_annot.output.phas21,
        phas24= rules.dcl_phasirna_annot.output.phas24,
    output:
        phas21PRE=    "phasirna/{specie}/5report/21phasPRE.txt",
        phas24PRE=    "phasirna/{specie}/5report/24phasPRE.txt",
        phas21NOT=    "phasirna/{specie}/5report/21phasNOT.txt",
        phas24NOT=    "phasirna/{specie}/5report/24phasNOT.txt",
        list21PRE=    temp("phasirna/{specie}/5report/list/21phasPRE.txt"),
        list24PRE=    temp("phasirna/{specie}/5report/list/24phasPRE.txt"),
        list21NOT=    temp("phasirna/{specie}/5report/list/21phasNOT.txt"),
        list24NOT=    temp("phasirna/{specie}/5report/list/24phasNOT.txt"),
    message: """ Extract, subset and classify phasiRNAs """
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/phasirna.subset.R"



rule dcl_phasirna_phas_coordinates:
    """
    Extract the coordinates of phasiRNA loci per categories
    """
    input:
        phas21PRE=  rules.dcl_phasirna_subset.output.list21PRE,
        phas24PRE=  rules.dcl_phasirna_subset.output.list24PRE,
        phas21NOT=  rules.dcl_phasirna_subset.output.list21NOT,
        phas24NOT=  rules.dcl_phasirna_subset.output.list24NOT,
    output:
        phas21PRE=    "phasirna/{specie}/5motif/1annot/21phasPRE.txt",
        phas24PRE=    "phasirna/{specie}/5motif/1annot/24phasPRE.txt",
        phas21NOT=    "phasirna/{specie}/5motif/1annot/21phasNOT.txt",
        phas24NOT=    "phasirna/{specie}/5motif/1annot/24phasNOT.txt",
    params:
        "phasirna/{specie}/2call/ShortStack_All.gff3"
    conda:
        "../envs/annotate.yaml"
    shell:
        """
        rg -w -f {input.phas21PRE} {params} | sed 's/ID=//g' | sed 's/;/\t/g' | awk '{{print $1 "\t" $4 "\t" $5 "\t" $9}}' | sort-bed --max-mem 20G - > {output.phas21PRE}
        rg -w -f {input.phas24PRE} {params} | sed 's/ID=//g' | sed 's/;/\t/g' | awk '{{print $1 "\t" $4 "\t" $5 "\t" $9}}' | sort-bed --max-mem 20G - > {output.phas24PRE}
        rg -w -f {input.phas21NOT} {params} | sed 's/ID=//g' | sed 's/;/\t/g' | awk '{{print $1 "\t" $4 "\t" $5 "\t" $9}}' | sort-bed --max-mem 20G - > {output.phas21NOT}
        rg -w -f {input.phas24NOT} {params} | sed 's/ID=//g' | sed 's/;/\t/g' | awk '{{print $1 "\t" $4 "\t" $5 "\t" $9}}' | sort-bed --max-mem 20G - > {output.phas24NOT}
        """

#############################################
####  EXTRACT PHAS LOCI BED COORDINATES  ####
#############################################
rule dcl_phasirna_phas_loci_regions_extend:
    """
    Extend the regions of putative PHAS loci, in BED files
    """
    input:
        phas21PRE=  rules.dcl_phasirna_phas_coordinates.output.phas21PRE,
        phas24PRE=  rules.dcl_phasirna_phas_coordinates.output.phas24PRE,
        phas21NOT=  rules.dcl_phasirna_phas_coordinates.output.phas21NOT,
        phas24NOT=  rules.dcl_phasirna_phas_coordinates.output.phas24NOT,
    output:
        phas21PRE=  "phasirna/{specie}/5motif/2coord/21phasPRE.txt",
        phas24PRE=  "phasirna/{specie}/5motif/2coord/24phasPRE.txt",
        phas21NOT=  "phasirna/{specie}/5motif/2coord/21phasNOT.txt",
        phas24NOT=  "phasirna/{specie}/5motif/2coord/24phasNOT.txt",
    params:
        CHR
    conda:
        "../envs/annotate.yaml"
    shell:
        """
        bedtools slop -i {input.phas21PRE} -g {params} -b 1000 > {output.phas21PRE}
        bedtools slop -i {input.phas24PRE} -g {params} -b 1000 > {output.phas24PRE}
        bedtools slop -i {input.phas21NOT} -g {params} -b 1000 > {output.phas21NOT}
        bedtools slop -i {input.phas24NOT} -g {params} -b 1000 > {output.phas24NOT}
        """


rule dcl_phasirna_phas_loci_regions_extend_seq:
    """
    Extract putative PHAS precursor transcript sequences from the BED file
    """
    input:
        phas21PRE=  rules.dcl_phasirna_phas_loci_regions_extend.output.phas21PRE,
        phas24PRE=  rules.dcl_phasirna_phas_loci_regions_extend.output.phas24PRE,
        phas21NOT=  rules.dcl_phasirna_phas_loci_regions_extend.output.phas21NOT,
        phas24NOT=  rules.dcl_phasirna_phas_loci_regions_extend.output.phas24NOT,
    output:
        phas21PRE=      "phasirna/{specie}/5motif/3fasta/21phasPRE.fasta",
        phas24PRE=      "phasirna/{specie}/5motif/3fasta/24phasPRE.fasta",
        phas21NOT=      "phasirna/{specie}/5motif/3fasta/21phasNOT.fasta",
        phas24NOT=      "phasirna/{specie}/5motif/3fasta/24phasNOT.fasta",
    params:
        FASTA
    conda:
        "../envs/annotate.yaml"
    shell:
        """
        bedtools getfasta -nameOnly -fi {params} -bed {input.phas21PRE}    | sed 's/Cluster_/{wildcards.specie}_PHAS_/g' > {output.phas21PRE}
        bedtools getfasta -nameOnly -fi {params} -bed {input.phas24PRE}    | sed 's/Cluster_/{wildcards.specie}_PHAS_/g' > {output.phas24PRE}
        bedtools getfasta -nameOnly -fi {params} -bed {input.phas21NOT} | sed 's/Cluster_/{wildcards.specie}_PHAS_/g' > {output.phas21NOT}
        bedtools getfasta -nameOnly -fi {params} -bed {input.phas24NOT} | sed 's/Cluster_/{wildcards.specie}_PHAS_/g' > {output.phas24NOT}
        """



####################################################################
####  IDENTIFY AND PRINT CONSERVED MOTIFS - SPECIES-BY-SPECIES  ####
####################################################################
rule dcl_phasirna_phas_motifs:
    """
    Identify and print the logo of conserve motifs enriched in putative PHAS precursor transcripts
    """
    input:
        phas21PRE=  rules.dcl_phasirna_phas_loci_regions_extend_seq.output.phas21PRE,
        phas24PRE=  rules.dcl_phasirna_phas_loci_regions_extend_seq.output.phas24PRE,
        phas21NOT=  rules.dcl_phasirna_phas_loci_regions_extend_seq.output.phas21NOT,
        phas24NOT=  rules.dcl_phasirna_phas_loci_regions_extend_seq.output.phas24NOT,
    output:
        phas21PRE=  "phasirna/{specie}/5motif/4motif/21phasPRE/meme.txt",
        phas24PRE=  "phasirna/{specie}/5motif/4motif/24phasPRE/meme.txt",
        phas21NOT=  "phasirna/{specie}/5motif/4motif/21phasNOT/meme.txt",
        phas24NOT=  "phasirna/{specie}/5motif/4motif/24phasNOT/meme.txt",
    conda:
        "../envs/visual.yaml"
    shell:
        """
        meme -oc phasirna/{wildcards.specie}/5motif/4motif/21phasPRE -objfun classic -minw 21 -maxw 22 -nmotifs 3 -revcomp -mod zoops -dna {input.phas21PRE}
        meme -oc phasirna/{wildcards.specie}/5motif/4motif/24phasPRE -objfun classic -minw 21 -maxw 22 -nmotifs 3 -revcomp -mod zoops -dna {input.phas24PRE}
        meme -oc phasirna/{wildcards.specie}/5motif/4motif/21phasNOT -objfun classic -minw 21 -maxw 22 -nmotifs 3 -revcomp -mod zoops -dna {input.phas21NOT}
        meme -oc phasirna/{wildcards.specie}/5motif/4motif/24phasNOT -objfun classic -minw 21 -maxw 22 -nmotifs 3 -revcomp -mod zoops -dna {input.phas24NOT}
        """



