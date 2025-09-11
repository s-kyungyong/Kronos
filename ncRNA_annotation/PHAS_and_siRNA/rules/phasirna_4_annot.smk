################################
####   SORT PHASIRNA LOCI   ####
################################
rule dcl_phasirna_annot:
    """
    Extract lines of PHAS loci significant at PhaseScore >= 40, subset 21-nt and 24-nt phasiRNAs
    and classify groups of pre-/mid-/post-meiotic phasiRNAs
    """
    input:
        all= rules.dcl_phasirna_cbind.output.srna,
    output:
        srna=           "phasirna/{specie}/4annot/summary.txt",
        phas=           "phasirna/{specie}/4annot/phasirna.txt",
        phas21=         "phasirna/{specie}/4annot/21phas.txt",
        phas24=         "phasirna/{specie}/4annot/24phas.txt",
        ambiguous=      "phasirna/{specie}/4annot/ambiguous.txt",
        ambiguous21=    "phasirna/{specie}/4annot/21ambiguous.txt",
        ambiguous24=    "phasirna/{specie}/4annot/24ambiguous.txt",
        sirna=          "phasirna/{specie}/4annot/sirna.txt",
    message: """ Extract, subset and classify phasiRNAs """
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/phasirna.annot.R"


rule dcl_phasirna_annot_count:
    """
    Extract raw or normalized count matrix for data visualization
    """
    input:
        srna=        rules.dcl_phasirna_annot.output.srna,
        phas21=      rules.dcl_phasirna_annot.output.phas21,
        phas24=      rules.dcl_phasirna_annot.output.phas24,
        ambiguous21= rules.dcl_phasirna_annot.output.ambiguous21,
        ambiguous24= rules.dcl_phasirna_annot.output.ambiguous24,
        sirna=       rules.dcl_phasirna_annot.output.sirna,
    output:
        srna=       "phasirna/{specie}/4count/all.txt",
        phas21=     "phasirna/{specie}/4count/21phas.txt",
        phas24=     "phasirna/{specie}/4count/24phas.txt",
        ambiguous21="phasirna/{specie}/4count/21ambiguous.txt",
        ambiguous24="phasirna/{specie}/4count/24ambiguous.txt",
        sirna=      "phasirna/{specie}/4count/sirna.txt",
    shell:
        """
        cut -f 2,76-129 {input.srna} | sed 's/.y//g' > {output.srna}
        cut -f 2,76-129 {input.phas21} | sed 's/.y//g' > {output.phas21}
        cut -f 2,76-129 {input.phas24} | sed 's/.y//g' > {output.phas24}
        cut -f 2,76-129 {input.sirna} | sed 's/.y//g' > {output.sirna}
        cut -f 2,76-129 {input.ambiguous21} | sed 's/.y//g' > {output.ambiguous21}
        cut -f 2,76-129 {input.ambiguous24} | sed 's/.y//g' > {output.ambiguous24}
        """



rule dcl_phasirna_annot_coord:
    """
    Extract coordinates of annotated loci
    """
    input:
        srna=       rules.dcl_phasirna_annot.output.srna,
        phas21=     rules.dcl_phasirna_annot.output.phas21,
        phas24=     rules.dcl_phasirna_annot.output.phas24,
        ambiguous21=rules.dcl_phasirna_annot.output.ambiguous21,
        ambiguous24=rules.dcl_phasirna_annot.output.ambiguous24,
        sirna=      rules.dcl_phasirna_annot.output.sirna,
    output:
        srna=        "phasirna/{specie}/4coord/all.txt",
        phas21=      "phasirna/{specie}/4coord/21phas.txt",
        phas24=      "phasirna/{specie}/4coord/24phas.txt",
        ambiguous21= "phasirna/{specie}/4coord/21ambiguous.txt",
        ambiguous24= "phasirna/{specie}/4coord/24ambiguous.txt",
        sirna=       "phasirna/{specie}/4coord/sirna.txt",
    shell:
        """
        cut -f 1-21 {input.srna} > {output.srna}
        cut -f 1-21 {input.phas21} > {output.phas21}
        cut -f 1-21 {input.phas24} > {output.phas24}
        cut -f 1-21 {input.sirna} > {output.sirna}
        cut -f 1-21 {input.ambiguous21} > {output.ambiguous21}
        cut -f 1-21 {input.ambiguous24} > {output.ambiguous24}
        """


#########################################
####  MDS PLOT - DATASET VALIDATION  ####
#########################################
rule dcl_phasirna_mds_plot:
    input:
        matrix= rules.dcl_phasirna_annot_count.output.srna,
        design= "data/{specie}_design.txt",
    output:
        mds=    "phasirna/{specie}/4plots/mds.pdf",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/phasirna.mds.R"



###############################################
####  HEATMAP - SIRNA TEMPORAL EXPRESSION  ####
###############################################
rule dcl_phasirna_heat_plot:
    input:
        phas21=     rules.dcl_phasirna_annot_count.output.phas21,
        phas24=     rules.dcl_phasirna_annot_count.output.phas24,
        ambiguous21=rules.dcl_phasirna_annot_count.output.ambiguous21,
        ambiguous24=rules.dcl_phasirna_annot_count.output.ambiguous24,
        sirna=      rules.dcl_phasirna_annot_count.output.sirna,
    output:
        phas21=      "phasirna/{specie}/4plots/heat_phas21.pdf",
        phas24=      "phasirna/{specie}/4plots/heat_phas24.pdf",
        ambiguous21= "phasirna/{specie}/4plots/heat_ambiguous21.pdf",
        ambiguous24= "phasirna/{specie}/4plots/heat_ambiguous24.pdf",
        sirna=       "phasirna/{specie}/4plots/heat_sirna.pdf",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/phasirna.heat.R"
