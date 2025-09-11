###########################
# Drop count file columns #
###########################
rule dcl_phasirna_dicer:
    """
    Make a list of loci with valid DicerCall
    """
    input:
        "phasirna/{specie}/2call/Results.txt"
    output:
        temp("phasirna/{specie}/3clean/ListScoredLoci.txt")
    message: """ Remove loci where DicerCall is 'N' and get the list of scored loci """
    conda:
        "../envs/format.yaml"
    shell:
        """
        rg --invert-match '\tN\t' {input} | awk '{{if($5 >= 0.5)print}}' |  awk '{{print $2}}' > {output}
        """



rule dcl_phasirna_count:
    """
    Select lines with a valid DicerCall
    """
    input:
        list=   rules.dcl_phasirna_dicer.output,
        counts= "phasirna/{specie}/2call/Counts.txt",
    output:
        counts=  "phasirna/{specie}/3clean/Counts.txt",
    message: """ Select lines with a valid DicerCall """
    conda:
        "../envs/format.yaml"
    shell:
        """
        awk 'NR == 1; NR > 1 {{print $0 | "rg -w -f {input.list}"}}' {input.counts} | cut -f 2,4-57 | sed 's/dcl5_RS_//g' > {output}
        """


#################################################################
# Transform the count matrix in count per million mapped (rpmm) #
#################################################################
rule dcl_phasirna_rpm:
    """
    Transform raw count in rpmm matrix (from code example: https://www.reneshbedre.com/blog/expression_units.html)
    """
    input:
        rules.dcl_phasirna_count.output.counts
    output:
        "phasirna/{specie}/3clean/RpmCounts.txt"
    conda:
        "../envs/bioinfokit.yaml"
    script:
        "../scripts/phasirna.bioinfokit.py"



###########################################################################
# Bind result and count files and remove loci that are not Dicer products #
###########################################################################
rule dcl_phasirna_cbind:
    """
    Bind result and count files and remove loci that are not Dicer products
    """
    input:
        result= "phasirna/{specie}/2call/Results.txt",
        raw=    rules.dcl_phasirna_count.output.counts,
        rpm=    rules.dcl_phasirna_rpm.output,
    output:
        srna= "phasirna/{specie}/3clean/srna.txt"
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/phasirna.cbind.R"


