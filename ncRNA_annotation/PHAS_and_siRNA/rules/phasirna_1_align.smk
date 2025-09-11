rule srna_align:
    """
    Align sRNA reads using Shortstack v4.0.2
    """
    output:
        dir= directory("phasirna/{specie}/1align/"),
        log= "srna/{specie}.log",
    params:
        FASTA
    message: """ Align sRNA reads using Shortstack v4.0.2 """
    shell:
        """
        ShortStack --genomefile {params} --readfile /Users/sbelanger/data/{wildcards.specie}/srna/clean/fastq/dcl5_RS_* --outdir {output.dir} --threads 20 --mmap u --mincov 1 --align_only --dicermin 20 --dicermax 24 --pad 75 > {output.log}
        """
