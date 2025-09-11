rule dcl_phasirna_call:
    """
    Annotate phasiRNA loci using Shortstack v3.8.5
    """
    output:
        directory("phasirna/{specie}/2call/"),
    params:
        FASTA
    message: """ Annotate known and novel microRNA Shortstack v3.8.5 """
    shell:
        """
        /Users/sbelanger/software/ShortStack-3.8.5/ShortStack \
            --bowtie_cores 20 \
            --sort_mem 20G \
            --bowtie_m 50 \
            --dicermin 20 \
            --dicermax 24 \
            --mmap u \
            --mincov 1.0rpm \
            --pad 75 \
            --outdir {output} \
            --bamfile /Users/sbelanger/analysis/kronos/phasirna/{wildcards.specie}/1align/merged_alignments.bam \
            --genomefile {params}
        """
