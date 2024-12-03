import os
import pysam
from collections import defaultdict

import os
import pysam
from collections import defaultdict

# Get list of BAM files
bam_files = [b for b in os.listdir() if b.endswith('.bam') and 'primary' not in b]

for bam_file in bam_files:
    # Dictionary to store counts of primary alignments per read
    primary_alignments_count = defaultdict(int)

    # First pass: Count primary alignments per read
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if not read.is_secondary and not read.is_supplementary:
                primary_alignments_count[read.query_name] += 1

    # Open output BAM file for writing
    output_file = bam_file.replace('.bam', '.UniqPrimary.bam')
    ct = 0  # Counter for written reads
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        header = bam.header  # Retrieve header
        with pysam.AlignmentFile(output_file, "wb", header=header) as out_bam:
            for read in bam:
                # Write reads with exactly one primary alignment
                if (
                    not read.is_secondary
                    and not read.is_supplementary
                    and primary_alignments_count[read.query_name] == 1
                ):
                    out_bam.write(read)
                    ct += 1

    # Summarize counts for this BAM file
    print(f"Processed {bam_file}. {ct} alignments were written to {output_file}.")
