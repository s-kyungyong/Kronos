from operator import itemgetter
from Bio import SeqIO
import logging
from collections import defaultdict
import sys

#Reference wheat genome sequence IDs
ref_order = [ f'{i}{letter}' for i in range(1,8) for letter in ['A', 'B', 'D']]

# Predefiend Fragment size (used to break scaffolds for minimap)
frag_size = 100000000
overlap   = frag_size//4

# long enough alignment size cutoff
cutoff = 2500000

def process_primary_alignments(line):
  # Process primary alignments only
  items = line.split()
  return items[12] == 'tp:A:P', int(items[3]) - int(items[2]) + 1, items[5] if 'NC' in items[5] else items[5].rsplit("_", 1)[0]

def process_long_alignments(alignment_type, alignment_length, ref_id, cutoff):
  return (alignment_type == 'Plasmid') or (alignment_length > cutoff and ref_id != 'Un')

def find_overlap(ranges, new_positions):
  # Check if the new positions have overhaps with the existing coordinates

  new_start, new_end = new_positions
  no_overlap = True

  #Iterate over the reversed list
  reversed = [ ranges[ len(ranges) - (i + 1) ] for i in range(0, len(ranges)) ]
  for idx in range(len(ranges) - 1, -1, -1):
    start, end = ranges[idx]

    # The new coordinates are within the existing coordinates
    if start <= new_start <= end and start <= new_end <= end:
      ranges[idx] = (start, end)
      no_overlap = False
      break

    # Extend to the left
    elif new_start < start and start <= new_end <= end:
      ranges[idx] = (new_start, end)
      no_overlap = False
      break

    # Extend to the right
    elif start <= new_start <= end and end < new_end:
      ranges[idx] = (start, new_end)
      no_overlap = False
      break

    # new coordinates contain previous ones
    elif new_start <= start <= new_end and new_start <= end <= new_end:
      ranges[idx] = (new_start, new_end)
      no_overlap = False
      break

  # Iterations are finished, no overlaps found
  if no_overlap:
    ranges.append( (new_start, new_end) )

  return sorted(ranges)

def merge_overlapping_ranges(ranges):

    if not ranges:
      return []

    # Sort the ranges based on the start value and initialize
    ranges.sort(key=lambda x: x[0])
    merged_ranges = [ranges[0]]

    # Iterate over the remaining ranges
    for current_range in ranges[1:]:
      prev_start, prev_end = merged_ranges[-1]

      if current_range[0] <= prev_end:
        # Ranges overlap, update the end value of the current range
        merged_ranges[-1] = (prev_start, max(prev_end, current_ranges[1]))

      else:
        # Ranges don't overlap, add the current range to the merged list and update the current range
        merged_ranges.append(current_range)

    return merged_ranges


def identify_plasmids(minimap_output, scaf_length):

  plasmid_match = process_alignment(minimap_output, 'Plasmid')
  plasmids = []
  plasmids_type = {'NC_002762.1': 'Chloroplast', 'NC_036024.1': 'Mitochondria'}

  for scaf in plasmid_match:
    q, h = scaf.split("+")
    final_ranges = merge_overlapping_ranges( plasmid_match[scaf] )
    cov = sum([end - start + 1 for start, end in final_ranges ])
    percov = round( cov / scaf_length[q] , 4)

    plasmid_type = plasmids_type.get(h)
    if plasmid_type and percov >= 0.7:
      plasmids.append((q, plasmid_type, percov))
      logging.info(f'{q} assigned to {plasmid_type}: {round(percov * 100, 3)}% coverage matches')

  chloroplasts = [q for q, pl_type, _ in plasmids if pl_type == 'Chloroplast']
  mitochondria = [q for q, pl_type, _ in plasmids if pl_type == 'Mitochondria']

  return chloroplasts, mitochondria

def process_alignment_line(line, alignment_type, match_coordinates):
  is_primary_alignment, alignment_length, ref_id = process_primary_alignments(line)

  if is_primary_alignment:
    if process_long_alignments(alignment_type, alignment_length, ref_id, cutoff):
       items = line.split()
       query, q_fag = items[0].rsplit("_", 1) if items[0].count("_") == 2 else (items[0], 0)
       q_h = f"{query}+{ref_id}"
       shift = max(0, int(q_fag) * frag_size - overlap)
       start = int(items[2]) + shift
       end   = int(items[3]) + shift
       cord  = (start, end)

       if q_h not in match_coordinates:
        match_coordinates[q_h] = [cord]

       else:
        match_coordinates[q_h] = find_overlap(match_coordinates[q_h], cord)

  return match_coordinates

def process_alignment(alignment_paf, alignment_type):
  match_coordinates = {}
  for line in open(alignment_paf, 'r'):
    match_coordinates = process_alignment_line(line, alignment_type, match_coordinates)
  return match_coordinates

def assign_chromosomes(minimap_output, scaf_length):

  # for the largest scaffolds
  ref_match = process_alignment(minimap_output, 'Reference')

  # Assign 14 largest scaffolds into chromosomes
  scaf_sorted = sorted(scaf_length.items(), key=itemgetter(1), reverse=True)
  largest_14 = [k for k, v in scaf_sorted[:14]]

  # Calculate coverage for the reference wheat genome (7 x 3 chromosomes)
  scaf_cov = {s: [0] * 21 for s in largest_14}

  for p in ref_match:
    scaf, ref = p.split("+")

    if scaf in largest_14:
      final_ranges = merge_overlapping_ranges(ref_match[p])
      cov = sum([end - start + 1 for start, end in final_ranges])
      percov = cov / scaf_length[scaf]
      scaf_cov[scaf][ref_order.index(ref)] = round(percov, 3)

  scaf_dict = {}
  for scaf in scaf_cov:
    max_coverage = max(scaf_cov[scaf])
    max_index = scaf_cov[scaf].index(max_coverage)
    chromosome = ref_order[max_index]

    print(f'{scaf} assigned to {chromosome}: {round(max(scaf_cov[scaf]) * 100,3)}% coverage matches')
    scaf_dict[scaf] = chromosome

  return scaf_dict


def write_sequences(fasta_file, chr_dict, chloroplasts, mitochondira, outfile_prefix):

  records = list(SeqIO.parse(fasta_file, 'fasta'))

  if chr_dict != None:
    chr2scaf = {}
  
    for record in records:
        #Replace sequence IDs for large scaffolds
        if str(record.id) in chr_dict:
          chr2scaf[ chr_dict[str(record.id)] ] = str(record.id)
          record.id = chr_dict[str(record.id)]

  categorized_records = defaultdict(list)
  for record in records:
    if record.id in ref_order:
      categorized_records['chr'].append(record)
    elif record.id in chloroplasts:
      categorized_records['chlo'].append(record)
    elif record.id in mitochondria:
      categorized_records['mito'].append(record)
    else:
      categorized_records['other'].append(record)

  if chr_dict != None:
    with open(f'{outfile_prefix}.chromosomes.fa', 'w') as output_file:
      for item in categorized_records['chr']:
        output_file.write(f">{item.id} source:{chr2scaf[item.id]}\n{item.seq}\n")
  
      unplaced = ('N' * 200).join( str(item.seq) for item in categorized_records['other'] )
      output_file.write(f">Un\n{unplaced}\n")
  else:
    with open(f'{outfile_prefix}.genomic.fa', 'w') as output_file:
      for item in categorized_records['other']:
        SeqIO.write(item, output_file, 'fasta')

  with open(f'{outfile_prefix}.chloroplasts.fa', 'w') as chlo_file, open(
            f'{outfile_prefix}.mitochondira.fa', 'w') as mito_file:
    SeqIO.write(categorized_records['chlo'], chlo_file, 'fasta')
    SeqIO.write(categorized_records['mito'], mito_file, 'fasta')


# inputs should be: 
# 1) minimap paf output against plasmid
# 2) minimap paf output against reference genome (could be None) 
# 3) Input genome
# 4) Output prefix 

inputs = sys.argv[1:]

if len(inputs) != 4:
  print('Required inputs are: ')
  print('1) minimap paf output against plasmid')
  print('2) minimap paf output against reference genome (could be None)') 
  print('3) Input genome')
  print('4) Output prefix')  
  print('Additional input named "scaf.length" that stores the sequence IDs and sequence lengths will be used')

#inputs = ['minimap.plasmid.sorted.paf',
#        'minimap.ref.sorted.paf',
#        '../YaHS_scaffolds_final.fa',
#        'Kronos.collapsed' ]

# get scaffold length
scaf_len = {}

for line in open('scaf.length', 'r'):
  scaf, length = line.split()
  scaf_len[scaf] = int(length)


# collect plasmid sequences here
plasmid_output, genome_output, fasta_file, outfile_prefix = inputs
chloroplasts, mitochondria  = identify_plasmids(plasmid_output, scaf_len)
if genome_output != 'None':
  chr_dict = assign_chromosomes(genome_output, scaf_len)
else:
  chr_dict = None 
plasmids = chloroplasts + mitochondria
plasmids_records = write_sequences(fasta_file, chr_dict, chloroplasts, mitochondria, outfile_prefix)
