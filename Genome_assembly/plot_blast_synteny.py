# -*- coding: utf-8 -*-
"""
Visualizing Synteny Plots for BLAST Output (Tabular Format 6)
"""

import argparse
import matplotlib.pyplot as plt

def resize(num, param):
    """Scale a number by a resize parameter."""
    return int(num) / param

parser = argparse.ArgumentParser(description='Read a BLAST output in tabular format (6) and produce a synteny plot in SVG format.')
parser.add_argument('blast_file', help='The path to a BLAST output file.')
parser.add_argument('query', help='The sequence name for a query to be visualized.')
parser.add_argument('hit', help='The sequence name for a hit to be compared to.')
parser.add_argument('--hstart', type=int, required=True, help='Starting coordinate of the hit visualization (Y-axis).')
parser.add_argument('--hend', type=int, required=True, help='Ending coordinate of the hit visualization (Y-axis).')
parser.add_argument('--alignment_length', type=int, default=500, help='Minimum alignment length to be processed. Default is 500.')
parser.add_argument('--resize_param', type=int, default=1000000, help='Parameter to resize the plot (in Mb).')


args = parser.parse_args()


query_start = 0
query_end = 20000000


start_x, start_y, end_x, end_y = [], [], [], []

with open(args.blast_file, 'r') as file:
    for line in file:
        items = line.split()

        
        qseqid = items[0]    # Query sequence ID
        sseqid = items[1]    # Subject sequence ID
        length = int(items[3])    # Alignment length
        qstart = int(items[6])    # Query start
        qend = int(items[7])      # Query end
        sstart = int(items[8])    # Subject start
        send = int(items[9])      # Subject end

        
        if qseqid == args.query and sseqid == args.hit and length >= args.alignment_length:
            start_x.append(resize(qstart, args.resize_param))
            end_x.append(resize(qend, args.resize_param))
            start_y.append(resize(sstart, args.resize_param))
            end_y.append(resize(send, args.resize_param))

fig, ax = plt.subplots()

# Set fixed visualization frame for the X-axis (0 - 20Mb)
ax.set_xlim(resize(query_start, args.resize_param), resize(query_end, args.resize_param))
ax.set_ylim(resize(args.hstart, args.resize_param), resize(args.hend, args.resize_param))

for sx, sy, ex, ey in zip(start_x, start_y, end_x, end_y):
    ax.plot([sx, ex], [sy, ey], color='k')

ax.set_xlabel(f'{args.query} (Mb)')
ax.set_ylabel(f'{args.hit} ({args.hstart}-{args.hend} bp)')

output_dir = '/global/scratch/users/jiaqitang0422/Kronos/Chinese_spring/start_end'
output_filename = f'{output_dir}/{args.query}_vs_{args.hit}.png'
plt.savefig(output_filename, format='png')

print(f"Plot successfully saved as: {output_filename}")
