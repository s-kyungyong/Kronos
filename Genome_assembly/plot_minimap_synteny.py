import argparse
import matplotlib.pyplot as plt

def resize(num, param):
    """Scale a number by a resize parameter."""
    return [int(num) / param]

parser = argparse.ArgumentParser(description='Read a PAF file from minimap and produce a synteny plot in SVG format.')
parser.add_argument('paf', help='The path to a PAF file from minimap.')
parser.add_argument('query', help='The sequence name for a query to be visualized.')
parser.add_argument('hit', help='The sequence name for a hit to be compared to.')
parser.add_argument('--start', type=int, default=0, help='Starting coordinate of a query to be visualized. Default is 0.')
parser.add_argument('--end', type=int, default=None, help='Ending coordinate of a query to be visualized. Default is the last position of the query.')
parser.add_argument('--hstart', type=int, default=None, help='Starting coordinate of a hit to be visualized. Defaults to query start.')
parser.add_argument('--hend', type=int, default=None, help='Ending coordinate of a hit to be visualized. Defaults to query end.')
parser.add_argument('--alignment_length', type=int, default=5000, help='Minimum alignment length to be processed. Default is 50,000.')
parser.add_argument('--agp', default=None, help='agp file containing scaffolding information')
args = parser.parse_args()

fig_type = 'global' if args.start == 0 and args.end is None else f'local.{args.start}-{args.end}'

alignment_type = "tp:A:P"
quality = 30
resize_param = 1000000  # Mb

start_x, start_y, end_x, end_y = [], [], [], []

with open(args.paf, 'r') as file:
    for line in file:
        items = line.split()
        q, h = items[0], items[5]

        if q == args.query and h == args.hit:
            if args.end is None:
                args.end = int(items[1])
                args.hend = int(items[6])

            start_y += resize(items[2], resize_param)  # Swapped x and y
            end_y += resize(items[3], resize_param)  # Swapped x and y

            if items[4] == "+":
                start_x += resize(items[7], resize_param)  # Swapped x and y
                end_x += resize(items[8], resize_param)  # Swapped x and y
            else:
                start_x += resize(items[8], resize_param)  # Swapped x and y
                end_x += resize(items[7], resize_param)  # Swapped x and y

vline = []
if not args.agp == None:
  for line in open(args.agp, 'r'):
    items = line.split('\t')
    if items[0] == args.query and items[6] == 'scaffold':
      vline.append( (int(items[1]) + 49) / resize_param )

if args.hstart is None: args.hstart = args.start
if args.hend is None: args.hend = args.end

fig, ax = plt.subplots()

ax.set_xlim(args.hstart / resize_param, args.hend / resize_param)  # Swapped x and y
ax.set_ylim(args.start / resize_param, args.end / resize_param)  # Swapped x and y

for sx, sy, ex, ey in zip(start_x, start_y, end_x, end_y):
    ax.plot([sx, ex], [sy, ey], color='k')

if len(vline) >= 1:
  plt.vlines(x=vline, ymin=args.start / resize_param, ymax=args.end / resize_param, color='g')  # Swapped x and y

ax.set_xlabel(f'{args.hit} (Mb)')  # Swapped x and y labels
ax.set_ylabel(f'{args.hit} (Mb)')  # Changed y-axis label to match x-axis
ax.grid(True)

plt.savefig(f'{args.query}_vs_{args.hit}.{fig_type}.svg', format='svg')
