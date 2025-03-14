#reformat gff files from apollo

in_gff  = 'Kronos_all.NLRs.gff3'
out_gff = 'Kronos_all.NLRs.reformatted.gff3'

gdict = {}
mdict = {}
cds_count = {}

with open(out_gff, 'w') as o:
  for line in open(in_gff, 'r'):
      if line.startswith('#'):
          continue
  
      fields = line.split("\t")
      chromosome, coordinates = fields[0].split("_")
      gstart, gend = map(int, coordinates.split("-"))
  
      fields[1] = 'Kronos_NLRs'
      fields[3] = int(fields[3]) + gstart - 1
      fields[4] = int(fields[4]) + gstart - 1
      fields[0] = chromosome
  
      if fields[2] == "gene":
          name = fields[8].split('Name=')[1].split(';')[0]
          description = fields[8].split('description=')[1].split(';')[0]
          gid = fields[8].split('ID=')[1].split(';')[0]
          gdict[gid] = name
  
          fields[8] = f'Name={name};ID={name};description={description}\n'
  
      elif fields[2] == "mRNA":
          name = fields[8].split('Name=')[1].split(';')[0]
          parent = fields[8].split('Parent=')[1].split(';')[0]
          mid = fields[8].split('ID=')[1].split(';')[0]
          mdict[mid] = name
  
          fields[8] = f'Name={name};ID={name};Parent={gdict[parent]}\n'
  
      elif fields[2] in ["exon", "CDS"]:
          parent = fields[8].split('Parent=')[1].split(';')[0]
  
          if parent not in cds_count:
              cds_count[parent] = [1, 1]
  
          if fields[2] == "exon":
              num = cds_count[parent][0]
              cds_count[parent][0] += 1
          elif fields[2] == "CDS":
              num = cds_count[parent][1]
              cds_count[parent][1] += 1
  
          fields[8] = f'Name={mdict[parent]}.{fields[2]}{num};ID={mdict[parent]}.{fields[2]}{num};Parent={mdict[parent]}\n'
      o.write("\t".join(map(str, fields)))
