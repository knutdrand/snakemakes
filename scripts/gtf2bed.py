import sys

all_parts = (line.split() for  line in sys.stdin if not line.startswith("#"))
tuples = (("chr"+parts[0], int(parts[3])-1, int(parts[4])-1, (parts[13] if parts[13]!='"ensembl";' else parts[9])[1:-2], parts[6])
          for parts in all_parts if parts[2]=="gene")
for t in tuples:
    print ("%s\t%s\t%s\t%s\t0\t%s" % t)
    
