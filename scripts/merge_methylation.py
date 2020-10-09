import heapq
from itertools import groupby, dropwhile
from operator import itemgetter
from functools import reduce, partial
import gzip

sort_order = """chr7
chr20
chr22
chr14
chrY
chr19
chr8
chr1
chr11
chr6
chr17
chr21
chr16
chr18
chr3
chr12
chr15
chrX
chr4
chrM
chr2
chr9
chr13
chr10
chr5
""".strip().split()
indexes = {chrom: i for i, chrom in enumerate(sort_order)}

def parse_file(lines):
    lines = dropwhile(lambda line: line.startswith("#"), lines)
    split_lines = (line.split()[:6] for line in lines)
    return (((indexes[chrom], int(pos)), (int(total), int(meth)))
            for chrom, pos, _, _, total, meth in split_lines)


entries_list = [parse_file(gzip.open(filename, "rt")) for filename in snakemake.input]
sorted_entries = heapq.merge(*entries_list)
grouped = groupby(sorted_entries, itemgetter(0))
new_entries = ((pos, reduce(lambda a, b: (a[0]+b[0], a[1]+b[1]), (e[1] for e in entries)))
               for pos, entries in grouped)

with open(snakemake.output[0], "w") as f:
    for (chrom, pos), (total, meth) in new_entries:
        f.write(f"{sort_order[chrom]}\t{pos}\t{total}\t{meth}\n")
