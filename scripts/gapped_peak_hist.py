import numpy as np
N = 3000
split_lines = (line.split() for line in open(snakemake.input[0]))
ds = ((int(s[2])-int(s[1]), [int(d) for d in s[10].split(",")]) for s in split_lines)
blocks = (block for p, blocks in ds for block in blocks if block not in (1, p))
hist = np.zeros(N, dtype="int")
for block in blocks:
    hist[min(block, N-1)] += 1

np.save(snakemake.output[0], hist)
