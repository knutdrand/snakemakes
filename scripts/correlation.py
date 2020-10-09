import numpy as np
import matplotlib.pyplot as plt
from pandas import Series
def get_counts(filename):
    return np.array([int(line.split()[-1]) for line in open(filename)], dtype="int")

counts_a, counts_b = [get_counts(f) for f in snakemake.input]
mask = np.logical_or(counts_a, counts_b)
masked_a, masked_b = (counts_a[mask], counts_b[mask])
plt.scatter(masked_a, masked_b)
plt.xlabel(snakemake.input[0])
plt.ylabel(snakemake.input[1])
r = Series(masked_a).corr(Series(masked_b))
plt.title(f"r={r}")
plt.savefig(snakemake.output[0])
