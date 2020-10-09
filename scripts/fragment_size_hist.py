import gzip
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
path = Path(snakemake.input[0])

if path.suffix == ".gz":
    f = gzip.open(path, "rb")
else:
    f = open(path)
parts = (line.split() for line in f)
sizes = (int(p[2])-int(p[1]) for p in parts)

n_bins = snakemake.params.get("nbins", 1000)
bin_size = snakemake.params.get("bin_size", 1)
size_bins = (min(s//bin_size, n_bins-1) for s in sizes)
bins = np.zeros(n_bins, dtype="float")

for sb in size_bins:
    bins[int(sb)] += 1
np.save(snakemake.output.data, bins)
plt.plot(bins)
plt.savefig(snakemake.output.fig)
