import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
get_name = lambda parts: f"{parts[0]}:{parts[1]}-{parts[2]}"
matches = {line.split("\t")[2] for line in open(snakemake.input[0]) if not line.startswith("#") and line.strip()}
hits = [get_name(line.split()) in matches for line in open(snakemake.input[1])]
ratio = np.cumsum(hits)/np.arange(1, len(hits)+1)
plt.plot(ratio)
plt.savefig(snakemake.output[0])
