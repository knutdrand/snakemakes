import numpy as np
import matplotlib.pyplot as plt
from pathlib import PurePath

names = snakemake.params["samples"]
arrays = [np.load(filename) for filename in snakemake.input]
print([a.shape for a in arrays])
lines = [plt.plot(array)[0] for array in arrays]
plt.legend(lines, names)
plt.savefig(snakemake.output[0])
