from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
names = [Path(i).stem for i in snakemake.input]
hists = [np.load(i) for i in snakemake.input]
min_xs = [np.flatnonzero(h["y"])[0] for h in hists]
min_x = min(min_xs)
lines = [plt.plot(h["x"][min_x:], h["y"][min_x:]/np.sum(h["y"]))[0] for h in hists]
plt.xscale("log")
plt.ylim(0, 0.03)
plt.legend(lines, names)
plt.savefig(snakemake.output[0])
