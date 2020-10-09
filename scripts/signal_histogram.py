from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
n_bins=200
def get_signals(filename):
    return (float(line.split()[6]) for line in open(filename))

def get_hist(values):
    h = np.zeros(n_bins, dtype="int")
    for v in values:
        h[min(int(v*10), n_bins-1)]+=1
    return h
signals = map(get_signals, snakemake.input)
hists = map(get_hist, signals)
lines = [plt.plot(h/np.sum(h))[0] for h in hists]
names = [Path(i).stem for i in snakemake.input]
plt.legend(lines, names)
plt.savefig(snakemake.output[0])
plt.savefig(snakemake.output[1])
