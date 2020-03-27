import numpy as np
import matplotlib.pyplot as plt

rule size_hist:
    input:
        "{path}/{name}.bed"
    output:
        "{path}_logsize_hist/{name}.npz"
    wildcard_constraints:
        name="[^/]+"
    run:
        n_bins = 250
        bin_size = 0.05
        parts = (line.split() for line in open(input[0]))
        log_sizes = (np.log(int(p[2])-int(p[1])) for p in parts)
        size_bins = (min(s//bin_size, n_bins-1) for s in log_sizes)
        bins = np.zeros(n_bins, dtype="float")
        for sb in size_bins:
            bins[int(sb)] += 1
        bins /= np.sum(bins)
        x = np.exp(bin_size*np.arange(n_bins))
        np.savez(output[0], x=x, y=bins)
