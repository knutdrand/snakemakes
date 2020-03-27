import gzip
import numpy as np
from pathlib import Path
path = Path(snakemake.input[0])

if path.suffix == ".gz":
    f = gzip.open(path, "rb")
    file_format = path.suffixes[-2]
else:
    f = open(path)
    file_format = path.suffix
if file_format == ".bed":
    parts = (line.split() for line in f)
    sizes = (int(p[2])-int(p[1]) for p in parts)
elif file_format == ".fastq":
    sizes = (len(line.strip()) for i, line in enumerate(f) if i % 4 == 1)
else:
    assert False, (path, file_format)

n_bins = snakemake.params.get("nbins", 250)
bin_size = snakemake.params.get("bin_size", 0.05)
trans = np.log if  snakemake.params.get("do_log", True) else lambda x: x

log_sizes = (trans(size) for size in sizes)
size_bins = (min(s//bin_size, n_bins-1) for s in log_sizes)
bins = np.zeros(n_bins, dtype="float")

for sb in size_bins:
    bins[int(sb)] += 1
    x = np.exp(bin_size*np.arange(n_bins))
    np.savez(snakemake.output[0], x=x, y=bins)
