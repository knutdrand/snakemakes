import gzip
sizes = (len(line.strip()) for i, line in enumerate(gzip.open(snakemake.input[0], "rb")) if i % 4 == 1)
