from snakemake.shell import shell
counts = [int(open(f).read().strip()) for f in snakemake.input.counts]
lowest_count = min(counts)
factors = [lowest_count/c for c in counts]

for fastq_file, output_file, factor in zip(snakemake.input.fastqs, snakemake.output, factors):
    shell(f"seqtk sample {fastq_file} {factor} > {output_file}")

