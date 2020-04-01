import pandas as pd
from snakemake.shell import shell
genome_sizes = {"mm10": "mm", "hg38": "hs"}
genomesize=genome_sizes.get(snakemake.wildcards.species, 2913022398)
analysis_info = pd.read_csv(snakemake.config["analysisinfo"]).set_index("Name")
extra = snakemake.params.get("extra", "")
def INPUT(name):
    if not name in analysis_info.index:
        return None
    return analysis_info.loc[name].at["Input"]

control="-c {snakemake.input[1]}" if not pd.isnull(INPUT(snakemake.wildcards.sample)) else ""

shell("""macs2 callpeak -t {snakemake.input[0]} %s {extra} --bdg --broad --outdir {snakemake.wildcards.species}/broadpeakcalling -n {snakemake.wildcards.sample}""" % control)
