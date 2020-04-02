import pandas as pd
from snakemake.shell import shell

NELS = "u1452@nelstor0.cbu.uib.no:/elixir-chr/nels/users/u1452/Projects/UiO_Dahl_Chromatin_2018/"
samples = pd.read_csv(config["nelsinfo"]).set_index("filename")
pools = set(samples["name"])

def get_pool(pool, read):
    return [sample for sample in sorted(samples[samples["name"]==pool].index) if sample.endswith(str(read))]
    

def REMOTE_ADDRESS(filename):
    entry = samples.loc[filename]
    return NELS+config["nelsfolder"]+entry["folder"] + "/" + filename + ".fastq.gz"
    
rule nels_import_data:
    output:
        "raw/{filename}.fastq.gz"
    run:
        remote_name = REMOTE_ADDRESS(wildcards.filename)
        shell("scp -i {config[nels_key]} {remote_name} {output}")

rule merge_nels:
    input:
        lambda wildcards: expand("raw/{sample}.fastq.gz", sample=get_pool(wildcards.pool, wildcards.read))
    output:
        "pe_reads/{pool}_{read}.fastq.gz"
    wildcard_constraints:
        pool="|".join(pools)
    shell:
        "zcat {input} | gzip > {output}"
    
rule nels_import_all:
    input:
        expand("raw/{filename}.fastq.gz", filename=samples.index)

def LANE_FILES(name):
    return sorted([i for i, e in  samples.loc[samples["name"]==name].iterrows()])

# rule merge_se:
#     input:
#         lambda wildcards: expand("raw/{sample}.fastq.gz", sample=LANE_FILES(wildcards.name))
#     output:
#         "reads/{name}.fastq.gz"
#     shell:
#         "zcat {input} | gzip > {output}"

