import pandas as pd
from snakemake.shell import shell

NELS = "u1452@nelstor0.cbu.uib.no:/elixir-chr/nels/users/u1452/Projects/UiO_Dahl_Chromatin_2018/"
samples = pd.read_csv(config["sampleinfo"]).set_index("filename")

def REMOTE_ADDRESS(filename):
    entry = samples.loc[filename]
    return NELS+config["project"]+entry["folder"] + "/" + filename + ".fastq.gz"
    
rule import_data:
    output:
        "raw/{filename}.fastq.gz"
    run:
        remote_name = REMOTE_ADDRESS(wildcards.filename)
        shell("scp -i {config[nels_key]} {remote_name} {output}")
