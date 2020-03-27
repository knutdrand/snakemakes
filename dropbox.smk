from snakemake.shell import shell
from snakemake.remote.dropbox import RemoteProvider as DropboxRemoteProvider

import pandas as pd
token = open(config["token"]).read().strip()
DBox = DropboxRemoteProvider(oauth2_access_token=token)
samples = pd.read_csv(config["sampleinfo"]).set_index("filename")
analysis = pd.read_csv(config["analysisinfo"]).set_index("Name")
print("in dropbox")

def get_export_files():
    return ["domain_coverage.csv", "human_domain_sizes.png", "species_domain_sizes.png", "gc_content_summary.tsv"]

def REMOTE_ADDRESS(filename):
    entry = samples.loc[filename]
    return config["project"]+entry["folder"] + "/" + filename + ".bed.gz"

rule import_data:
    input:
        lambda wildcards: DBox.remote(REMOTE_ADDRESS(wildcards.filename))
    output:
        "{species}/dedup/{filename}.bed.gz"
    shell:
        "cp '{input}' {output}"

rule export_results_to_dropbox:
    input:
        get_export_files()
    output:
        [DBox.remote(config["project"] + "snakemake/" + export_file)
         for export_file in get_export_files()]
    run:
        for i, o in zip(input, output):
            shell("cp '{i}' '{o}'")
