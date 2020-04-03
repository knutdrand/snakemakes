from snakemake.shell import shell
from snakemake.remote.dropbox import RemoteProvider as DropboxRemoteProvider

import pandas as pd
token = open(config["token"]).read().strip()
DBox = DropboxRemoteProvider(oauth2_access_token=token)
samples = pd.read_csv(config["dropboxinfo"]).set_index("filename")
analysis = pd.read_csv(config["analysisinfo"]).set_index("Name")
sample_regexp = "|".join(samples.index).replace("+" ,"\+")
species_dict = {"mouse": "mm10", "bovine": "bosTau8", "human": "hg38", "porcine": "susScr11", "zebra": "danRer11", "rat": "rn6"}
get_species = lambda sample: species_dict[analysis.loc[sample].at["Species"]]

def get_export_files():
    samples = analysis[analysis["GB"]==1].index
    domains = [get_species(sample)+f"/domains/{sample}.bed.gz" for sample in samples]
    regions = ["domain_flanks", "tss_containing_domains", "non_tss_containing_domains"]
    regions = [get_species(sample)+f"/regions/{region}/{sample}.bed.gz" for sample in samples for region in regions]
    reads = ["rn6/dedup_pe/rMII_K4.bed.gz"]
    return ["domain_coverage.csv", "human_domain_sizes.svg", "species_domain_sizes.svg", "gc_content_summary.tsv"] +domains + regions + reads

def REMOTE_ADDRESS(filename):
    print(filename)
    entry = samples.loc[filename]
    return config["project"]+entry["folder"] + "/" + filename + ".bed.gz"

def remote_output_address(filename):
    print(DBox.remote(config["project"] + f"snakemake/{filename}"))
    return DBox.remote(config["project"] + f"snakemake/{filename}")

rule import_data:
    input:
        lambda wildcards: DBox.remote(REMOTE_ADDRESS(wildcards.filename))
    output:
        "{species}/dedup/{filename}.bed.gz"
    wildcard_constraints:
        filename=sample_regexp
    shell:
        "cp '{input}' {output}"

# rule export_results_to_dropbox:
#     input:
#         get_export_files()
#     output:
#         [DBox.remote(config["project"] + "snakemake/" + export_file)
#          for export_file in get_export_files()]
#     run:
#         for i, o in zip(input, output):
#             shell("cp '{i}' '{o}'")

rule export_file:
    input:
        "{filename}"
    output:
        DBox.remote(config["project"] + "snakemake/{filename}")
    shell:
        "cp '{input}' '{output}'"
        
rule input_version_export:
    input:
        [remote_output_address(f) for f in get_export_files()]
        # remote_output_address("rn6/domains/rMII_K4.bed")

rule new_mouse_pool:
    input:
        "mm10/dedup/mm_MII_{type}_r1.bed.gz",
        "mm10/dedup/mm_MII_{type}_r2.bed.gz"
    output:
        "mm10/dedup/mm_MII_newpool_{type}.bed.gz",
    shell:
        "zcat {input} | gzip > {output}"
