include: "rules/common.smk"
from snakemake.shell import shell
from snakemake.remote.dropbox import RemoteProvider as DropboxRemoteProvider
import os
import pandas as pd

token = open(config["token"]).read().strip()
DBox = DropboxRemoteProvider(oauth2_access_token=token)
if "dropboxinfo" in config:
    samples = pd.read_csv(config["dropboxinfo"]).set_index("filename")
else:
    samples = pd.DataFrame()
sample_regexp = "|".join(samples.index).replace("+" ,"\+")

print("HEREE")
def REMOTE_ADDRESS(filename):
    print("#############", filename)
    return os.path.join(config["project"],
                        samples.loc[filename]["folder"],
                        f"{filename}.bed.gz")

def remote_output_address(filename):
    print("++++++++++++++", filename)
    return DBox.remote(os.path.join(config["project"], f"snakemake/{filename}"))

rule import_data:
    input:
        lambda wildcards: print(wildcards) or DBox.remote(REMOTE_ADDRESS(wildcards.filename))
    output:
        "{species}/dedup/{filename}.bed.gz"
    wildcard_constraints:
        filename=sample_regexp
    log:
        "logs/dropboximport/{species}/{filename}.log"
    shell:
        "cp '{input}' {output} 2> {log}"

rule export_file:
    input:
        "{filename}"
    output:
        DBox.remote(config["project"] + "snakemake/{filename}")
    shell:
        "cp '{input}' '{output}'"
        

# domains = expand_analysis("{species}/domains/{sample}.bed.gz", "GB")
# regions = [filename.format(region=region) for filename in expand_analysis("{species}/regions/{{region}}/{sample}.bed.gz", "GB")
#            for region in ("domain_flanks", "tss_containing_domains", "non_tss_containing_domains")]
# export_files = ["domain_coverage.csv", "human_domain_sizes.svg", "species_domain_sizes.svg", "gc_content_summary.tsv", "rn6/dedup_pe/rMII_K4.bed.gz"] + domains + regions
# 
# #TODO: Hack
# export_files =["mapping_stats_wo_alts.csv"]
# export_files = [f"hg38wo_alts/dedup/{sample}.bed.gz" for sample in analysis_info.index]
# export_files = get_all_corr(None)
#bad_pools = ["Day2B", "Day3B", "BlastB"]
#export_files = [f"hg38/clean_domains/{pool}_Easeq.bed" for pool in bad_pools]
# export_files = [f"hg38/dedup/BlastB_K4_catpool.bed.gz"] #  for pool in bad_pools for t in ("In", "K4")]
export_files =["domain_coverage.csv", "peak_coverage.csv"]
# 
# 
rule input_version_export:
    input:
        [remote_output_address(f) for f in export_files]
