from snakemake.remote.dropbox import RemoteProvider as DropboxRemoteProvider
import pandas as pd
token = open(config["token"]).read().strip()
DBox = DropboxRemoteProvider(oauth2_access_token=token)
samples = pd.read_csv(config["sampleinfo"]).set_index("filename")
analysis = pd.read_csv(config["analysisinfo"]).set_index("Name")
print("in dropbox")
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
