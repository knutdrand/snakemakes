import pandas as pd

if "methylationinfo" in config:
    methyl_samples = pd.read_csv(config["methylationinfo"])
else:
    methyl_samples = pd.DataFrame.from_dict({"Accession": [], "Title": []})
methyl_samples = methyl_samples.set_index("Accession")

def get_celltype_filenames(wildcards):
    mask = [wildcards.celltype in title for title in methyl_samples["Title"]]
    samples = methyl_samples[mask]
    return [f"hg19/methylation/{gsm}_{name}.Cmet.bed.gz" 
            for gsm, name in zip(samples.index, samples["Title"])]

def get_download_link(gsm, name):
    folder = gsm[:7]
    return f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{folder}nnn/{gsm}/suppl/{gsm}_{name}.Cmet.bed.gz"

rule aggregate_cell_type:
    input:
        get_celltype_filenames
    output:
        "hg19/methylation_aggr/{celltype}.tsv"
    script:
        "scripts/merge_methylation.py"

rule convert_to_bdg:
    input:
        "hg19/methylation_aggr/{celltype}.tsv"
    output:
        "hg19/methylation_aggr/{celltype}.bdg"
    shell:
        """awk '{{OFS="\t"}}{{if ($3>4) print $1,$2,$2+1,$4/$3}}' {input} > {output}"""

rule download_methyl:
    output:
        "hg19/methylation/{gsm}_{name}.Cmet.bed.gz"
    run:
        link = get_download_link(wildcards.gsm, wildcards.name)
        shell("wget %s -O {output}" % link)
