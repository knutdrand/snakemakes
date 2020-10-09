import pandas as pd

wildcard_constraints:
    species="|".join(config.get("species", ["mm10", "hg38"])),
    from_species="[^/]+",
    to_species="[^/]+",
    read="1|2",
    sample="[^/]+",
    sra="[A-Z]+\d+"

if "analysisinfo" in config:
    analysis_info = pd.read_csv(config["analysisinfo"]).set_index("Name")
else:
    analysis_info = pd.DataFrame({"Name": [], "GB": [], "Species": []})
print(analysis_info)
species_dict = {"mouse": "mm10", "bovine": "bosTau8", "human": "hg38", "porcine": "susScr11", "zebra": "danRer11", "rat": "rn6"}
inverted_species_dict = {v: k for k, v in species_dict.items()}
chromosome_grep = "grep -Ew -e 'chr[0-9]{{1,2}}' -e chrX -e chrY"

def get_samples_for_analysis(analysis_name):
    if not analysis_name in analysis_info:
        return analysis_info.index
    return list(analysis_info[analysis_info[analysis_name]==1].index)

def get_species(sample):
    return species_dict[analysis_info.loc[sample].at["Species"]]

def expand_analysis(format_str, analysis_name):
    samples = get_samples_for_analysis(analysis_name)
    return [format_str.format(sample=sample, species=get_species(sample)) for sample in samples]

def expand_subsampled_analysis(format_str, analysis_name):
    samples = get_samples_for_analysis(analysis_name)
    return [format_str.format(sample=sample, species=get_species(sample)) for sample in samples]


def get_species_tracks(species):
    d = analysis_info[analysis_info["Species"]==inverted_species_dict[species]]
    print("Species")
    print(d)
    d = d[d["GB"]==1]
    print(d)
    return d.index

def get_all_corr(wildcards):
    if not "correlation" in analysis_info:
        return []
    values = [(species_dict[e["Species"]], i, e["correlation"])
              for i, e in analysis_info.iterrows() if e["correlation"]]
    return [f"{species}/correlations/windows/10000/{sample_a}_VS_{sample_b}.png"
            for species, sample_a, sample_b in values]
#    return [f"{species}/correlations/tss/500/{sample_a}_VS_{sample_b}.png"
#            for species, sample_a, sample_b in values]

def pyplot(code):
    return """python -c "import numpy as np; import matplotlib.pyplot as plt; %s" """ % code
