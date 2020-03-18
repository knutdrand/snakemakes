include: "mapping.smk"
include: "peakcalling.smk"
include: "nels.smk"
include: "trackhub.smk"

analysis_info = pd.read_csv(config["analysisinfo"]).set_index("Name")
species_dict = {"mouse": "mm10", "bovine": "bosTau8"}

def get_samples_for_analysis(analysis_name):
    if not analysis_name in analysis_info:
        return []
    return list(analysis_info[analysis_info[analysis_name]==1].index)

def get_species(sample):
    return species_dict[analysis_info.loc[sample].at["Species"]]

trackhub_species = [species_dict[s] for s in pd.unique(analysis_info[analysis_info["GB"]==1]["Species"])]

rule all:
    input:
        expand("qc/fastq_screen/{sample}.png", sample=get_samples_for_analysis("Screen")),
        [get_species(sample) + f"/motif_plots/{sample}.png" for sample in get_samples_for_analysis("Motif")],
        expand("trackhub/{species}/trackDb.txt", species=trackhub_species)

rule all_screen:
    input:
        expand("qc/fastq_screen/{sample}.png", sample=get_samples_for_analysis("Screen"))
        
        
