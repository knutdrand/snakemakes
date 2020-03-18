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
        
rule all_motifs:
    input:    
        [get_species(sample) + f"/motif_plots/{sample}.png" for sample in get_samples_for_analysis("Motif")]

rule full_screen_table:
    input:
        expand("qc/fastq_screen/{sample}.txt", sample=get_samples_for_analysis("Screen"))
    output:
        "full_fastq_table.txt"
    script:
        "scripts/combinefiles.R"

rule screen_table:
    input:
        expand("qc/fastq_screen/{sample}.txt", sample=get_samples_for_analysis("Screen"))
    output:
        "fastq_table.txt"
    run:
        names = ["Species"]+ [i.split("/")[-1].split(".")[0] for i in input]
        humans_row = ["Human"] 
        mouse_row = ["Mouse"]
        for i in input:
            parts = [line.split() for line in open(i)]
            values = {p[0]: p[5] for p in parts  if p and p[0] in ("Mouse", "Human")}
            humans_row.append(values["Human"])
            mouse_row.append(values["Mouse"])
        with open(output[0], "w") as f:
            for line in [names, humans_row, mouse_row]:
                f.write("\t".join(line)+"\n")
