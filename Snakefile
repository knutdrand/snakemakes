from pathlib import Path
for name in config["analyses"]:
    include: f"{name}.smk"
include:
    "analysis.smk"
# if False:
#     include: "dropbox.smk"
# else:
#     include: "mapping.smk"
#     include: "nels.smk"
# include: "peakcalling.smk"
# include: "trackhub.smk"
# include: "donwloads.smk"
# include: "regions.smk"
# include: "gc.smk"


analysis_info = pd.read_csv(config["analysisinfo"]).set_index("Name")
print(analysis_info)
species_dict = {"mouse": "mm10", "bovine": "bosTau8", "human": "hg38", "porcine": "susScr11", "zebra": "danRer11"}

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
        expand("trackhub/{species}/trackDb.txt", species=trackhub_species),
        "mapping_stats.png",
        expand("trimming_effects/{sample}.png", sample=analysis_info.index)

rule all_screen:
    input:
        expand("qc/fastq_screen/{sample}.png", sample=get_samples_for_analysis("Screen"))
        
rule all_motifs:
    input:    
        [get_species(sample) + f"/motif_plots/{sample}.png" for sample in get_samples_for_analysis("Motif")]

rule all_gc:
    input:
        [get_species(sample) + f"/gc_content/{sample}.txt" for sample in get_samples_for_analysis("GC")],
    output:
        "gc_content_summary.tsv",
    shell:
        """
        echo "domains_minus_tss500\tnon_tss_containing_domains\tdomain_flanks\ttss_containing_domains\tname" > {output}
        cat {input} >> {output}
        """

rule all_coverage:
    input:
        [get_species(name) + f"/domain_coverage/{name}.txt" for name in get_samples_for_analysis("GB")]
    output:
        "domain_coverage.csv"
    run:
        with open(output[0], "w") as out_file:
            out_file.write(",".join(("Sample", "DomainCoverage", "Genomesize", "Prct"))+"\n")
            for in_file in input:
                c, g = open(in_file).read().strip().split()
                out_file.write(",".join([Path(in_file).stem, c, g, str(int(c)/int(g))])+"\n")

rule full_screen_table:
    input:
        expand("qc/fastq_screen/{sample}.txt", sample=get_samples_for_analysis("Screen"))
    output:
        "full_fastq_table.txt"
    script:
        "scripts/combinefiles.R"

rule human_domain_size_hist:
    input:
        [f"hg38/domains_logsize_hist/{sample}.npz" for sample in get_samples_for_analysis("human_comparison")]
    output:
        "human_domain_sizes.png"
    script:
        "scripts/peak_histograms.py"

rule species_domain_size_hist:
    input:
        [get_species(sample)+f"/domains_logsize_hist/{sample}.npz" for sample in get_samples_for_analysis("species_comparison")]
    output:
        "species_domain_sizes.png"
    script:
        "scripts/peak_histograms.py"


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

rule gzip_bed:
    input:
        "{filename}.bed"
    output:
        "{filename}.bed.gz"
    shell:
        "gzip {input} --keep"

names = get_samples_for_analysis("GB") # analysis_info.index # 
rule mapping_stats:
    input:
        reads=[f"reads/{sample}.fastq.gz.count" for sample in names],
        trimmed=[f"trimmed/{sample}.fastq.gz.count" for sample in names],
        mapped=[get_species(sample)+f"/mapped_filtered/{sample}.bam.count" for sample in names],
        dedup=[get_species(sample)+f"/dedup/{sample}.bam.count" for sample in names]
    output:
        "mapping_stats.csv"
    run:
        types = ["reads", "trimmed", "mapped", "dedup"]
        get_count = lambda f: int(open(f).read().strip())
        data = {t: [get_count(f) for f in getattr(input, t)] for t in types}
        pd.DataFrame(data, index=names).to_csv(output[0])

rule trimming_effects_all:
    input:
        expand("trimming_effects/{sample}.png", sample=names)

rule barchart:
    input:
        "{sample}.csv"
    output:
        report("{sample}.png", category="Counts")
    script:
        "scripts/barplot.R"
