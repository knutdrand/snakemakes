include: "rules/common.smk"

from pathlib import Path
for name in config.get("analyses", ["mapping", "bdgplot", "trackhub", "peakcalling", "donwloads", "regions"]):
    include: f"{name}.smk"

trackhub_species = [species_dict[s] for s in pd.unique(analysis_info[analysis_info["GB"]==1]["Species"])]

rule all:
    input:
        expand_analysis("qc/fastq_screen/{sample}.png", "Screen"),
        expand_analysis("{species}/motif_plots/{sample}.png", "Motif"),
        expand("trackhub/{species}/trackDb.txt", species=trackhub_species),
        "mapping_stats.png",
        expand("trimming_effects/{sample}.png", sample=analysis_info.index)

rule all_screen:
    input:
        expand_analysis("qc/fastq_screen/{sample}.png","Screen")
        
rule all_motifs:
    input:    
        expand_analysis("{species}/motif_plots/{sample}.png", "Motif")

rule all_gc:
    input:
        expand_analysis("{species}/gc_content/{sample}.txt", "GC")
    output:
        "gc_content_summary.tsv",
    shell:
        """
        echo "domains_minus_tss500\tnon_tss_containing_domains\tdomain_flanks\ttss_containing_domains\tname" > {output}
        cat {input} >> {output}
        """

rule all_coverage:
    input:
        expand_analysis("{species}/{{region}}_coverage/{sample}.txt", "GB")
    output:
        "{region}_coverage.csv"
    run:
        with open(output[0], "w") as out_file:
            out_file.write(",".join(("Sample", "DomainCoverage", "Genomesize", "Prct"))+"\n")
            for in_file in input:
                c, g = open(in_file).read().strip().split()
                out_file.write(",".join([Path(in_file).stem, c, g, str(int(c)/int(g))])+"\n")

rule full_screen_table:
    input:
        expand_analysis("qc/fastq_screen/{sample}.txt", "Screen")
    output:
        "full_fastq_table.txt"
    script:
        "scripts/combinefiles.R"

rule human_domain_size_hist:
    input:
        expand_analysis("hg38/domains_logsize_hist/{sample}.npz", "human_comparison")
    output:
        multiext("human_domain_sizes.svg", ".svg", ".png")
    script:
        "scripts/peak_histograms.py"

rule human_domain_size_hist_adj:
    input:
        expand_analysis("hg38/domains_adjsize_hist/{sample}.npz", "human_comparison")
    output:
        multiext("human_adj_domain_sizes", ".svg", ".png")
    script:
        "scripts/peak_histograms.py"

rule species_domain_size_hist_adj:
    input:
        expand_analysis("{species}/domains_adjsize_hist/{sample}.npz", "species_comparison")
    output:
        multiext("species_adj_domain_sizes", ".svg", ".png")
    script:
        "scripts/peak_histograms.py"


rule species_domain_size_hist:
    input:
        expand_analysis("{species}/domains_logsize_hist/{sample}.npz", "species_comparison")
    output:
        multiext("species_domain_sizes", ".svg", ".png")
    script:
        "scripts/peak_histograms.py"

rule species_peak_size_hist:
    input:
        expand_analysis("{species}/broadpeakcalling_hist_peaks/{sample}.npz", "species_comparison")
    output:
        multiext("species_peak_sizes", ".svg", ".png")
    script:
        "scripts/peak_histograms.py"

rule human_peak_size_hist:
    input:
        expand_analysis("{species}/broadpeakcalling_hist_peaks/{sample}.npz", "human_comparison")
    output:
        multiext("human_peak_sizes", ".svg", ".png")
    script:
        "scripts/peak_histograms.py"

rule species_gapped_peak_size_hist:
    input:
        expand_analysis("{species}/broadpeakcalling_hist_gapped_peaks/{sample}.npz", "species_comparison")
    output:
        multiext("species_gapped_peak_sizes", ".svg", ".png")
    script:
        "scripts/peak_histograms.py"

rule human_gapped_peak_size_hist:
    input:
        expand_analysis("{species}/broadpeakcalling_hist_gapped_peaks/{sample}.npz", "human_comparison")
    output:
        multiext("human_gapped_peak_sizes", ".svg", ".png")
    script:
        "scripts/peak_histograms.py"

rule human_signal_hist:
    input:
        expand_analysis("{species}/broadpeakcalling/{sample}_peaks.broadPeak", "human_comparison")
    output:
        multiext("human_signal_hist", ".svg", ".png")
    script:
        "scripts/signal_histogram.py"

rule species_signal_hist:
    input:
        expand_analysis("{species}/broadpeakcalling/{sample}_peaks.broadPeak", "species_comparison")
    output:
        multiext("species_signal_hist", ".svg", ".png")
    script:
        "scripts/signal_histogram.py"

rule species_heatplots:
    input:
        [f"{get_species(sample)}/heatplots/{sample}_{kind}.png" for sample in get_samples_for_analysis("species_comparison") for kind in ("treat_pileup", "control_lambda", "qvalues")]

rule human_heatplots:
    input:
        [f"hg38/heatplots/{sample}_{kind}.png" for sample in get_samples_for_analysis("human_comparison") for kind in ("treat_pileup", "control_lambda", "qvalues")]


rule screen_table:
    input:
        expand_analysis("qc/fastq_screen/{sample}.txt", "Screen")
    output:
        "fastq_table.txt"
    script:
        "scripts/screen_table.py"
        
rule gzip_bed:
   input:
       "{species}/{filename}.bed"
   output:
       "{species}/{filename}.bed.gz"
   shell:
       "gzip {input} --keep"

names = get_samples_for_analysis("GB") # analysis_info.index # 
print(names)
rule mapping_stats:
    input:
        reads=[f"reads/{sample}.fastq.gz.count" for sample in names],
        trimmed=[f"trimmed/{sample}.fastq.gz.count" for sample in names],
        mapped=expand_analysis("{species}/mapped_filtered/{sample}.bam.count", "GB"),
        dedup=expand_analysis("{species}/dedup/{sample}.bam.count", "GB")
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
        report("{sample}_bars.png", category="Counts")
    script:
        "scripts/barplot.R"

rule bedgz_count:
    input:
        "{filename}.bed.gz"
    output:
        "{filename}.bed.gz.count"
    shell:
        "z{chromosome_grep} {input} | wc -l > {output}"

rule genome_size:
    input:
        "{species}/data/chrom.sizes.txt"
    output:
        "{species}/data/genome_size.txt"
    shell:
        "awk '{{t+=$2}}END{{print t}}' {input} > {output}"
