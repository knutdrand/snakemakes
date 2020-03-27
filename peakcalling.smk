from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

macs_output=["treat_pileup.bdg", "control_lambda.bdg", "peaks.narrowPeak"]
broad_output=["treat_pileup.bdg", "control_lambda.bdg", "peaks.broadPeak"]
genome_sizes = {"mm10": "mm", "hg38": "hs"}
jaspar_address = "http://jaspar.genereg.net/api/v1/matrix/"


analysis_info = pd.read_csv(config["analysisinfo"]).set_index("Name")

def MOTIF_FILE(tf):
    motifs = {"oct4": "MA0142.1",
              "nanog": "UN0383.1",
              "sox2": "MA0143.1"}
    return "motives/{m}.meme".format(m=motifs[tf.lower()])

def INPUT(name):
    if not name in analysis_info:
        return None
    return analysis_info.loc[name].at["Input"]

# rule callpeak_wo_control:
#     input:
#         "{species}/dedup/{sample}.bed"
#     output:
#         expand("{{species}}/peakcalling/{{sample}}_{filetype}", filetype=macs_output)
#     run:
#         genomesize=genome_sizes.get(wildcards.species, 2913022398)
#         shell("""macs2 callpeak -t {input} -g {genomesize} --bdg --outdir {wildcards.species}/peakcalling -n {wildcards.sample}""")
# 

def macs_input(wildcards):
    print("HERE")
    i = ["{species}/dedup/{sample}.bed.gz"]
    print(i)
    if not pd.isnull(INPUT(wildcards.sample)):
        i.append("{species}/dedup/%s.bed.gz" % INPUT(wildcards.sample))
    print(i)
    return i

rule callpeak:
    input:
        macs_input
    output:
        expand("{{species}}/peakcalling/{{sample}}_{filetype}", filetype=macs_output)
    run:
        genomesize=genome_sizes.get(wildcards.species, 2913022398)
        control="-c {input[1]}" if not pd.isnull(INPUT(wildcards.sample)) else ""
        shell("""macs2 callpeak -t {input[0]} %s -g {genomesize} --bdg --outdir {wildcards.species}/peakcalling -n {wildcards.sample}""" % control)

rule call_broad_peak:
    input:
        macs_input
    output:
        expand("{{species}}/broadpeakcalling/{{sample}}_{filetype}", filetype=broad_output)
    run:
        genomesize=genome_sizes.get(wildcards.species, 2913022398)
        control="-c {input[1]}" if not pd.isnull(INPUT(wildcards.sample)) else ""
        shell("""macs2 callpeak -t {input[0]} %s -g {genomesize} --bdg --broad --outdir {wildcards.species}/broadpeakcalling -n {wildcards.sample}""" % control)

rule merge_domains:
    input:
        "{species}/broadpeakcalling/{name}_peaks.broadPeak"
    output:
        "{species}/domains/{name}.bed"
    shell:
        "bedtools merge -d 5000 -i {input} > {output}"

rule sort_peaks:
    input:
        "{species}/peakcalling/{sample}_peaks.narrowPeak"
    output:
        "{species}/sorted_peaks/{sample}.narrowPeak"
    shell:
        "sort -nr -k5 {input} > {output}"

rule get_peak_sequences:
    input:
        peaks="{species}/peakcalling/{sample}_peaks.narrowPeak",
        reference=config["data_dir"]+"/{species}/{species}.fa"
    output:
        "{species}/peak_fasta/{sample}.fa"
    shell:
        "bedtools getfasta -fi {input.reference} -bed {input.peaks} > {output}"

rule motif_enrichment:
    input:
        "{species}/peak_fasta/{sample}.fa", 
        lambda wildcards: MOTIF_FILE(analysis_info.loc[wildcards.sample].at["AB"])
    output:
        multiext("{species}/motif_matches/{sample}/fimo", ".html", ".xml", ".tsv", ".gff") 
    shell:
        "fimo --oc {wildcards.species}/motif_matches/{wildcards.sample}/ {input[1]} {input[0]} "
    
rule motif_plot:
    input:
        "{species}/motif_matches/{sample}/fimo.tsv",
        "{species}/sorted_peaks/{sample}.narrowPeak"
    output:
        report("{species}/motif_plots/{sample}.png", category="Motif_plots")
    run:
        get_name = lambda parts: f"{parts[0]}:{parts[1]}-{parts[2]}"
        matches = {line.split("\t")[2] for line in open(input[0]) if not line.startswith("#") and line.strip()}
        hits = [get_name(line.split()) in matches for line in open(input[1])]
        #hits = [(line[1:].strip() in matches) for line in open(input[1]) if line.startswith(">")]
        ratio = np.cumsum(hits)/np.arange(1, len(hits)+1)
        plt.plot(ratio)
        print(output[0])
        plt.savefig(output[0])

rule get_meme:
    output:
        "motives/{name}"
    shell:
        'wget {jaspar_address}{wildcards.name} -O {output} --user-agent="Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.8.1.6) Gecko/20070725 Firefox/2.0.0.6"'

rule get_domain_coverage:
    input:
        "{species}/domains/{name}.bed",
        "{species}/data/chrom.sizes.txt"
    output:
        "{species}/domain_coverage/{name}.txt"
    shell:
        """
        awk '{{t+=($3-$2)}}END{{print t}}' {input[0]} > {output}
        awk '{{t+=$2}}END{{print t}}' {input[1]} >> {output}
        """
