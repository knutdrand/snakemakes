include: "rules/common.smk"
import pandas as pd
import os
macs_output=["treat_pileup.bdg", "control_lambda.bdg", "peaks.narrowPeak"]
broad_output=["treat_pileup.bdg", "control_lambda.bdg", "peaks.broadPeak"]
genome_sizes = {"mm10": "mm", "hg38": "hs"}
jaspar_address = "http://jaspar.genereg.net/api/v1/matrix/"

def MOTIF_FILE(tf):
    motifs = {"oct4": "MA0142.1",
              "nanog": "UN0383.1",
              "sox2": "MA0143.1"}
    return "motives/{m}.meme".format(m=motifs[tf.lower()])


INPUT = lambda name: analysis_info.loc[name].at["Input"] if name in analysis_info.index else None

def macs_input(wildcards):
    i = ["{species}/dedup/{sample}.bed.gz"]
    name = wildcards.sample.split("subsampled_")[-1]
    if not pd.isnull(INPUT(name)):
        i.append("{species}/dedup/%s.bed.gz" % INPUT(name))
    return i

def subsampled_macs_input(wildcards):
    i = ["{species}/subsampled/{n_reads}/{sample}.bed.gz"]
    if not pd.isnull(INPUT(wildcards.sample)):
        i.append("{species}/dedup/%s.bed.gz" % INPUT(wildcards.sample))
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
    conda:
        "envs/oldmacs.yaml"
    script:
        "scripts/macscall.py"

rule call_subsampled_broad_peak:
    input:
        subsampled_macs_input
    output:
        expand("{{species}}/subsampled_broadpeakcalling/{{n_reads}}/{{sample}}_{filetype}", filetype=broad_output)
    conda:
        "envs/oldmacs.yaml"
    script:
        "scripts/macscall.py"


rule frag_peak_call:
    input:
        macs_input
    output:
        expand("{{species}}/broadpeakcalling/{{sample}}_f{{fraglen}}bp_{filetype}", filetype=broad_output)
    conda:
        "envs/oldmacs.yaml"
    params:
        extra="--extsize {snakemake.wildcards.fraglen}"
    script:
        "scripts/macscall.py"

# rule call_broadpeak_pe:
#     input:
#         lambda wildcards: [f"{{species}}/dedup_pe/{t}.bam" for t in [wildcards.sample, INPUT(wildcards.sample)] if not pd.isnull(t) and not t is None]
#     output:
#         expand("{{species}}/broadpeakcalling/{{sample}}_{filetype}", filetype=broad_output)
#     params:
#         extra="-f BAMPE"
#     script:
#         "scripts/macscall.py"

rule get_qvalues:
    input:
        "{macs_folder}/{sample}_treat_pileup.bdg",
	"{macs_folder}/{sample}_control_lambda.bdg"
    output:
        "{macs_folder}/{sample}_qvalues.bdg"
    shell:
        "macs2 bdgcmp -t {input[0]} -c {input[1]} -m qpois -o {output}"

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
        reference=os.path.join(config["data_dir"],"{species}/{species}.fa")
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
    script:
        "scripts/motifplot.py"

rule get_meme:
    output:
        "motives/{name}"
    shell:
        'wget {jaspar_address}{wildcards.name} -O {output} --user-agent="Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.8.1.6) Gecko/20070725 Firefox/2.0.0.6"'

rule get_domain_coverage:
    input:
        "{species}/clean_domains/{name}.bed",
        "{species}/data/chrom.sizes.txt"
    output:
        "{species}/domain_coverage/{name}.txt"
    shell:
        """
        awk '{{t+=($3-$2)}}END{{print t}}' {input[0]} > {output}
        awk '{{t+=$2}}END{{print t}}' {input[1]} >> {output}
        """

rule get_peak_coverage:
    input:
        "{species}/broadpeakcalling/{name}_peaks.broadPeak",
        "{species}/data/chrom.sizes.txt"
    output:
        "{species}/peak_coverage/{name}.txt"
    shell:
        """
        awk '{{t+=($3-$2)}}END{{print t}}' {input[0]} > {output}
        awk '{{t+=$2}}END{{print t}}' {input[1]} >> {output}
        """
