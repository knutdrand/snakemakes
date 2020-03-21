from snakemake.shell import shell
import pandas as pd

samples = pd.read_csv(config["sampleinfo"]).set_index("filename")

def LANE_FILES(name):
    return sorted([i for i, e in  samples.loc[samples["name"]==name].iterrows()])

rule merge_se:
    input:
        lambda wildcards: expand("raw/{sample}.fastq.gz", sample=LANE_FILES(wildcards.name))
    output:
        "reads/{name}.fastq.gz"
    shell:
        "zcat {input} | gzip > {output}"

rule cutadapt_se:
    input:
        "reads/{sample}.fastq.gz"
    output:
        fastq="trimmed/{sample}.fastq.gz",
        qc="trimmed/{sample}.qc.txt"
    params:
        '--nextseq-trim=20 -m 10 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a AATGATACGGCGACCACCGAGATCTACAC'
    log:
        "logs/cutadapt/{sample}.log"
    threads: 16 
    wrapper:
        "0.50.0/bio/cutadapt/se"
#         cutadapt {params} -j {threads} -o {output.fastq} {input[0]} > {output.qc} 2> {log}


rule bwa_index:
    input:
        "{name}.fa.gz"
    output:
        multiext("{name}.fa.gz", ".ann", ".amb", ".bwt", ".sa", ".pac")
    shell:
        "bwa index {input}"

rule bwa_mem_se:
    input:
        reads="trimmed/{sample}.fastq.gz"
    output:
        "{species}/mapped/{sample}.bam"
    log:
        "logs/bwa_mem/{species}/{sample}.log"
    params:
        index=config["data_dir"]+"{species}/{species}.fa.gz",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
    threads: 16
    wrapper:
        "0.49.0/bio/bwa/mem"
#        (bwa mem -t {threads} {params.index} {input.reads} | samtools sort -o {output}) 2> {log}

rule filter:
    input:
        "{folder}/{sample}.bam"
    output:
        "{folder}_filtered/{sample}.bam"
    params:
        "-Bb -q 30" 
    wildcard_constraints:
        sample="[^/]+"
    wrapper:
        "0.50.0/bio/samtools/view"

rule samtools_remove_duplicates:
    input:
        "{species}/mapped_filtered/{sample}.bam"
    output:
        bam="{species}/dedup/{sample}.bam",
    log:
        "logs/dedup/{species}/{sample}.log"
    shell:
        "samtools markdup -rs {input} {output} 2> {log}"

rule bamtobed:
    input:
        "{name}.bam"
    output:
        "{name}.bed"
    shell:
        "bedtools bamtobed -i {input} > {output}"

rule fastq_screen:
    input:
        "reads/{sample}.fastq.gz"
    output:
        txt="qc/fastq_screen/{sample}.txt",
        png=report("qc/fastq_screen/{sample}.png", category="SpeciesScreen")
    params:
        fastq_screen_config=config["data_dir"]+"fastq_screen_config.txt",
        subset=1000000,
        aligner='bwa'
    threads: 8
    wrapper:
        "0.50.0/bio/fastq_screen"

#rule mapping_stats:
#    input:
#        [f"{{species}}/{folder}/{{sample}}.{ext}" for folder, ext in [("reads", "fastq.gz"), ("trimmed", "fastq.gz"), ("mapped_filtered", "bam"), ("dedup", "bam")]]
#    output:
#        "{species}/mappingstats/{sample.txt}"
#    shell:
#        """
#        touch {output}
#        zgrep ^@ {input[0]} | wc -l > {output}
#        zgrep ^@ {input[1]} | wc -l > {output}
#        samtools view {input[2]} | wc -l > {output}
#        samtools view {input[3]} | wc -l > {output}
#        """
