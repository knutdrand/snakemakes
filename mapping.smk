from snakemake.shell import shell
import pandas as pd

samples = pd.read_csv(config["sampleinfo"]).set_index("filename")
analysis_info = pd.read_csv(config["analysisinfo"]).set_index("Name")
def LANE_FILES(name):
    return sorted([i for i, e in  samples.loc[samples["name"]==name].iterrows()])

# rule merge_se:
#     input:
#         lambda wildcards: expand("raw/{sample}.fastq.gz", sample=LANE_FILES(wildcards.name))
#     output:
#         "reads/{name}.fastq.gz"
#     shell:
#         "zcat {input} | gzip > {output}"

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

rule bwa_mem_pe:
    input:
        reads=expand("merged_sra_reads/{{sample}}_{read}.fastq.gz", read=[1,2])
    output:
        "{species}/mapped_pe/{sample}.bam"
    log:
        "logs/bwa_mem/{species}/{sample}.log"
    params:
        index=config["data_dir"]+"{species}/{species}.fa.gz",
        sort="samtools",
        sort_order="queryname"
        # Can be 'none', 'samtools' or 'picard'.
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

rule filter_pe:
    input:
        "{folder}/{sample}.bam"
    output:
        "{folder}_pe_filtered/{sample}.bam"
    params:
        "-Bb -q 30 -F 1804 -f 2 {input} > {output[1]}" 
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

rule tmp_sort:
    input:
        "{species}/mapped_pe_filtered/{sample}.bam"
    output:
        "{species}/mapped_pe_filtered_namesort/{sample}.bam"
    shell:
        "samtools sort -n {input} -o  {output}"

rule fixmate:
    input:
        "{species}/mapped_pe_filtered_namesort/{sample}.bam"
    output:
        bam="{species}/mapped_pe_matefix/{sample}.bam"
    shell:
        "samtools fixmate --threads {threads} -mr {input} {output}"

rule samtools_remove_duplicates_pe:
    input:
        "{species}/mapped_pe_matefix/{sample}.bam"
    output:
        bam="{species}/dedup_pe/{sample}.bam",
    log:
        "logs/dedup/{species}/{sample}.log"
    threads:
        4
    shell:
        "samtools sort {input} -o - | samtools markdup -rs --threads {threads} - {output} 2> {log}"

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

rule count_fastq:
    input:
        "{name}.fastq.gz"
    output:
        "{name}.fastq.gz.count"
    shell:
        "zcat {input} | awk '{{s++}}END{{print s/4}}' > {output}"

rule count_bam:
    input:
        "{name}.bam"
    output:
        "{name}.bam.count"
    shell:
        "samtools view -c {input} > {output}"

rule fastq_size_hist:
    input:
        "{folder}/{sample}.fastq.gz"
    output:
        "{folder}_sizehist/{sample}.txt"
    params:
        nbins=102,
        bin_size=1,
        do_log=False
    wildcard_constraints:
        sample="[^/]+"
    shell:
        "zcat {input} | awk '{{if (NR % 4 == 2) ++a[length()]}} END{{for (i in a) print i, a[i]}}' > {output}"
#         "scripts/size_hist.py"

rule trimming_effect_plot:
    input:
        "reads_sizehist/{sample}.txt",
        "trimmed_sizehist/{sample}.txt"
    output:
        report("trimming_effects/{sample}.png", category="Trimming hist")
    run:
        get_values = lambda f: [(int(l), int(count)) for l, count in [line.split() for line in open(f)]]
        values = [get_values(f) for f in input]
        size = max(max(v) for v in values)
        lines = []
        
        for i, (f, lable) in enumerate(zip(input, ("reads", "trimmed"))):
            l, count = zip(*get_values(f))
            plt.bar([p+i*0.5 for p in l], list(count), label=lable)
            #a = np.zeros(size+1, dtype="int")
            #for l, count in get_values(f).items():
            #    a[l] = count
            
            #lines.append(plt.plot(a)[0])
        #plt.legend(lines, ("reads", "trimmed"))
        plt.savefig(output[0])

rule fastqc_all:
    input:
        expand("qc/fastqc/{sample}.html", sample=analysis_info.index)
        

rule fastqc:
    input:
        "reads/{sample}.fastq.gz"
    output:
        html="qc/fastqc/{sample}.html",
        zip="qc/fastqc/{sample}_fastqc.zip" 
    params: ""
    log:
        "logs/fastqc/{sample}.log"
    wrapper:
        "0.50.4/bio/fastqc"


#rule mapping_tats:
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
