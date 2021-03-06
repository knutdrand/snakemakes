from snakemake.shell import shell
import pandas as pd
include: "rules/common.smk"

rule bamtobedpe2se:
    input:
        "{species}/dedup_pe/{sample}.bam",
    output:
        "{species}/dedup/{sample}_r1.bed.gz",
    shell:
        "samtools view -b -f 64 {input} | bedtools bamtobed -i - | %s | gzip > {output}" % chromosome_grep

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
    script:
        "scripts/cutadaptse.py"
    #wrapper:
    #    "0.50.0/bio/cutadapt/se"

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
        sort="samtools",
    threads: 16
    script:
        "scripts/bwamem.py"
    #wrapper:
    #    "0.49.0/bio/bwa/mem"

rule bwa_mem_pe:
    input:
        reads=expand("pe_reads/{{sample}}_{read}.fastq.gz", read=[1,2])
    output:
        "{species}/mapped_pe/{sample}.bam"
    log:
        "logs/bwa_mem/{species}/{sample}.log"
    params:
        index=config["data_dir"]+"{species}/{species}.fa.gz",
        # sort="samtools",
        # sort_order="queryname"
    wildcard_constraints:
        species="[^/]+"
    threads:
        16
    wrapper:
        "0.49.0/bio/bwa/mem"
#        (bwa mem -t {threads} {params.index} {input.reads} | samtools sort -o {output}) 2> {log}

rule filter:
    input:
        "{folder}mapped/{sample}.bam"
    output:
        "{folder}mapped_filtered/{sample}.bam"
    params:
        "-Bb -q 30" 
    wildcard_constraints:
        sample="[^/]+"
    script:
        "scripts/samtoolsview.py"
#    wrapper:
#        "0.50.4/bio/samtools/view"

rule filter_pe:
    input:
        "{species}/{folder}_pe/{sample}.bam"
    output:
        "{species}/{folder}_pe_filtered/{sample}.bam"
    params:
        "-Bb -q 30 -F 1804 -f 2"
    wildcard_constraints:
        sample="[^/]+",
    shell:
        "samtools view {params} {input} > {output}"

rule filter_pe_bad:
    input:
        "{species}/{folder}_pe/{sample}.bam"
    output:
        "{species}/{folder}_pe_filtered/unfiltered_{sample}.bam"
    params:
        "-Bb -F 1804 -f 2"
    wildcard_constraints:
        sample="[^/]+",
    shell:
        "samtools sort -n {input} -T {wildcards.sample}.tmp | samtools view {params} > {output}"

#     wrapper:
#         "0.50.4/bio/samtools/view"


rule samtools_remove_duplicates:
    input:
        "{species}/mapped_filtered/{sample}.bam"
    output:
        bam="{species}/dedup/{sample}.bam",
    log:
        "logs/dedup/{species}/{sample}.log"
    shell:
        "samtools markdup -rs {input} {output} 2> {log}"

rule fixmate:
    input:
        "{species}/mapped_pe_filtered/{sample}.bam"
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
        "{species}/dedup/{sample}.bam",
    output:
        "{species}/dedup/{sample}.bed",
    shell:
        "bedtools bamtobed -i {input} > {output}"

rule collatebam:
    input:
        "{species}/dedup_pe/{sample}.bam",
    output:
        "{species}/dedup_pe/{sample}.collated.bam",
    shell:
        "samtools collate {input} -o {output}"

rule bamtobedpe:
    input:
        "{species}/dedup_pe/{sample}.collated.bam",
    output:
        "{species}/dedup_pe/{sample}.bed",
    shell:
        "bedtools bamtobed -i {input} | chiptools pairbed | bedtools sort > {output}"


rule petoeaseq:
    input:
        "{species}/dedup_pe/{sample}.bed"
    output:
        "{species}/dedup_pe/{sample}.easeq.bed"
    shell:
        chromosome_grep + """ {input} | awk '{{OFS="\t"}}{{print $1,$2,$3,".",".","+"; print $1,$2,$3,".",".","-"}}' > {output}"""

rule fastq_screen:
    input:
        "{folder}/{sample}.fastq.gz"
    output:
        txt="qc/fastq_screen/{folder}/{sample}.txt",
        png=report("qc/fastq_screen/{folder}/{sample}.png", category="SpeciesScreen")
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
        
rule fragment_size_hist:
    input:
        "{species}/dedup_pe/{sample}.bed"
    output:
        data="{species}/pe_frag_size/{sample}.npy",
        fig="{species}/pe_frag_size/{sample}.png"
    script:
        "scripts/fragment_size_hist.py"

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

rule clean_dedup:
    input:
        "{folder}/dedup/{sample}.bed.gz"
    output:
        "{folder}/clean_dedup/{sample}.bed.gz"
    shell:
        "z%s {input} | gzip > {output}" % chromosome_grep

rule subsample:
    input:
        "{species}/dedup/{sample}.bed.gz"
    output:
        "{species}/dedup/{sample}_{n_reads}reads.bed.gz"
    wildcard_constraints:
        n_reads = "\d*"
    shell:
        "zcat {input} | shuf -n {wildcards.n_reads} | gzip > {output}"

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
