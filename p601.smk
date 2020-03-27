import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re

compliments = {"A": "T", "C": "G", "G": "C", "T": "A"}
p601 = "CCGGATGCCCTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGTTCCAGTGCCGG"
reverse_p601 = p601.translate(compliments)[::-1]
analysis_info = pd.read_csv(config["analysisinfo"]).set_index("Name")

downsample_names = list(analysis_info[analysis_info["DownsamplePool"]==1].index)
all_names = ["acNuc", "acPreNuc", "27acNuc"]
rule cutadapt_p601:
    input:
        "reads/{sample}.fastq.gz"
    output:
        fastq="trimmed_601/{sample}.fastq.gz",
        qc="trimmed_601/{sample}.qc.txt"
    params:
        f'--nextseq-trim=20 -m 10 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a AATGATACGGCGACCACCGAGATCTACAC -a {p601} -a {reverse_p601}'
    log:
        "logs/cutadapt/{sample}.log"
    threads: 16
    wrapper:
        "0.50.0/bio/cutadapt/se"

rule cutp601_all:
    input:
        fastq=expand("trimmed_601/{sample}.fastq.gz", sample=downsample_names)

rule move_p601:
    input:
        "trimmed_601/{sample}.fastq.gz"
    output:
        "trimmed/{sample}_p601.fastq.gz"
    shell:
        "mv {input} {output}"

rule donwsample_to_lowest:
    #TODO: expand on pool
    input:
        fastqs=expand("reads/{sample}.fastq.gz", sample=downsample_names),
        counts=expand("reads/{sample}.fastq.gz.count", sample=downsample_names)
    output:
        expand("downsampled/{sample}.fastq", sample=downsample_names)
    script:
        "scripts/downsample.py"

rule move_downsampled:
    input:
        "downsampled/{sample}.fastq.gz"
    output:
        "reads/downsampled_{sample}.fastq.gz"
    shell:
        "mv {input} {output}"
    
rule all_p601:
    input:
        expand("trackhub/mm10/downsampled_{name}_p601.bw", name=downsample_names),
        expand("trackhub/mm10/{name}_p601.bw", name=all_names)

rule p601_trackhub:
    input:
        expand("trackhub/mm10/downsampled_{name}_p601.bw", name=downsample_names),
        expand("trackhub/mm10/{name}_p601.bw", name=all_names),
    output:
        "trackhub/{species}/trackDbp601.txt"
    shell:
        'chiptools trackdb single {input} > {output}'

rule p601fa:
    output:
        "data/p601.fa"
    run:
        open(output[0],"w").write(f">p601\n{p601}")

rule p601stats:
    input:
        "data/p601.fa.gz",
        "reads/downsampled_{sample}.fastq.gz"
    output:
        "p601stats/{sample}.npy"
    shell:
        "bwa mem {input} | chiptools alignscores 101 {output}"

rule p601summary:
    input:
        [f"p601stats/{sample}.npy" for sample in downsample_names]
    output:
        "p601summary.csv"
    run:
        hists = [np.load(i) for i in input]
        fractions = [np.sum(h[80:])/np.sum(h)*100 for h in hists]
        pd.DataFrame({"sample": downsample_names, "p601_prcnt": fractions}).to_csv(output[0], index=False)

rule p601plot:
    input:
        [f"p601stats/{sample}.npy" for sample in downsample_names]
    output:
        "p601fraction.png"
    run:
        import numpy as np
        hists = [np.load(i) for i in input]
        lines = [plt.plot(h/np.sum(h))[0] for h in hists]
        plt.legend(lines, downsample_names)
        plt.savefig(output[0])

rule trim_stats:
    input:
        qcs=expand("trimmed_601/{sample}.qc.txt", sample=downsample_names),
        counts=expand("reads/{sample}.fastq.gz.count", sample=downsample_names)
    output:
        "trim_stats.csv"
    run:
        lines = []
        for qc_file, count_file in zip(input.qcs, input.counts):
            text = open(qc_file).read()
            count= int(open(count_file).read().strip())
            f, r = [int(re.compile(f"{seq}.*Trimmed: (\\d+) times.*").findall(text)[0])
                    for seq in (p601, reverse_p601)]
            t = f+r
            p = t/count*100
            lines.append((count, f, r, t, p))
        with open(output[0], "w") as out_file:
            out_file.write(",".join(("sample", "read_count", "trimmedF", "trimmedR", "trimmed", "percent"))+"\n")
            for sample, (c, f, r, t, p) in zip(downsample_names, lines):
                out_file.write(f"{sample},{c},{f},{r},{t},{p}\n")
            
