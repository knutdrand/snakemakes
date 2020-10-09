import pandas as pd
import numpy as np
# include: "rules/common.smk"
sra_info = pd.read_csv(config["srainfo"]).set_index("sra")
#pools = set(sra_info["name"])
#
#get_sra_file_names = lambda wildcards: [f"sra_data/{sra}_{wildcards.read}.fastq.gz" 
#                                        for sra in sra_info[sra_info["name"]==wildcards.sample].index]
#gsm_info = pd.read_csv(config["gsminfo"]).set_index("GSM")
# gsms = [str(gsm).strip() for gsm in gsm_info.index if str(gsm).startswith("GSM")]
# sra_info = pd.read_csv("prototype_samples.csv")
# se_samples = sra_info[sra_info["LibraryLayout"]=="SINGLE"]["Run"]
# pe_samples = sra_info[sra_info["LibraryLayout"]=="PAIRED"]["Run"]
pe_samples=[]
se_samples=[]
gsms=[]
p12_samples = sra_info[sra_info["name"]=="mm10_P12_K4"].index
print(pe_samples)
print(se_samples)
pools = []
rule import_sra:
    output:
        "sra_data/{sra}_1.fastq.gz",
        "sra_data/{sra}_2.fastq.gz"
    shell:
        """
        fastq-dump --split-files --gzip {wildcards.sra} -O sra_data/
        """

# rule import_sra_pe:
#     output:
#         "pe_reads/{sra}_1.fastq.gz",
#         "pe_reads/{sra}_2.fastq.gz"
#     wildcard_constraints:
#         sra="|".join(pe_samples)
#     shell:
#         """
#         fastq-dump --split-files --gzip {wildcards.sra} -O sra_data/
#         """

rule prefetch:
    output:
        "{sra}/{sra}.sra"
    shell:
        "prefetch {wildcards.sra}"

# rule fastq_dump_se:
#     input:
#         "{sra}/{sra}.sra"
#     output:
#         "{sra}.fastq.gz"
#     shell:
#         "fastq-dump --gzip {input}"

rule fastq_dump_pe:
    input:
        "{sra}/{sra}.sra"
    output:
        "{sra}_1.fastq.gz",
        "{sra}_2.fastq.gz"
    shell:
        "fastq-dump --split-files --gzip {input}"

rule fastq_dump_se:
    input:
        "{sra}/{sra}.sra"
    output:
        "{sra}.fastq.gz",
    shell:
        "fastq-dump --gzip {input}"
# 
# 
# rule import_sra_se:
#     output:
#         "reads/{sra}.fastq.gz",
#     wildcard_constraints:
#         sra="|".join(se_samples)
#     shell:
#         """
#         fastq-dump --gzip {wildcards.sra} -O sra_data/
#         """

# rule merge_sra_data:
#     input:
#         get_sra_file_names
#     output:
#         "pe_reads/{sample}_{read}.fastq.gz",
#     wildcard_constraints:
#         read="1|2",
#         sample="|".join(pools)
#     shell:
#         "zcat {input} | gzip > {output}"

rule esearch:
    output:
        "infos/{sample}.csv"
    shell:
        "esearch -db sra -query {wildcards.sample} | efetch -format runinfo > {output}"

rule dump_all_paired:
    input:
        expand("{sample}_1.fastq.gz", sample=pe_samples)

rule dump_all_single:
    input:
        expand("{sample}.fastq.gz", sample=se_samples)


rule map_all_paired:
    input:
        expand("mm10/dedup_pe/{sample}.bed", sample=pe_samples)

rule map_all_p12:
    input:
        expand("mm10/mapped_pe/{sample}.bam", sample=p12_samples)

rule map_all_single:
    input:
        expand("mm10/dedup/{sample}.bed", sample=se_samples)


rule move_paired:
    input:
        "{sra}_1.fastq.gz",
        "{sra}_2.fastq.gz"
    output:
        "pe_reads/{sra}_1.fastq.gz",
        "pe_reads/{sra}_2.fastq.gz"
    shell:
        "mv {input} pe_reads/"

rule move_single:
    input:
        "{sra}.fastq.gz",
    output:
        "trimmed/{sra}.fastq.gz",
    shell:
        "mv {input} trimmed/"


rule dump_gsm_se:
    input:
        expand("infos/{sample}.csv", sample=gsms)
    output:
        "gsm_samples.csv"
    run:
        table = pd.concat([pd.read_csv(i).assign(GSM=gsm) for i, gsm in zip(input, gsms)]).set_index("GSM")
        table.to_csv(output[0])
                   
        
