sra_info = pd.read_csv(config["srainfo"]).set_index("sra")
pools = set(sra_info["name"])

def get_sra_file_names(wildcards):
    sra_idxs = sra_info[sra_info["name"]==wildcards.sample].index
    return [f"sra_data/{sra}_{wildcards.read}.fastq.gz" for sra in sra_idxs]

rule import_sra:
    output:
        "sra_data/{sra}_1.fastq.gz",
        "sra_data/{sra}_2.fastq.gz"
    shell:
        """
        fastq-dump --split-files --gzip {wildcards.sra} -O sra_data/
        """
rule merge_sra_data:
    input:
        get_sra_file_names
    output:
        "pe_reads/{sample}_{read}.fastq.gz",
    wildcard_constraints:
        read="1|2",
        sample="|".join(pools)
    shell:
        "zcat {input} | gzip > {output}"
