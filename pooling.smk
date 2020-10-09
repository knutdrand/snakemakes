pools = pd.read_csv(config["poolinginfo"]).set_index("pool")
pools_regexp = "|".join(pools.index).replace("+" ,"\+")
rule pool:
    input:
        lambda wildcards: [f"{{species}}/dedup/{sample}.bed.gz" for sample in
                           pools.loc[wildcards.stage]["samples"].split(";")]
    output:
        "{species}/dedup/{stage}_catpool.bed.gz"
    wildcard_constraints:
        stage=pools_regexp
    shell:
        "zcat {input} | gzip > {output}"
