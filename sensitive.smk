include: "rules/common.smk"

rule copy_non_sensitive:
    input:
        "{species}/dedup/{sample}.bed.gz"
    output:
        config["nonsensitive"] + "/{species}/dedup/{sample}.bed.gz"
    shell:
        "cp {input} {output}"

rule copy_info:
    input:
        "mapping_stats.csv"
    output:
        config["nonsensitive"] + "/mapping_stats_wo_alts.csv"
    shell:
        "cp {input} {output}"

rule all_non_sensitive:
    input:
        expand_analysis(config["nonsensitive"] + "/{species}/dedup/{sample}.bed.gz", "GB")
        
