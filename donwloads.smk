include: "rules/common.smk"
annotations = {"susScr11": "ensembl"}

rule download_reference:
    output:
        config["data_dir"] + "{species}/{species}.fa.gz"
    shell:
        "wget http://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.species}/bigZips/{wildcards.species}.fa.gz -O {output}"

rule unzip_reference:
    input:
        config["data_dir"] + "{species}/{species}.fa.gz"
    output:
        config["data_dir"] + "{species}/{species}.fa"
    shell:
        "gunzip {input}"


rule download_genes:
    output:
        "{species}/data/refGene.txt.gz"
    shell:
        "wget https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.species}/database/refGene.txt.gz -O {output}"

rule get_genes_bed:
    input:
        lambda wildcards: "{species}/data/%sGene.txt.gz" % annotations.get(wildcards.species, "ref")
    output:
        "{species}/data/genes.bed"
    shell:
        """z%s {input} | awk '{{OFS="\t"}}{{print $3, $5, $6, ".", ".", $4}}' | uniq > {output}""" % chromosome_grep

rule download_chrom_sizes:
    output:
        "{species}/data/chromInfo.txt.gz"
    shell:
        "wget https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.species}/database/chromInfo.txt.gz -O {output}"

rule clean_chrom_sizes:
    input:
        "{species}/data/chromInfo.txt.gz"
    output:
        "{species}/data/chrom.sizes.txt"
    shell:
        f"z{chromosome_grep} {{input}} > {{output}}"
