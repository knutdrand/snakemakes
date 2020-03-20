chromosome_grep = "grep -Ew -e 'chr[0-9]{{1,2}}' -e chrX -e chrY"

rule download_genes:
    output:
        "{species}/data/refGene.txt.gz"
    shell:
        "wget https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.species}/database/refGene.txt.gz -O {output}"

rule get_genes_bed:
    input:
        "{species}/data/refGene.txt.gz"
    output:
        "{species}/data/genes.bed"
    shell:
        """z%s {input} | awk '{{OFS="\t"}}{{print $3, $5, $6, ".", ".", $4}}' > {output}""" % chromosome_grep

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
        "z"+chromosome_grep + " {input} > {output}"
