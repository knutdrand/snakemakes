chromosome_grep = "grep -Ew -e 'chr[0-9]{{1,2}}' -e chrX -e chrY"

rule get_clean_domains:
    input:
        "{species}/domains/{name}.bed",
    output:
        "{species}/clean_domains/{name}.bed"
    shell:
        chromosome_grep + " {input}>{output}"

rule get_tss:
    input:
        "{species}/data/genes.bed"
    output:
        "{species}/regions/tss/1.bed"
    shell:
        """awk '{{OFS="\t"}}{{if ($6=="+") {{print $1,$2,$2+1}} else {{print $1, $3-1, $3}}}}' {input} | uniq > {output}"""

rule get_tss_500:
    input:
        "{species}/regions/tss/1.bed",
        "{species}/data/chrom.sizes.txt"
    output:
        "{species}/regions/tss/500.bed"
    shell:
        "bedtools slop -i {input[0]} -g {input[1]} -b 500 > {output}"

rule get_domain_flanks:
    input:
        "{species}/clean_domains/{name}.bed",
        "{species}/data/chrom.sizes.txt"
    output:
        "{species}/regions/domain_flanks/{name}.bed"
    shell:
        "bedtools flank -i {input[0]} -g {input[1]}  -b 5000 | awk '{{if ($3>$2) print}}'  > {output}"

rule get_tss_containing_domains:
    input:
        "{species}/clean_domains/{name}.bed",
        "{species}/regions/tss/1.bed"
    output:
        "{species}/regions/tss_containing_domains/{name}.bed"
    shell:
        """bedtools intersect -a {input[0]} -b {input[1]} -c | awk '{{if ($NF>0) print}}' > {output}"""

rule get_non_tss_containing_domains:
    input:
        "{species}/clean_domains/{name}.bed",
        "{species}/regions/tss/1.bed"
    output:
        "{species}/regions/non_tss_containing_domains/{name}.bed"
    shell:
        "bedtools intersect -a {input[0]} -b {input[1]} -c | awk '{{if ($NF==0) print}}' > {output}"

rule get_tss_containing_domains_minus_tss500:
    input:
        "{species}/regions/tss_containing_domains/{name}.bed",
        "{species}/regions/tss/500.bed"
    output:
        "{species}/regions/domains_minus_tss500/{name}.bed"
    shell:
        "bedtools subtract -a {input[0]} -b {input[1]} > {output}"

