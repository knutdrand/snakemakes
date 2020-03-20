regiontypes = ["domains_minus_tss500", "non_tss_containing_domains", "domain_flanks",
                           "tss_containing_domains"]

rule get_gc:
    input:
        config["data_dir"] + "{species}/{species}.fa",
        "{species}/regions/{folder}/{name}.bed"
    output:
        "{species}/regions/{folder}/gc/{name}.txt"
    shell:
        "bedtools nuc -fi {input[0]} -bed {input[1]} | grep -v chrM | awk '{{gc+=$(NF-4)+$(NF-5);t+=$NF}}END{{print gc/t}}' > {output}"

rule summarize_gc:
    input:
        expand("{{species}}/regions/{regiontype}/gc/{{name}}.txt",
               regiontype=regiontypes
               )
    output:
        "{species}/gc_content/{name}.txt"
    shell:
        'paste {input} <(echo "{wildcards.name}") > {output}'

# rule all_gc:
#     input:
#         expand("hg38/v3/gc_content/{name}.txt",
#                name=pools)
#     output:
#         "hg38/v3/gc_content/ALL.tsv"
#     shell:
#         """
#         echo "domains_minus_tss500\tnon_tss_containing_domains\tdomain_flanks\ttss_containing_domains\tname" > {output}
#         cat {input} >> {output}
#         """
