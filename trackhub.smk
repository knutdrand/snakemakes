include: "rules/common.smk"
track_types = ["domains.bb", "peaks.bb", "treat_pileup.bw", "control_lambda.bw", "qvalues.bw"]

rule get_bdg2bw:
    output:
        "src/bdg2bw"
    shell:
        """
        wget https://gist.githubusercontent.com/taoliu/2469050/raw/34476e91ebd3ba9d26345da88fd2a4e7b893deea/bdg2bw -O {output}
        chmod a+x {output}
        """

rule move_to_trackhub:
    input:
        "{species}/broadpeakcalling/{name}_treat_pileup.bw"
    output:
        "trackhub/{species}/{name}.bw"
    wildcard_constraints:
        species="[^/]+"
    shell:
        "mv {input} {output}"

rule trackhub:
    input:
        lambda wildcards: expand("trackhub/{{species}}/{sample}.bw", sample=get_species_tracks(wildcards.species))
    output:
        "trackhub/{species}/trackDb.txt"
    shell:
        'chiptools trackdb single {input} > {output}'

# rule adv_trackhub:
#     input:
#         lambda wildcards: expand("trackhub/{{species}}/{sample}_{filetype}",
#                                  sample=get_species_tracks(wildcards.species),
#                                  filetype=track_types)
#     output:
#         "trackhub/{species}/trackDb.txt"
#     shell:
#         'chiptools trackdb {input} > {output}'

rule clip_bw:
    input:
        "{species}/{name}.bdg",
        "{species}/data/chrom.sizes.txt"
    output:
        "{species}/{name}.bdg.clip"
    wildcard_constraints:
        species="[^/]+"
    shell:
        "chiptools clipbed {input} > {output}"

rule ucsc_sort:
    input:
        "{species}/{name}.bdg.clip"
    output:
        "{species}/{name}.bdg.clip.uscssort"
    wildcard_constraints:
        species="[^/]+"
    shell:
        "LC_COLLATE=C sort -k1,1 -k2,2n {input} -T tmp/ > {output}"


rule remove_overlap:
    input:
        "{fullpath}.bdg.clip.uscssort"
    output:
        "{fullpath}.bdg.clip.uscssort.disjoint"
    shell:
        "bedRemoveOverlap {input} {output}"

rule create_bw_track:
    input:
        bedGraph="{species}/{name}.bdg.clip.uscssort.disjoint",
        chromsizes="{species}/data/chrom.sizes.txt"
    output:
        "{species}/{name}.bw"
    wildcard_constraints:
        species="[^/]+"
    shell:
        "bedGraphToBigWig {input} {output}"
#        "0.50.3/bio/ucsc/bedGraphToBigWig"
