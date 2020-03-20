track_types = ["domains.bb", "peaks.bb", "treat_pileup.bw", "control_lambda.bw", "qvalues.bw"]
analysis_info = pd.read_csv(config["analysisinfo"]).set_index("Name")
chromosome_grep = "grep -Ew -e 'chr[0-9]{{1,2}}' -e chrX -e chrY"

def get_species_tracks(species):
    species_dict = {"mm10": "mouse", "bosTau8": "bovine"}
    d = analysis_info[analysis_info["Species"]==species_dict[species]]
    d = d[d["GB"]==1]
    return d.index

rule get_bdg2bw:
    output:
        "src/bdg2bw"
    shell:
        """
        wget https://gist.githubusercontent.com/taoliu/2469050/raw/34476e91ebd3ba9d26345da88fd2a4e7b893deea/bdg2bw -O {output}
        chmod a+x {output}
        """

# rule create_bw_track:
#     input:
#         "src/bdg2bw",
#         "{species}/{name}.bdg",
#         "{species}/data/chrom.sizes.txt"
#     output:
#         temp("{species}/{name}.clean.bdg"),
#         "{species}/{name}.clean.bw"
#     wildcard_constraints:
#         species="[^/]+"
#     shell:
#         """
#         %s {input[1]} > {output[0]}
#         {input[0]} {output[0]} {input[2]}
#         """ % chromosome_grep



rule move_to_trackhub:
    input:
        "{species}/peakcalling/{name}_treat_pileup.bw"
    output:
        "trackhub/{species}/{name}.bw"
    wildcard_constraints:
        species="[^/]+"
    shell:
        "mv {input} {output}"

rule trackhub:
    input:
        lambda wildcards: expand("trackhub/{{species}}/{name}.bw", name=get_species_tracks(wildcards.species))
    output:
        "trackhub/{species}/trackDb.txt"
    shell:
        'chiptools trackdb single {input} > {output}'

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
     
# "%s {input[0]} | bedtools slop -i - -g {input[1]} -b 0 | bedClip stdin {input[1]} {output}" % chromosome_grep

rule ucsc_sort:
    input:
        "{species}/{name}.bdg.clip"
    output:
        "{species}/{name}.bdg.clip.uscssort"
    wildcard_constraints:
        species="[^/]+"
    shell:
        "LC_COLLATE=C sort -k1,1 -k2,2n {input} -T tmp/ > {output}"

rule create_bw_track:
    input:
        bedGraph="{species}/{name}.bdg.clip.uscssort",
        chromsizes="{species}/data/chrom.sizes.txt"
    output:
        "{species}/{name}.bw"
    wildcard_constraints:
        species="[^/]+"
    shell:
        "bedGraphToBigWig {input} {output}"
#        "0.50.3/bio/ucsc/bedGraphToBigWig"
