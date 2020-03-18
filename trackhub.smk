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

rule create_bw_track:
    input:
        "src/bdg2bw",
        "{species}/{name}.bdg",
        "{species}/data/chrom.sizes.txt"
    output:
        "{species}/{name}.bw"
    wildcard_constraints:
        species="[^/]+"
    shell:
        "{input}"

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
