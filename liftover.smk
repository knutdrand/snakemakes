
def get_liftover_file(from_species, to_species):
    to_species = to_species.capitalize()
    return f"https://hgdownload.soe.ucsc.edu/goldenPath/{from_species}/liftOver/{from_species}To{to_species}.over.chain.gz"

rule download_liftover:
    output:
        config["data_dir"]+"{from_species}/{from_species}To{to_species}.over.chain.gz"
    run:
        link = get_liftover_file(wildcards.from_species, wildcards.to_species)
        shell("wget %s -O {output}" % link)

rule liftover_bdg:
    input:
        "{from_species}/{path}.bdg", 
        lambda wc: config["data_dir"]+f"{wc.from_species}/{wc.from_species}To{wc.to_species.capitalize()}.over.chain.gz"
    output:
        "{to_species}/from_{from_species}/{path}.bdg",
        "{to_species}/from_{from_species}/{path}.unmapped.bdg"
    shell:
        "liftOver {input} {output}"
