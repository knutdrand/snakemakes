include: "rules/common.smk"

region_folders = expand("{flag}tss_containing_domains/sub{size}", flag=("non_", ""), size=(10000, 100000))
plottypes=["v", "heat", "tss", "average"]
wildcard_constraints:
    kind="treat_pileup|control_lambda|qvalues",
    plottype="|".join(plottypes),

rule regionplot:
    input:
        bedgraph="{species}/broadpeakcalling/{sample}_{kind}.bdg",
        regions="{species}/regions/{regionfolder}.bed"
    output:
        "{species}/regionplots/{regionfolder}/{sample}_{kind}_{plottype}.png",
        "{species}/regionplots/{regionfolder}/{sample}_{kind}_{plottype}.pkl"
    shell:
        "bdgplot {wildcards.plottype} {input.bedgraph} {input.regions} -o {output[0]} -od {output[1]}"

rule joinplots:
    input:
        expand_analysis("{species}/regionplots/{{regionfolder}}/{sample}_{{kind}}_{{plottype}}.pkl", "bdgplot")
    output:
        report("{species}/regionplots/{regionfolder}/all_{kind}_{plottype}.png", category="JoinedPlots")
    shell:
        "bdgtools joinfigs {wildcards.plottype} {input} -o {output} --name {wildcards.regionfolder}"


rule alldomainplots:
    input:
        expand("mm10/regionplots/{regionfolder}/2016_domains/all_treat_pileup_{plottype}.png",
               regionfolder=region_folders, plottype=["average", "v"]),
        [path for regionfolder in region_folders for path in  expand_analysis(f"{{species}}/regionplots/{regionfolder}/2016_domains/{{sample}}_treat_pileup_v.png", "bdgplot")]
        
