rule size_hist:
    input:
        "{path}/{name}.bed"
    output:
        "{path}_logsize_hist/{name}.npz"
    wildcard_constraints:
        name="[^/]+"
    script:
        "scripts/size_hist.py"

rule adj_size_hist:
    input:
        "{path}/{name}.bed"
    output:
        "{path}_adjsize_hist/{name}.npz"
    wildcard_constraints:
        name="[^/]+"
    script:
        "scripts/adj_size_hist.py"

rule peak_size_hist:
    input:
        "{path}/{name}_peaks.broadPeak"
    output:
        "{path}_hist_peaks/{name}.npz"
    wildcard_constraints:
        name="[^/]+"
    script:
        "scripts/size_hist.py"

rule gapped_peak_size_hist:
    input:
        "{path}/{name}_peaks.gappedPeak"
    output:
        "{path}_hist_gapped_peaks/{name}.npz"
    wildcard_constraints:
        name="[^/]+"
    params:
        n_bins=100
    script:
        "scripts/size_hist.py"

rule heatplot:
    input:
        bedgraph="{root}/broadpeakcalling/{sample}_{kind}.bdg",
        domains="{root}/clean_domains/{sample}.bed"
    output:
        multiext("{root}/heatplots/{sample}_{kind}", ".pkl", ".png")
    wildcard_constraints:
        kind="treat_pileup|control_lambda|qvalues"
    shell:
        "chipplots heat {input.bedgraph} {input.domains} -od {output[0]} -o {output[1]}"
        #"cat {input.bedgraph} | chiptools heatplot {input.domains} {output}"

rule combined_heatplot:
    input:
        bedgraph="{root}/broadpeakcalling/{sample}_{kind}.bdg",
        domains="{root}/clean_domains/{domainsample}.bed"
    output:
        multiext("{root}/heatplots/{domainsample}/{sample}_{kind}", ".pkl", ".png")
    shell:
        "chipplots heat {input.bedgraph} {input.domains} -od {output[0]} -o {output[1]}"

rule region_averageplots:
    input:
        bedgraph="{species}/broadpeakcalling/{sample}_treat_pileup.bdg",
        domains="{species}/regions/{folder}/sub{size}/{other_sample}.bed"
    output:
        multiext("{species}/average_plots/regions/{folder}/sub{size}/{other_sample}/{sample}", ".pkl", ".png")
    shell:
        "chipplots average {input.bedgraph} {input.domains} -od {output[0]} -o {output[1]}"

rule all_region_average_plots:
    input:
        expand("hg38/average_plots/regions/tss_containing_domains/sub10000/Oocyte_K4_pool_II/{sample}.pkl",
               sample=get_samples_for_analysis("GB"))
    output:
        "hg38/combined_average_plots/regions/tss_containing_domains/sub10000/Oocyte_K4_pool_II/all.png"
    params:
        samples=get_samples_for_analysis("GB")
    shell:
        "plotjoin {input} -o {output}"


rule all_region_average_plots2:
    input:
        expand("hg38/average_plots/regions/tss_containing_domains/sub100000/Oocyte_K4_pool_II/{sample}.pkl",
               sample=get_samples_for_analysis("GB"))
    output:
        "hg38/combined_average_plots/regions/tss_containing_domains/sub100000/Oocyte_K4_pool_II/all.png"
    params:
        samples=get_samples_for_analysis("GB")
    shell:
        "plotjoin {input} -o {output}"

rule all_region_average_plots3:
    input:
        expand("hg38/average_plots/regions/non_tss_containing_domains/sub10000/Oocyte_K4_pool_II/{sample}.pkl",
               sample=get_samples_for_analysis("GB"))
    output:
        "hg38/combined_average_plots/regions/non_tss_containing_domains/sub10000/Oocyte_K4_pool_II/all.png"
    params:
        samples=get_samples_for_analysis("GB")
    shell:
        "plotjoin {input} -o {output}"

rule all_region_average_plots4:
    input:
        expand("hg38/average_plots/regions/non_tss_containing_domains/sub100000/Oocyte_K4_pool_II/{sample}.pkl",
               sample=get_samples_for_analysis("GB"))
    output:
        "hg38/combined_average_plots/regions/non_tss_containing_domains/sub100000/Oocyte_K4_pool_II/all.png"
    params:
        samples=get_samples_for_analysis("GB")
    shell:
        "plotjoin {input} -o {output}"


rule gene_vplots:
    input:
        "{species}/broadpeakcalling/{sample}_treat_pileup.bdg",
        "{species}/data/genes.bed"
    output:
        "{species}/vplots/genes/{sample}.pkl",
        "{species}/vplots/genes/{sample}.png",
    shell:
        "chipplots v {input[0]} {input[1]} -od {output[0]} -o {output[1]}"

rule region_matrixplots:
    input:
        "{species}/broadpeakcalling/{sample}_treat_pileup.bdg",
        "{species}/regions/{folder}/sub{size}/{other_sample}.bed"
    output:
        "{species}/{plottype}plots/regions/{folder}/sub{size}/{other_sample}/{sample}.pkl",
        "{species}/{plottype}plots/regions/{folder}/sub{size}/{other_sample}/{sample}.png",
    wildcard_constraints:
        other_sample="[^/]+",
        plottype="v|heat"
    shell:
        "chipplots {wildcards.plottype} {input[0]} {input[1]} -od {output[0]} -o {output[1]} --regionsize {wildcards.size}"
        # "chiptools vplot {input[1]} {input[0]} {output} --undirected -s {wildcards.size}"

rule all_region_vplots:
    input:
        expand("hg38/{plottype}plots/regions/{pred}tss_containing_domains/sub{size}/Oocyte_K4_pool_II/{sample}.png",
               sample=get_samples_for_analysis("GB"), pred=["non_", ""], size=[10000, 100000], plottype=["v", "heat"])

rule logvplot:
    input:
        "{species}/vplots/{place}/{name}.npy",
    output:
        report("{species}/logvplots/{place}/{name}.png", category="VPlots")
    shell:
        pyplot("plt.imshow(np.log(np.load('{input}')+1));plt.savefig('{output}')")

rule all_gene_logvplot:
    input:
        expand_analysis("{species}/logvplots/genes/{sample}.png", "vplot")
    
        
rule binned_reads:
    input:
        windows="{species}/regions/{folder}.bed",
        reads="{species}/dedup/{sample}.bed.gz"
    output:
        "{species}/region_counts/{folder}/{sample}.bed"
    shell:
        "z{chromosome_grep} {input.reads} | bedtools coverage -a {input.windows} -b stdin -counts > {output}"

rule domain_reads:
    input:
        windows="{species}/clean_domains/{domains}.bed",
        reads="{species}/dedup/{sample}.bed.gz"
    output:
        "{species}/domain_counts/{domains}/{sample}.bed"
    shell:
        "z{chromosome_grep} {input.reads} | bedtools coverage -a {input.windows} -b stdin -counts > {output}"

# rule binned_correlation:
#     input:
#         "{species}/binned_counts/{windowsize}/{sample_a}.bed",
#         "{species}/binned_counts/{windowsize}/{sample_b}.bed"
#     output:
#         "{species}/correlations/{windowsize}/{sample_a}_VS_{sample_b}.png"
#     script:
#         "scripts/correlation.py"

rule region_correlation:
    input:
        "{species}/region_counts/{folder}/{sample_a}.bed",
        "{species}/region_counts/{folder}/{sample_b}.bed"
    output:
        "{species}/correlations/{folder}/{sample_a}_VS_{sample_b}.png"
    script:
        "scripts/correlation.py"

rule check_domain_enrichment:
    input:
        domains="{species}/domain_counts/{domains}/{sample}.bed",
        count="{species}/dedup/{sample}.bed.gz.count",
        genome_size="{species}/data/genome_size.txt"
    output:
        "{species}/domain_counts/{domains}/{sample}.csv"
    run:
        count = int(open(input.count).read().strip())
        genome_size = int(open(input.genome_size).read().strip())
        parts = (line.split() for line in open(input.domains))
        sizes, counts = zip(*((int(p[2])-int(p[1]), int(p[-1])) for p in parts))
        t_size = sum(sizes)
        t_counts = sum(counts)
        open(output[0], "w").write("%s,%s,%s,%s" % (t_counts, t_size, count-t_counts, genome_size-t_size))

rule all_correlation:
    input:
        get_all_corr
