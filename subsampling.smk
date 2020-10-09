rule get_N_for_subsample:
    input:
        lambda wildcards: expand_analysis("{species}/dedup/{sample}.bed.gz.count", wildcards.analysis)
    output:
        "subsampling/{analysis}/count.txt"
    run:
        counts = [int(open(input_file).read().strip()) for input_file in input]
        for name, c in zip(input, counts):
            assert c>0, (name, c)
        open(output[0], "w").write(str(min(counts)))

rule subsample_from_count:
    input:
        "{species}/dedup/{sample}.bed.gz",
        "{species}/dedup/{sample}.bed.gz.count",
        "subsampling/{analysis}/count.txt"
    output:
        "{species}/dedup/{analysis}_subsampled_{sample}.bed.gz"
    run:
        N = int(open(input[2]).read().strip())
        my_N = int(open(input[1]).read().strip())
        shell("zcat {input[0]} | awk '{{if (rand()*%s<%s)print}}' |  gzip > {output}" % (my_N, N))

rule all_subsample_coverage:
    input:
        lambda wildcards: expand_analysis("{species}/{{region}}_coverage/{{analysis}}_subsampled_{sample}.txt", wildcards.analysis)
    output:
        "{region}_coverage_subsampled_{analysis}.csv"
    run:
        with open(output[0], "w") as out_file:
            out_file.write(",".join(("Sample", "DomainCoverage", "Genomesize", "Prct"))+"\n")
            for in_file in input:
                c, g = open(in_file).read().strip().split()
                out_file.write(",".join([Path(in_file).stem, c, g, str(int(c)/int(g))])+"\n")
