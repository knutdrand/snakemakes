names = ["Species"]+ [i.split("/")[-1].split(".")[0] for i in snakemake.input]
humans_row = ["Human"] 
mouse_row = ["Mouse"]
for i in snakemake.input:
    parts = [line.split() for line in open(i)]
    values = {p[0]: p[5] for p in parts  if p and p[0] in ("Mouse", "Human")}
    humans_row.append(values["Human"])
    mouse_row.append(values["Mouse"])
with open(snakemake.output[0], "w") as f:
    for line in [names, humans_row, mouse_row]:
        f.write("\t".join(line)+"\n")
