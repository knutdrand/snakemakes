library(tidyverse)
get_table = function(filename) filename %>% read_tsv(skip=1) %>% mutate(SampleName=basename(filename))
print(snakemake@output[1])
snakemake@input %>% lapply(get_table) %>% bind_rows -> table
print(table)
filename = snakemake@output[[1]]
print(filename)
write_tsv(table, filename)
