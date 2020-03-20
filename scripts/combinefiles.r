library(tidyverse)
get_table = function(filename) filename %>% read_tsv(skip=1) %>% mutate(SampleName=filename)
snakemake@input %>% lapply(get_table) %>% bind_rows %>% write_tsv(snakemake@output)
