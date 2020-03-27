library(ggplot2)
library(tidyverse)
print(snakemake@input[[1]])
snakemake@input[[1]] %>% read_csv %>% pivot_longer(-X1, names_to="type", values_to="count") -> data
ggplot() + geom_col(data=data, aes(x = X1, y = count, fill=type), position=position_dodge())
print(snakemake@output[[1]])
ggsave(snakemake@output[[1]], width=20)
