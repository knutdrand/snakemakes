library(ggplot2)
library(tidyverse)
snakemake@input %>% read_csv %>% ggplot() + geom_col(data = new_data, aes(x = X1, y = count, fill=type), position=position_dodge())
ggsave(snakemake@output, width=20)
