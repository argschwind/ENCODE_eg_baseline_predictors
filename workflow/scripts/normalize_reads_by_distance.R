## Normalize reads by distance predictor to sum up to 1 for each gene

# required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

# column types in reads by distance prediction file
pred_cols <- cols(
  .default = col_double(),
  chr = col_character(),
  name = col_character(),
  class = col_character(),
  TargetGene = col_character(),
  TargetGeneEnsemblID = col_character(),
  CellType = col_character()
)

# load prediction file
pred <- read_tsv(snakemake@input[[1]], col_types = pred_cols)

# compute normalized prediction scores
pred <- pred %>% 
  group_by(TargetGene, TargetGeneTSS) %>% 
  mutate(Score = Score / sum(Score)) %>% 
  replace_na(replace = list(Score = 0))

# write to output file
write_tsv(pred, file = snakemake@output[[1]])
