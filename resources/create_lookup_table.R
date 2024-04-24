## Add ABC prediction files to lookup table

# required packages
library(tidyverse)
library(here)

# lookup table for baseline predictors
lookup <- read_tsv(here("resources/ENCODE_lookup.combined.unique.scratch.tsv"))

# DNase + H3K27ac ABC samples
abc <- read_tsv(read_tsv("resources/biosample_old_new_lookup.tsv", show_col_types = FALSE))

# only retain one abc sample per sample...
abc <- abc %>% 
  group_by(sample) %>% 
  slice_head(n = 1)

# add ABC samples to lookup table
lookup <- left_join(lookup, select(abc, sample, ABC_sample = Biosample), by = "sample")

# directories containing ABC files
abc_pred_dir <- "/oak/stanford/groups/akundaje/kmualim/02212022_hg38Preds_QNORM-avg_hic_track2"
abc_peak_dir <- "/oak/stanford/groups/akundaje/kmualim/02212022_hg38Preds_QNORM-avg_hic_track1"

# infer paths to ABC files for each sample that could be assign to an ABC sample
output <- lookup %>%
  rowwise() %>% 
  mutate(ABC_preds = if_else(!is.na(ABC_sample),
                             true = file.path(abc_pred_dir, ABC_sample, "Predictions", "EnhancerPredictionsAllPutative.txt.gz"),
                             false = "none"),
         ABC_peaks = if_else(!is.na(ABC_sample),
                             true = file.path(abc_peak_dir, ABC_sample, "Peaks", "macs2_peaks.narrowPeak.sorted.candidateRegions.bed"),
                             false = "none")
  )

# reformat for output
output <- replace_na(output, replace = list(ABC_sample = "none"))

# save output to file
write_tsv(output, file = here("resources/ENCODE_lookup.combined.unique.scratch.abc.tsv"))
