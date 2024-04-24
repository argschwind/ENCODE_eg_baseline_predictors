## Add ABC prediction files to lookup table

# required packages
library(tidyverse)
library(here)

# lookup table for baseline predictors
lookup <- read_tsv(here("resources/ENCODE_lookup.combined.unique.scratch.tsv"))

# DNase + H3K27ac AB samples
abc <- read_tsv("resources/H3K27ac_DNase_qc_samples.tsv",
                col_names = c("ABC_sample", "DHS_tagAlign", "H3K27ac_tagAlign"))

# base directory for all abc predictions
abc_dir <- "/oak/stanford/groups/akundaje/kmualim/02212022_hg38Preds_QNORM-avg_hic_track2"

## process lookup table ----------------------------------------------------------------------------

# get DHS file ids for each biosample
lookup <- lookup %>%
  mutate(DHS_ids = strsplit(DHS_tagAlign, split = ",", fixed = TRUE)) %>% 
  rowwise() %>% 
  mutate(DHS_ids = list(basename(DHS_ids))) %>% 
  mutate(DHS_ids = list(gsub("[.se]*.tagAlign.gz", "", DHS_ids))) %>% 
  mutate(DHS_ids_string = list(paste(sort(DHS_ids), collapse = "_")))

# get H3K27ac file ids for each biosample
lookup <- lookup %>%
  mutate(H3K27ac_ids = strsplit(H3K27ac_tagAlign, split = ",", fixed = TRUE)) %>% 
  rowwise() %>% 
  mutate(H3K27ac_ids = list(basename(H3K27ac_ids))) %>% 
  mutate(H3K27ac_ids = list(gsub("[.se]*.tagAlign.gz", "", H3K27ac_ids))) %>% 
  mutate(H3K27ac_ids_string = list(paste(sort(H3K27ac_ids), collapse = "_")))
 
# combine DHS and H3K27ac files into new files id
lookup <- unite(lookup, col = "files_string", DHS_ids_string, H3K27ac_ids_string)

## process abc table -------------------------------------------------------------------------------

# get DHS file ids for each biosample
abc <- abc %>%
  mutate(DHS_ids = strsplit(DHS_tagAlign, split = ",", fixed = TRUE)) %>% 
  rowwise() %>% 
  mutate(DHS_ids = list(basename(DHS_ids))) %>% 
  mutate(DHS_ids = list(gsub("[.se]*.bam", "", DHS_ids))) %>% 
  mutate(DHS_ids_string = list(paste(sort(DHS_ids), collapse = "_")))

# get H3K27ac file ids for each biosample
abc <- abc %>%
  mutate(H3K27ac_ids = strsplit(H3K27ac_tagAlign, split = ",", fixed = TRUE)) %>% 
  rowwise() %>% 
  mutate(H3K27ac_ids = list(basename(H3K27ac_ids))) %>% 
  mutate(H3K27ac_ids = list(gsub("[.se]*.bam", "", H3K27ac_ids))) %>% 
  mutate(H3K27ac_ids_string = list(paste(sort(H3K27ac_ids), collapse = "_")))

# combine DHS and H3K27ac files into new files id
abc <- unite(abc, col = "files_string", DHS_ids_string, H3K27ac_ids_string)

## add abc sample ids to lookup table --------------------------------------------------------------

# add abc samples whenever the exact same files are used as input
output <- abc %>% 
  select(ABC_sample, files_string) %>% 
  left_join(lookup, ., by = "files_string")

# add path to abc predictions file for matched samples
output <- output %>% 
  mutate(
    ABC_file = if_else(!is.na(ABC_sample),
      true = file.path(abc_dir, ABC_sample, "Predictions/EnhancerPredictionsAllPutative.txt.gz"),
      false = NA_character_)
    )

# remove unused columns and reformat
output <- output %>% 
  select(-c(DHS_ids, H3K27ac_ids, files_string)) %>% 
  replace_na(replace = list(ABC_sample = "none", ABC_file = "none"))

# save output to file
write_tsv(output, file = here("resources/ENCODE_lookup.combined.unique.scratch.abc.tsv"))
