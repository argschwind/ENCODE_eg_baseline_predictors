## Create sample table for baseline predictors containing paths or ENCODE accession ids for all
## input files  

library(tidyverse)
library(here)

# Load input metadata ------------------------------------------------------------------------------

# download ENCODE-rE2G metadata table and load into R
url <- "https://docs.google.com/spreadsheets/d/1aJ44UU5S09ZVPBmwkaucm8Zxc7_30uWzkKt4p_eGRfI/export?format=tsv"
download.file(url, destfile = here("config/encode_re2g_metadata.tsv"))
e2g_metadata <- read_tsv(here("config/encode_re2g_metadata.tsv"), show_col_types = FALSE)

# get all ABC candidate element files for all ENCODE-rE2G samples
e2g_dir <- "/oak/stanford/groups/engreitz/Users/atan5133/encode_dataset_processing/results/dhs_only"
e2g_samples <- list.files(e2g_dir, pattern = "_ENCSR", full.names = TRUE)
abc_elements <- file.path(e2g_samples, "Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed")
names(abc_elements) <- basename(e2g_samples)

# get all ABC prediction files to define expressed genes
abc_predictions <- file.path(e2g_samples, "Predictions/EnhancerPredictionsAllPutative.tsv.gz")
names(abc_predictions) <- basename(e2g_samples)

# create a table with all ABC candidate elements and predictions per sample
abc_elements <- enframe(abc_elements, name = "local_sample", value = "ABC_peaks")
abc_predictions <- enframe(abc_predictions, name = "local_sample", value = "ABC_preds")
abc_files <- inner_join(abc_elements, abc_predictions, by = "local_sample")

# load default DHS peaks for each ENCODE DNase-seq experiment
dhs <- read_tsv(here("config/encode_default_dhs_files.tsv"), show_col_types = FALSE)

# Create DNase-only baseline predictor input files table -------------------------------------------

# add local sample name to e2g_metadata table
e2g_metadata <- e2g_metadata %>% 
  mutate(local_sample = gsub("'", "", Sample)) %>% 
  mutate(local_sample = gsub("/", "_", local_sample))

# add ABC files to ENCODE-rE2G metadata table based on local sample id
e2g_metadata <- left_join(e2g_metadata, abc_files, by = "local_sample")

# add default DHS peaks file to ENCODE-rE2G metadata table based on DNase-seq experiment id
e2g_metadata <- dhs %>% 
  mutate(Dataset = basename(Dataset),
         `Download URL` = paste0("https://www.encodeproject.org", `Download URL`)) %>% 
  select(DHS_peaks = Accession, `DNase Experiment accession` = Dataset) %>% 
  left_join(e2g_metadata, ., by = "DNase Experiment accession")

# select relevant columns for baseline predictors pipeline
e2g_input_files <- e2g_metadata %>% 
  select(sample = local_sample, ABC_peaks, ABC_preds, DHS_peaks, cell_type = `Biosample term name`)

# save to .tsv file to use in pipeline
write_tsv(e2g_input_files, file = here("config/e2g_baseline_preds_input_files.tsv"))
