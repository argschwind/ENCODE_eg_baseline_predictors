## Manually add DNase-seq and H3K27ac bam files for K256 and GM12878 to metadata table

library(tidyverse)
library(here)

# load sample metadata table
e2g_input_files <- read_tsv(here("config/e2g_baseline_preds_input_files.tsv"),
                            show_col_types = FALSE)

# load table containing manually currated DNase-seq and bam files
bam_files <- read_tsv(here("config/bam_files_encode_portal.tsv"), show_col_types = FALSE)

# add to meta data table
e2g_input_files <- e2g_input_files %>% 
  left_join(bam_files, by = "sample") %>% 
  replace_na(replace = list(DHS_bam = "Null", H3K27ac_bam = "Null"))

# save to output file
write_tsv(e2g_input_files, file = here("config/e2g_baseline_preds_input_files_bam.tsv"))
