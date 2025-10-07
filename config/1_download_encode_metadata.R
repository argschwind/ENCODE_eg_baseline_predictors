## Query the ENCODE portal for default DNase-seq peak files and save to metadata table

library(tidyverse)
library(here)

# query ENCODE portal for all default DNase-seq peak files
url <- "https://www.encodeproject.org/report.tsv?type=File&file_format=bed&preferred_default=true&assay_title=DNase-seq&file_format_type=narrowPeak&status=released&field=analyses.title&field=accession&field=dataset&field=assembly&field=biosample_ontology.term_name&field=href"
dhs <- read_tsv(url, skip = 1, show_col_types = FALSE)

# filter for human ENCODE4 files
dhs <- filter(dhs, `Genome assembly` == "GRCh38", grepl("^ENCODE4", analyses.title))

# write to output table
write_tsv(dhs, here("config/encode_default_dhs_files.tsv"))
