## Get 1bp TSS annotations based on extended TSS window annotations such as those used by ABC

suppressPackageStartupMessages({
  library(readr)
  library(rtracklayer)
})

# load tss annotations
tss <- read_tsv(snakemake@input[[1]], show_col_types = FALSE, col_select = 1:6)

# resize annotations to get centers of windows
tss <- makeGRangesFromDataFrame(tss, seqnames.field = "#chr", keep.extra.columns = TRUE)
tss_centers <- resize(tss, width = 1, fix = "center")

# export to bed file
export(tss_centers, con = snakemake@output[[1]], format = "bed")
