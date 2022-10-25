## Get 1bp TSS annotations based on extended TSS window annotations such as those used by ABC

suppressPackageStartupMessages(library(rtracklayer))

# load tss annotations
tss <- import(snakemake@input[[1]])

# resize annotations to get centers of windows
tss_centers <- resize(tss, width = 1, fix = "center")

# export to bed file
export(tss_centers, con = snakemake@output[[1]], format = "bed")
