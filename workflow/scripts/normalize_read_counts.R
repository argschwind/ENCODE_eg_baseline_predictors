## Normalize read counts within DHS to 1 million total reads

# required packages
library(data.table)

# load read counts
counts <- fread(snakemake@input[[1]], header = FALSE)

# normalize read counts to 1 million total reads
total_reads <- sum(counts$V7)
counts$V7 <- counts$V7 / sum(counts$V7) * 1e6

# write to output file
fwrite(counts, file = snakemake@output[[1]], sep = "\t", col.names = FALSE,  quote = FALSE,
       na = ".")
