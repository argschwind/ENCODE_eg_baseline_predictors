# filter a gene/TSS bed file for expressed genes from ABC output for a given cell type

# required packages
suppressPackageStartupMessages({
  library(data.table)
})

# load gene/TSS annotations
annot <- fread(snakemake@input$annot, header = FALSE)

# load abc results
abc <- fread(snakemake@input$abc)

# extract expressed genes from abc results
expr_genes <- unique(abc[abc$TargetGeneIsExpressed == TRUE, TargetGene])

# filter annotations for these genes
annot_filt <- annot[annot$V4 %in% expr_genes, ]

# write to output file
fwrite(annot_filt, file = snakemake@output[[1]], sep = "\t", col.names = FALSE)
