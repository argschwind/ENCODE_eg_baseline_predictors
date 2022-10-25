
# global wildcard constraints
wildcard_constraints:
  type = "tss|gene"

## Process input files -----------------------------------------------------------------------------

# get 1bp TSS coordinates from TSS annotations
rule get_tss_coords:
  input: config["annot"]["tss_window"]
  output: temp(config["annot"]["tss"])
  conda: "../envs/baseline_predictors.yml"
  script:
    "../scripts/get_tss_coords.R"
    
# filter genome annotations for regular chromosomes only
rule filt_annot:
  input:
    annot = lambda wildcards: config["annot"][wildcards.type],
    reg_chrs = config["annot"]["regular_chrs"]
  output: temp("resources/genome_annotations/{type}_annot.bed")
  conda: "../envs/baseline_predictors.yml"
  shell:
    "bedtools intersect -a {input.annot} -b {input.reg_chrs} -wa > {output}"

# sort genome annotation files
rule sort_annot:
  input:
    annot = "resources/genome_annotations/{type}_annot.bed",
    chrs = config["annot"]["chrs"]
  output: "resources/genome_annotations/{type}_annot.sorted.bed"
  conda: "../envs/baseline_predictors.yml"
  shell:
    "bedtools sort -i {input.annot} -g {input.chrs} > {output}"

# sort DHS narrowPeaks files, merge overlapping peaks and extract unique DHS
rule extract_unique_dhs:
  input:
    peaks = lambda wildcards: list(str.split(config["samples"][wildcards.sample]['DHS_narrowPeak'], ",")),
    chrs = config["annot"]["chrs"]
  output: temp("results/DHS/{sample}/unique_dhs_narrowPeak.bed")
  conda: "../envs/baseline_predictors.yml"
  shell:
    "zcat {input.peaks} | bedtools sort -i stdin -g {input.chrs} | bedtools merge -i stdin | "
    """awk 'BEGIN{{OFS="\\t"}} {{print $1, $2, $3, $1":"$2"-"$3, "0", "."}}' > {output}"""
    
# get expressed genes for a given cell type from ABC output
rule get_expressed_genes:
  input:
    abc = lambda wildcards: config["samples"][wildcards.sample]["ABC_file"],
    annot = "resources/genome_annotations/{type}_annot.sorted.bed"
  output: temp("results/{sample}/expressed_{type}.bed")
  conda: "../envs/baseline_predictors.yml"
  script:
    "../scripts/get_expressed_genes.R"

## Create CRE - TSS/gene pairs ---------------------------------------------------------------------

# build CRE-gene pairs based on TSS or gene annotations
rule build_cre_gene_pairs:
  input:
    cres = "results/DHS/{sample}/unique_dhs_narrowPeak.bed",
    annot = "resources/genome_annotations/{type}_annot.sorted.bed"
  output: temp("results/DHS/{sample}/dhs_{type}_pairs.txt.gz")
  params:
    max_dist = 1e6
  conda: "../envs/baseline_predictors.yml"
  shell:
    "bedtools window -a {input.cres} -b {input.annot} -w {params.max_dist} | "
    "python workflow/scripts/annotate_pairs_distance.py | "
    "gzip > {output}"

## Count chromatin assay reads in CREs -------------------------------------------------------------

# count reads overlapping cres
rule count_reads:
  input:
    cres = "results/DHS/{sample}/unique_dhs_narrowPeak.bed",
    reads = lambda wildcards: list(str.split(config["samples"][wildcards.sample][wildcards.assay + "_tagAlign"], ",")),
    chrs = config["annot"]["chrs"]
  output: temp("results/DHS/{sample}/unique_dhs_narrowPeak.{assay}_read_counts.bed")
  conda: "../envs/baseline_predictors.yml"
  shell:
    "bedtools coverage -a {input.cres} -b {input.reads} -g {input.chrs} -counts -sorted > {output}"

# normamlize read counts to 1 million total reads mapping to DHS
rule normalize_counts:
  input: "results/DHS/{sample}/unique_dhs_narrowPeak.{assay}_read_counts.bed"
  output: temp("results/DHS/{sample}/unique_dhs_narrowPeak.{assay}_read_counts.normalized.bed")
  conda: "../envs/baseline_predictors.yml"
  script:
    "../scripts/normalize_read_counts.R"

## Create baseline predictors ----------------------------------------------------------------------

# create 'distance to TSS' or 'distance to gene body' predictor
rule distance_pred:
  input:
    pairs = "results/DHS/{sample}/dhs_{type}_pairs.txt.gz",
  output: "results/DHS/{sample}/dist_to_{type}.tsv.gz"
  conda: "../envs/baseline_predictors.yml"
  shell:
    "python workflow/scripts/compute_distance_predictor.py -i {input.pairs} -c {wildcards.sample} "
    "-t {wildcards.type} | gzip > {output}"

# compute boolean 'within distance to TSS' or 'within distance to gene body' predictor
rule within_distance_pred:
  input:
    pairs = "results/DHS/{sample}/dhs_{type}_pairs.txt.gz",
  output: "results/DHS/{sample}/within_{dist}kb_of_{type}.tsv.gz"
  conda: "../envs/baseline_predictors.yml"
  shell:
    "python workflow/scripts/compute_within_distance_predictor.py -i {input.pairs} "
    "-c {wildcards.sample} -d {wildcards.dist} -t {wildcards.type} | gzip > {output}"

# compute boolean 'is nearest TSS' or 'is nearest gene' predictor
rule nearest_tss_or_gene_pred:
  input:
    pairs = "results/DHS/{sample}/dhs_{type}_pairs.txt.gz",
    annot = "resources/genome_annotations/{type}_annot.sorted.bed",
    chrs = config["annot"]["chrs"]
  output: "results/DHS/{sample}/nearest_{type}.tsv.gz"
  conda: "../envs/baseline_predictors.yml"
  shell:
    "bedtools closest -a {input.pairs} -b {input.annot} -g {input.chrs} | cut -f1-15,19 | "
    "python workflow/scripts/compute_nearest_tss_or_gene.py -c {wildcards.sample} | gzip > {output}"
    
# compute boolean 'is nearest expressed TSS' or 'is nearest expressed gene' predictor
rule nearest_expressed_tss_or_gene_pred:
  input:
    pairs = "results/DHS/{sample}/dhs_{type}_pairs.txt.gz",
    annot = "results/{sample}/expressed_{type}.bed",
    chrs = config["annot"]["chrs"]
  output: "results/DHS/{sample}/nearest_expressed_{type}.tsv.gz"
  conda: "../envs/baseline_predictors.yml"
  shell:
    "bedtools closest -a {input.pairs} -b {input.annot} -g {input.chrs} | cut -f1-15,19 | "
    "python workflow/scripts/compute_nearest_tss_or_gene.py -c {wildcards.sample} | gzip > {output}"

# compute 'reads by distance to TSS' or 'reads by distance to gene' predictor
rule reads_by_dist:
  input:
    reads = "results/DHS/{sample}/unique_dhs_narrowPeak.{assay}_read_counts.normalized.bed",
    annot = "resources/genome_annotations/{type}_annot.sorted.bed"
  output: "results/DHS/{sample}/{assay}_reads_by_dist_to_{type}.tsv.gz"
  params:
    max_dist = 1e6,
    cre_cols = 7
  conda: "../envs/baseline_predictors.yml"
  shell:
    "bedtools window -a {input.reads} -b {input.annot} -w {params.max_dist} | "
    "python workflow/scripts/annotate_pairs_distance.py -c {params.cre_cols} | "
    "python workflow/scripts/compute_reads_by_distance.py -c {wildcards.sample} -t {wildcards.type} | "
    "gzip > {output}"

# normalize reads by distance predictor per gene
rule normalize_reads_by_dist:
  input: "results/DHS/{sample}/{assay}_reads_by_dist_to_{type}.tsv.gz"
  output: "results/DHS/{sample}/{assay}_reads_by_dist_to_{type}_norm.tsv.gz"
  conda: "../envs/baseline_predictors.yml"
  script:
    "../scripts/normalize_reads_by_distance.R"
