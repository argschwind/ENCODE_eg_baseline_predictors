
# output directory
outdir = config["outdir"]

# set rule orders
ruleorder: sort_bam > download_bam

# global wildcard constraints
wildcard_constraints:
  annot = "tss|gene"

# input function to get peak files for a given sample  
def get_peak_files(wildcards):
  sample = wildcards.sample
  peaks = wildcards.peaks
  if peaks == "dhs":
    file = outdir + "/" + sample + "/" + config["samples"][sample]["DHS_peaks"] + ".bed.gz"
  elif peaks == "abc":
    file = config["samples"][sample]["ABC_peaks"]
  else:
    print("Invalid 'type' wildcard")
    sys.exit()
  return file

# input function to get bam files for a given sample
def get_bam_files(wildcards):
  file_id = config["samples"][wildcards.sample][wildcards.assay + "_bam"]
  bam = outdir + "/" + wildcards.sample + "/" + file_id + ".sorted.bam"
  bai = outdir + "/" + wildcards.sample + "/" + file_id + ".sorted.bam.bai"
  return({"bam": bam, "bai": bai})


## Process annotation files ------------------------------------------------------------------------

# get 1bp TSS coordinates from TSS annotations
rule get_tss_coords:
  input: config["annot"]["tss"]
  output: temp("resources/genome_annotations/tss_annot.bed")
  conda: "../envs/baseline_predictors.yml"
  script:
    "../scripts/get_tss_coords.R"
    
# get primary columns from gene annotations
rule get_gene_annots:
  input: config["annot"]["gene"]
  output: temp("resources/genome_annotations/gene_annot.bed")
  conda: "../envs/baseline_predictors.yml"
  shell:
    """
    awk 'BEGIN {{OFS = "\\t"}} (NR>1) {{print $1, $2, $3, $4, $5, $6}}' {input} > {output}
    """

# filter genome annotations for regular chromosomes only and sort according to position
rule sort_annot:
  input:
    annot = "resources/genome_annotations/{annot}_annot.bed",
    reg_chrs = config["annot"]["regular_chrs"],
    chrs = config["annot"]["chrs"],
  output: temp("resources/genome_annotations/{annot}_annot.sorted.bed")
  conda: "../envs/baseline_predictors.yml"
  shell:
    """
    bedtools intersect -a {input.annot} -b {input.reg_chrs} -wa | \
    bedtools sort -i stdin -g {input.chrs} > {output}
    """
    
# get expressed genes for a given cell type from ABC output
rule get_expressed_genes:
  input:
    abc = lambda wildcards: config["samples"][wildcards.sample]["ABC_preds"],
    annot = "resources/genome_annotations/{annot}_annot.sorted.bed"
  output: temp(outdir + "/{sample}/expressed_{annot}.sorted.bed")
  conda: "../envs/baseline_predictors.yml"
  resources:
    mem = "8G"
  script:
    "../scripts/get_expressed_genes.R"
    
    
## Download ENCODE data ----------------------------------------------------------------------------

# download a DNase-seq narrowPeak DHS file from the ENCODE portal
rule download_encode_dhs:
  output: temp(outdir + "/{sample}/{file}.bed.gz")
  params:
    url = "https://www.encodeproject.org/files/{file}/@@download/{file}.bed.gz"
  conda: "../envs/baseline_predictors.yml"
  shell:
    "wget {params.url} -O {output}"
    
# download bam files
rule download_bam:
  output: temp(outdir + "/{sample}/{file}.bam")
  params:
    url = "https://www.encodeproject.org/files/{file}/@@download/{file}.bam"
  conda: "../envs/baseline_predictors.yml"
  shell:
    "wget {params.url} -O {output}"
    
# sort bam file and create index file
rule sort_bam:
  input: outdir + "/{sample}/{file}.bam"
  output:
    bam = temp(outdir + "/{sample}/{file}.sorted.bam"),
    bai = temp(outdir + "/{sample}/{file}.sorted.bam.bai")
  conda: "../envs/baseline_predictors.yml"
  resources:
    mem = "32G"
  shell:
    "samtools sort -o {output.bam} {input}; samtools index {output.bam}"
  

## Create CRE - TSS/gene pairs ---------------------------------------------------------------------

# sort peak files, merge overlapping peaks and create list of unique candidate elements
rule create_cres:
  input:
    peaks = get_peak_files,
    chrs = config["annot"]["chrs"]
  output: temp(outdir + "/{sample}/{peaks}_elements.bed")
  conda: "../envs/baseline_predictors.yml"
  resources:
    mem = "8G"
  shell:
    """
    zcat -f {input.peaks} | cut -f 1,2,3 | \
    bedtools sort -i stdin -g {input.chrs} | bedtools merge -i stdin | \
    awk 'BEGIN{{OFS="\\t"}} {{print $1, $2, $3, $1":"$2"-"$3, "0", "."}}' > {output}
    """

# build CRE-gene pairs based on TSS or gene annotations
rule build_cre_gene_pairs:
  input:
    cres = outdir + "/{sample}/{peaks}_elements.bed",
    annot = "resources/genome_annotations/{annot}_annot.sorted.bed"
  output: temp(outdir + "/{sample}/{peaks}_{annot}_pairs.txt.gz")
  params:
    max_dist = 1e6
  conda: "../envs/baseline_predictors.yml"
  resources:
    mem = "8G"
  shell:
    "bedtools window -a {input.cres} -b {input.annot} -w {params.max_dist} | "
    "python workflow/scripts/annotate_pairs_distance.py | "
    "gzip > {output}"
 
 
## Count chromatin assay reads in CREs -------------------------------------------------------------

# count reads overlapping cres
rule count_reads:
  input:
    unpack(get_bam_files),
    cres = outdir + "/{sample}/{peaks}_elements.bed",
    chrs = config["annot"]["chrs"]
  output: temp(outdir + "/{sample}/{peaks}_elements.{assay}_read_counts.bed")
  conda: "../envs/baseline_predictors.yml"
  resources:
    mem = "4G"
  shell:
    "bedtools coverage -a {input.cres} -b {input.bam} -g {input.chrs} -counts -sorted > {output}"

# normamlize read counts to 1 million total reads mapping to DHS
rule normalize_counts:
  input: outdir + "/{sample}/{peaks}_elements.{assay}_read_counts.bed"
  output: temp(outdir + "/{sample}/{peaks}_elements.{assay}_read_counts.normalized.bed")
  conda: "../envs/baseline_predictors.yml"
  resources:
    mem = "4G"
  script:
    "../scripts/normalize_read_counts.R"


## Create baseline predictors ----------------------------------------------------------------------

# create 'distance to TSS' or 'distance to gene body' predictor
rule distance_pred:
  input:
    pairs = outdir + "/{sample}/{peaks}_{annot}_pairs.txt.gz",
  output: outdir + "/{sample}/{peaks}/dist_to_{annot}.tsv.gz"
  params:
    cell_type = lambda wildcards: config["samples"][wildcards.sample]["cell_type"]
  conda: "../envs/baseline_predictors.yml"
  resources:
    mem = "4G"
  shell:
    """
    python workflow/scripts/compute_distance_predictor.py -i {input.pairs} -c "{params.cell_type}" \
    -t {wildcards.annot} | gzip > {output}
    """

# compute boolean 'within distance to TSS' or 'within distance to gene body' predictor
rule within_distance_pred:
  input:
    pairs = outdir + "/{sample}/{peaks}_{annot}_pairs.txt.gz",
    annot = "resources/genome_annotations/{annot}_annot.sorted.bed"
  output: outdir + "/{sample}/{peaks}/within_{dist}kb_of_{annot}.tsv.gz"
  params:
    cell_type = lambda wildcards: config["samples"][wildcards.sample]["cell_type"]  
  conda: "../envs/baseline_predictors.yml"
  resources:
    mem = "4G"
  shell:
    """
    python workflow/scripts/compute_within_distance_predictor.py -i {input.pairs} \
    -c "{params.cell_type}" -d {wildcards.dist} -a {input.annot} -t {wildcards.annot} | \
    gzip > {output}
    """
    
# compute boolean 'within distance to expressed TSS' or 'within distance to expressed gene' predictor
rule within_distance_expressed_pred:
  input:
    pairs = outdir + "/{sample}/{peaks}_{annot}_pairs.txt.gz",
    annot = outdir + "/{sample}/expressed_{annot}.sorted.bed"
  output: outdir + "/{sample}/{peaks}/within_{dist}kb_of_expressed_{annot}.tsv.gz"
  params:
    cell_type = lambda wildcards: config["samples"][wildcards.sample]["cell_type"]  
  conda: "../envs/baseline_predictors.yml"
  resources:
    mem = "4G"
  shell:
    """
    python workflow/scripts/compute_within_distance_predictor.py -i {input.pairs} \
    -c "{params.cell_type}" -d {wildcards.dist} -a {input.annot} -t {wildcards.annot} | \
    gzip > {output}
    """

# compute boolean 'is nearest TSS' or 'is nearest gene' predictor
rule nearest_tss_or_gene_pred:
  input:
    pairs = outdir + "/{sample}/{peaks}_{annot}_pairs.txt.gz",
    annot = "resources/genome_annotations/{annot}_annot.sorted.bed",
    chrs = config["annot"]["chrs"]
  output: outdir + "/{sample}/{peaks}/nearest_{annot}.tsv.gz"
  params:
    cell_type = lambda wildcards: config["samples"][wildcards.sample]["cell_type"]   
  conda: "../envs/baseline_predictors.yml"
  resources:
    mem = "4G"
  shell:
    """
    bedtools closest -a {input.pairs} -b {input.annot} -g {input.chrs} | cut -f1-15,19 | \
    python workflow/scripts/compute_nearest_tss_or_gene.py -c "{params.cell_type}" | gzip > {output}
    """

# compute boolean 'is nearest expressed TSS' or 'is nearest expressed gene' predictor
rule nearest_expressed_tss_or_gene_pred:
  input:
    pairs = outdir + "/{sample}/{peaks}_{annot}_pairs.txt.gz",
    annot = outdir + "/{sample}/expressed_{annot}.sorted.bed",
    chrs = config["annot"]["chrs"]
  output: outdir + "/{sample}/{peaks}/nearest_expressed_{annot}.tsv.gz"
  params:
    cell_type = lambda wildcards: config["samples"][wildcards.sample]["cell_type"]     
  conda: "../envs/baseline_predictors.yml"
  resources:
    mem = "4G"
  shell:
    """
    bedtools closest -a {input.pairs} -b {input.annot} -g {input.chrs} | cut -f1-15,19 | \
    python workflow/scripts/compute_nearest_tss_or_gene.py -c "{params.cell_type}" | gzip > {output}
    """

# compute 'reads by distance to TSS' or 'reads by distance to gene' predictor
rule reads_by_dist:
  input:
    reads = outdir + "/{sample}/{peaks}_elements.{assay}_read_counts.normalized.bed",
    annot = "resources/genome_annotations/{annot}_annot.sorted.bed"
  output: outdir + "/{sample}/{peaks}/{assay}_reads_by_dist_to_{annot}.tsv.gz"
  params:
    max_dist = 1e6,
    cell_type = lambda wildcards: config["samples"][wildcards.sample]["cell_type"]
  conda: "../envs/baseline_predictors.yml"
  resources:
    mem = "4G"
  shell:
    """
    bedtools window -a {input.reads} -b {input.annot} -w {params.max_dist} | \
    python workflow/scripts/annotate_pairs_distance.py -c 7 | \
    python workflow/scripts/compute_reads_by_distance.py -c "{params.cell_type}" -t {wildcards.annot} | \
    gzip > {output}
    """    
 
# normalize reads by distance predictor per gene
rule normalize_reads_by_dist:
  input: outdir + "/{sample}/{peaks}/{assay}_reads_by_dist_to_{annot}.tsv.gz"
  output: outdir + "/{sample}/{peaks}/{assay}_reads_by_dist_to_{annot}_norm.tsv.gz"
  conda: "../envs/baseline_predictors.yml"
  resources:
    mem = "4G"
  script:
    "../scripts/normalize_reads_by_distance.R"
