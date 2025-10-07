# ENCODE4 distal regulation E-G baseline predictors
Snakemake workflow to generate simple baseline predictors for 1,458 DNase-seq experiments used in
the ENCODE distal regulation paper.

## Baseline predictors
This workflow computes simple baseline predictors for element-gene universes based on either ENCODE
DHS peaks of ABC candidate elements. For both universes, element-gene (E-G) pairs are created by
pairing gene TSSs with all candidate elements within 1Mb.

The workflow computes following baseline predictors for all 1,458 DNase-seq experiments (both DHS
and ABC element universes):
- Distance to TSS
- Nearest expressed gene

For a subset of 24 selected experiments, a more comprehensive list of baseline predictors are
computed (both DHS and ABC element universes):
- Distance to TSS
- Distance to gene
- Nearest TSS
- Nearest gene
- Nearest expressed TSS
- Nearest expressed gene
- Within 100kb of TSS
- Within 100kb of gene
- Within 100kb of expressed TSS
- Within 100kb of expressed gene
- DNase-seq reads * 1/distance to TSS
- DNase-seq reads * 1/distance to TSS, normalized to 1 per gene
- DNase-seq reads * 1/distance to gene
- DNase-seq reads * 1/distance to gene, normalized to 1 per gene

For 20 of the 24 selected experiments, additional baseline predictors based on H3K27ac ChIP-seq are
computed (both DHS and ABC element universes):
- H3K27ac ChIP-seq reads * 1/distance to TSS
- H3K27ac ChIP-seq reads * 1/distance to TSS, normalized to 1 per gene
- H3K27ac ChIP-seq reads * 1/distance to gene
- H3K27ac ChIP-seq reads * 1/distance to gene, normalized to 1 per gene

## Required inputs
The workflow handles downloading of ENCODE DHS peak and DNase-seq and H3L27ac ChIP-seq bam files. To
generate baseline predictors for ABC candidate elements, ABC output files are required.

All input files including ABC predictions are listed in the
`config/e2g_baseline_preds_input_files_bam.tsv` config table. To run the workflow, we suggest to
download ABC candidate elements and ABC prediction files from the ENCODE portal and update the table
with paths to the downloaded files.

The `config/2_create_sample_table.R` helper script is used to create the config table and can be
modified to update file paths of ABC predictions (lines 14-22).

## Executing the workflow
To run the snakemake workflow, make sure that both snakemake and conda are installed. Snakemake can
uses conda to create environments containing all other required dependencies. Due to the large
number of steps, the workflow should be ran on an HPC cluster that allows parallelization of
snakemake jobs (see snakemake documentation for instructions).

The workflow will download and generate a large number of temporary files. It's recommended to use a
dedicated directory for all download and temporary files, for example on the scratch space of the
used HPC. This can be configured by modifying the `outdir` entry in the `config/config.yml` file.

To generate all baseline predictors, simply run the full workflow by executing:

```sh
# execute full workflow to generate all output files
snakemake --use-conda --profile <cluster_profile> -j100 -n
```

Edit the `all` rule in the main Snakefile (`workfow/Snakefile`) to generate custom subsets of
possible baseline predictors.

## Pre-computed baseline predictors
Files for all baseline predictors computed by this workflow can be found on synapse.org under: 
https://www.synapse.org/Synapse:syn58896208
