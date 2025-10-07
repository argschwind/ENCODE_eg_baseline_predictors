## Rules to copy baseline predictors from scratch to oak and organize files

# directory containing created baseline predictors
outdir = config["outdir"]

# directory on oak where predictors are stored
destdir = "/oak/stanford/groups/engreitz/Projects/Benchmarking/Revisions/Predictors/BaselinePredictors"

# copy all required baseline predictors to oak
rule all_oak:
  input:
    expand(destdir + "/{peaks}/{sample}/{pred}.tsv.gz", peaks = ["DHS", "ABC"],
           sample = config["dnase_pred_samples"] + config["h3ka7ac_pred_samples"],
           pred = ["dist_to_tss", "dist_to_gene", "within_100kb_of_tss", "within_100kb_of_gene",
                   "within_100kb_of_expressed_tss", "within_100kb_of_expressed_gene", "nearest_tss",
                   "nearest_gene", "nearest_expressed_tss", "nearest_expressed_gene"]),
    expand(destdir + "/{peaks}/{sample}/DHS_reads_by_dist_to_{annot}{ext}",
           sample = config["dnase_pred_samples"], peaks = ["DHS", "ABC"], annot = ["tss", "gene"],
           ext = [".tsv.gz", "_norm.tsv.gz"]),
    expand(destdir + "/{peaks}/{sample}/H3K27ac_reads_by_dist_to_{annot}{ext}",
           sample = config["h3ka7ac_pred_samples"], peaks = ["DHS", "ABC"], annot = ["tss", "gene"],
           ext = [".tsv.gz", "_norm.tsv.gz"])        

rule copy_dhs_baseline_predictor:
  input: outdir + "/{sample}/dhs/{pred}.tsv.gz"
  output: destdir + "/DHS/{sample}/{pred}.tsv.gz"
  shell:
    "cp {input} {output}"

rule copy_abc_baseline_predictor:
  input: outdir + "/{sample}/abc/{pred}.tsv.gz"
  output: destdir + "/ABC/{sample}/{pred}.tsv.gz"
  shell:
    "cp {input} {output}"
