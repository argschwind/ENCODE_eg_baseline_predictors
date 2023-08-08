# ENCODE4 distal regulation E-G baseline predictors
Build distal regulation baseline E-G predictors using ENCODE DHSs (89 biosamples) or ABC candidates (58 biosamples) as candidate element universes.

**Following baseline predictors are created for each biosample:**
* Distance to TSS
* Distance to gene
* Within 100kb of TSS (binary predictor)
* Within 100kb of gene (binary predictor)
* Nearest TSS (binary predictor)
* Nearest gene (binary predictor)
* DNase-seq reads * 1/distance to TSS
* DNase-seq reads * 1/distance to TSS, normalized to 1 per gene
* H3K27ac ChIP-seq reads * 1/distance to TSS
* H3K27ac ChIP-seq reads * 1/distance to TSS, normalized to 1 per gene
* Nearest expressed TSS (binary predictor, only for biosamples with ABC predictions)
* Nearest expressed gene (binary predictor, only for biosamples with ABC predictions)

Files containing all created baseline predictors can be found on synapse.org under: https://www.synapse.org/#!Synapse:syn52234396
