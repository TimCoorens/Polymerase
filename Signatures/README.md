# Workflow

For this paper, the signature extraction and assignment consisted roughly of three steps: extracting mutational signatures using the hdp packages, deconvoluting the obtained signatures into reference PCAWG signatures (where appropriate) and fitting the resulting final set of signatures to the observed counts of SNVs per branch or per bulk sample. 

### 1. Running HDP

The starting point of these analyses is a simple matrix of the 96 trinucleotide counts per branch. A function to convert a list of mutations (containing information on chromosome, position, reference genotype and alternative genotype) is contained in the R script labeled 'hdp_workflow.R' (mutlist_to_96_contexts()). If a branch contained more than 3000 SNVs, we've randomly sampled it down to that number to prevent difficulties in running HDP. In addition, we've imposed a minimum of 100 SNVs per branch to be included in the signature extraction step, to avoid overfitting and noise influencing the resulting signatures.

We then set off 20 independent chains of hdp (part 2 in hdp_workflow.R). The only hierarchy imposed is a per-patient parent node for all branches.

The last part of running HDP is combining the results of the 20 independent chains, plotting diagnostic figures and generating matrices for the mutational signatures and the mean exposure per signature per branch.

### 2. Deconvoluting HDP signatures into reference signatures

### 3. Re-fitting signatures to counts.



