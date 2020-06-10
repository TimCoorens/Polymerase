# Workflow

For this paper, the signature extraction and assignment consisted roughly of three steps: extracting mutational signatures using the hdp package (https://github.com/nicolaroberts/hdp), deconvoluting the obtained signatures into reference PCAWG signatures (where appropriate) and fitting the resulting final set of signatures to the observed counts of SNVs per branch or per bulk sample. 

### 1. Running HDP

The starting point of these analyses is a simple matrix of the 96 trinucleotide counts per branch. A function to convert a list of mutations (containing information on chromosome, position, reference genotype and alternative genotype) is contained in the R script labeled 'hdp_workflow.R' (mutlist_to_96_contexts()). If a branch contained more than 3000 SNVs, we've randomly sampled it down to that number to prevent difficulties in running HDP. In addition, we've imposed a minimum of 100 SNVs per branch to be included in the signature extraction step, to avoid overfitting and noise influencing the resulting signatures.

We then set off 20 independent chains of hdp (part 2 in hdp_workflow.R). The only hierarchy imposed is a per-patient parent node for all branches.

The last part of running HDP is combining the results of the 20 independent chains, plotting diagnostic figures and generating matrices for the mutational signatures and the mean exposure per signature per branch.

### 2. Deconvoluting HDP signatures into reference signatures

The next step is to compare the obtained hdp signatures to a set of reference signatures. Here, we used the mutational signatures from PCAWG (Alexandrov et al, 2020) and a set of novel signatures found in normal colon (Lee-Six et al, 2019). Both sets of references are uploaded here (PCAWG_sigProfiler_SBS_signatures.csv and Signatures_trinucleotide_composition_SBS_ABCD.txt)

For our hdp signatures that are combinations of known signatures, we deconvolve using an EM-algorithm, largely based on the approach previously described in Lee-Six et al, 2019. The code for this step can be found in deconvolute_hdp_sigs.R. Of course, the deconvolution step will work for any set of mutational signatures obtained after an extraction, not just via hdp. 

### 3. Re-fitting signatures to counts.

Lastly, to more accurately assess the exposures of each mutational signature in the final set, we used the R package sigfit (https://github.com/kgori/sigfit). This step is carried out as part of the figure generation and its precise application can be found in the appropriate R scripts. 



