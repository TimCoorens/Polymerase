# HDP Flow
# Tim Coorens, Feb 2020

# PART 1: Generate input matrix

options(stringsAsFactors = F)
library(data.table)

patients=read.table("patients.txt")[,1]
all_muts=c()
for (p in patients){
  if(file.exists(paste0(p,"/snp_assigned_to_branches.txt"))){
    all_muts=rbind(all_muts,fread(paste0(p,"/snp_assigned_to_branches.txt"),header=T,sep="\t",data.table=F))
    print(p)
  }
}

mutlist_to_96_contexts=function(mutlist){
  genomeFile = "/nfs/cancer_ref01/Homo_sapiens/37/genome.fa"
  library("GenomicRanges")
  library("Rsamtools")
  library("MASS")
  samples=unique(mutlist$SampleID)
  trinuc_mut_mat=matrix(0,ncol=96,nrow=length(samples))
  for (n in 1:length(samples)){
    s=samples[n]
    mutations=as.data.frame(mutlist[mutlist$SampleID==s,c("Chr","Pos","Ref","Alt")])
    colnames(mutations) = c("chr","pos","ref","mut")
    mutations$pos=as.numeric(mutations$pos)
    mutations = mutations[(mutations$ref %in% c("A","C","G","T")) & (mutations$mut %in% c("A","C","G","T")) & mutations$chr %in% c(1:22,"X","Y"),]
    mutations$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(mutations$chr, IRanges(as.numeric(mutations$pos)-1, 
                                                                                       as.numeric(mutations$pos)+1))))
    ntcomp = c(T="A",G="C",C="G",A="T")
    mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
    mutations$trinuc_ref_py = mutations$trinuc_ref
    for (j in 1:nrow(mutations)) {
      if (mutations$ref[j] %in% c("A","G")) { # Purine base
        mutations$sub[j] = paste(ntcomp[mutations$ref[j]],ntcomp[mutations$mut[j]],sep=">")
        mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
      }
    }
    freqs = table(paste(mutations$sub,paste(substr(mutations$trinuc_ref_py,1,1),substr(mutations$trinuc_ref_py,3,3),sep="-"),sep=","))
    sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
    ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
    full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
    freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
    trinuc_mut_mat[n,]=freqs_full
    print(s)
  }
  colnames(trinuc_mut_mat)=full_vec
  rownames(trinuc_mut_mat)=samples
  return(trinuc_mut_mat)
}

# Start with a list of mutations per sample (or per branch, if phylogeny)
# make sure colnames are Chr Pos Ref Alt SampleID (order shouldn't matter)

## If required/desired, subset max. mutations per sample here
freqs=table(all_muts$SampleID)
samples=names(freqs)
max_mut_num=3000
mut_list_subsampled=c()
for(s in samples){
  mut_list_sub=all_muts[all_muts$SampleID==s,]
  if(nrow(mut_list_sub)<=max_mut_num){
    mut_list_subsampled=rbind(mut_list_subsampled,mut_list_sub)
  }else{
    random_select=sample(size=max_mut_num,x=1:nrow(mut_list_sub),replace=F)
    mut_list_subsampled=rbind(mut_list_subsampled,mut_list_sub[random_select,])
  }
  print(s)
}
##

trinuc_mut_mat=mutlist_to_96_contexts(all_muts[,c(1:4,7)])
#if subsampled:
trinuc_mut_mat=mutlist_to_96_contexts(mut_list_subsampled[,c(1:4,7)])

samples=rownames(trinuc_mut_mat)
key_table=data.frame(Sample=samples,
                     Patient=substr(samples,1,7))
write.table(trinuc_mut_mat,"trinuc_mut_mat_2020_03_17.txt")
write.table(key_table,"key_table_2020_03_17.txt")

## PART 2: Run multiple HDP chains

# I run this using the following simple bash script:

#for n in $(seq 1 20);
#do
#bsub -o $PWD/log.%J -e $PWD/err.%J -q normal -R'select[mem>10000] rusage[mem=10000]' -M10000 -J $n Rscript hdp_single_chain.R $n
#done


###

options(stringsAsFactors = F)
library(hdp)
lower_threshold=100

n=as.numeric(commandArgs(T)[1])
mutations=read.table("trinuc_mut_mat.txt")
key_table=read.table("key_table.txt")

#If requiring a minimum number of mutations:
sample_remove=rownames(mutations)[rowSums(mutations)<lower_threshold]
mutations=mutations[!rownames(mutations)%in%sample_remove,]
key_table=key_table[!key_table$Sample%in%sample_remove,]

#Hierarchy is set per patient, can change if wanted
freq=table(key_table$Patient)

hdp_mut <- hdp_init(ppindex = c(0, rep(1,length(freq)),rep(2:(length(freq)+1), times=freq)), # index of parental node
                    cpindex = c(1, rep(2,length(freq)),rep(3:(length(freq)+2), times=freq)), # index of the CP to use
                    hh = rep(1, 96), # prior is uniform over 96 categories
                    alphaa = rep(1,length(freq)+2), # shape hyperparameters for 2 CPs
                    alphab = rep(1,length(freq)+2))  # rate hyperparameters for 2 CPs

hdp_mut <- hdp_setdata(hdp_mut, 
                       dpindex = (length(freq)+2):numdp(hdp_mut), # index of nodes to add data to
                       mutations)

hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=10,seed=n*300)

chain=hdp_posterior(hdp_activated,
                    burnin=20000,
                    n=100,
                    seed=n*1000,
                    space=200,
                    cpiter=3)
saveRDS(chain,paste0("hdp_chain_",n,".Rdata"))

### PART 3: combine the results
options(stringsAsFactors = F)
library(hdp)

chlist <- vector("list", 20)
for (i in 1:20){
  if(file.exists(paste0("hdp_chain_",i,".Rdata"))){
    chlist[[i]] <- readRDS(paste0("hdp_chain_",i,".Rdata"))
  }
}
if(any(unlist(lapply(chlist,is.null)))) chlist=chlist[-which(unlist(lapply(chlist,is.null)))]

mut_example_multi <- hdp_multi_chain(chlist)
pdf("QC_plots_chain.pdf") 
par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
p1 <- lapply(chains(mut_example_multi), plot_lik, bty="L", start=1000)
p2 <- lapply(chains(mut_example_multi), plot_numcluster, bty="L")
p3 <- lapply(chains(mut_example_multi), plot_data_assigned, bty="L")
dev.off()

mut_example_multi <- hdp_extract_components(mut_example_multi) #This step can take a while. If too long, submit R script as job
saveRDS(mut_example_multi,"HDP_multi_chain.Rdata")

pdf("muts_attributed.pdf")
plot_comp_size(mut_example_multi, bty="L")
dev.off()

trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              each=16))
mut_colours=c("dodgerblue","black","red","grey70","olivedrab3","plum2")

#dev.new(width=12,height=4)
#par(mfrow=c(3,4))


for (i in 0:mut_example_multi@numcomp){
  pdf(paste0("hdp_component_",i,".pdf"),width=12,height=4)
  
  plot_comp_distn(mut_example_multi, cat_names=trinuc_context,
                  grouping=group_factor, col=mut_colours,comp=i,
                  col_nonsig="grey80", show_group_labels=TRUE)
  dev.off()
}

plot_dp_comp_exposure(mut_example_multi,
                      dpindices=2:4, incl_numdata_plot=FALSE,
                      col=c(RColorBrewer::brewer.pal(12, "Paired"),"magenta","firebrick"),
                      incl_nonsig=TRUE, cex.names=0.8,
                      ylab_exp = 'Signature exposure', leg.title = 'Signature')

pdf("signature_attribution.pdf",width=10,height=8)

#key_table=read.table("key_table.txt")

plot_dp_comp_exposure(mut_example_multi, dpindices=(length(unique(key_table$Patient))+2):length(mut_example_multi@comp_dp_counts), incl_nonsig = T,ylab_exp = 'Signature exposure', leg.title = 'Signature',
                      col=c(RColorBrewer::brewer.pal(12, "Set3"),"magenta","firebrick"))
dev.off()

mean_assignment=as.data.frame(comp_dp_distn(mut_example_multi)$mean)
write.table(mean_assignment,"mean_assignment_hdp.txt")
mean_sigs=as.data.frame(t(comp_categ_distn(mut_example_multi)$mean))
write.table(mean_sigs,"hdp_sigs.txt")


