#Deconvolute HDP signatures into reference signatures and arrive at final set
options(stringsAsFactors = F)
library(hdp)
library(RColorBrewer)
library(lsa)
library(lattice)

mut.cols = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)

#Load HDP signatures
hdp_sigs=read.table("hdp_sigs.txt")
#Load reference signatures
colon_sigs=read.table("Signatures_trinucleotide_composition_SBS_ABCD.txt",header=T) #Reference for SBSA-D
colon_sigs=colon_sigs[,grepl("N",colnames(colon_sigs))]
colnames(colon_sigs)=c("SBSC","SBSA","SBSB","SBSD")
ref=cbind(read.csv("PCAWG_sigProfiler_SBS_signatures.csv", header=T, stringsAsFactors = F, row.names = 1),colon_sigs)

#Assess cosine similarities for all reference signatures
cosine_matrix=data.frame(matrix(nrow=ncol(hdp_sigs), ncol=ncol(ref)))
rownames(cosine_matrix)=colnames(hdp_sigs)
colnames(cosine_matrix)=colnames(ref)

for (n in 1:nrow(cosine_matrix)) {
  for (m in 1:ncol(cosine_matrix)) {
    cosine_matrix[n,m] <- cosine(x=hdp_sigs[,rownames(cosine_matrix)[n]],
                                 y=ref[,colnames(cosine_matrix)[m]])
  }
}

write.table(cosine_matrix, "Cosine_similarities.txt",sep="\t",quote=F)

#plot output
pdf("HDP_2020_03/cosine_similarities.pdf", height=5, width=15)
color.palette = colorRampPalette(c("white", "orange", "purple"))
levelplot(t(cosine_matrix[dim(cosine_matrix)[1]:1,]),col.regions=color.palette, aspect="fill", scales=list(x=list(rot=90)))
dev.off()

#select which signatures to decompose into reference signatures
sigs_to_decompose=rowSums(cosine_matrix>0.85)
for(n in 1:nrow(cosine_matrix)){
  print(paste0(rownames(cosine_matrix)[n],": ",paste(colnames(cosine_matrix)[cosine_matrix[n,]>0.7],collapse=",")))
}

#First iteration; decomposed hdp sigs into all suspected sigs 
gdsigs=c("SBS1", "SBS5", "SBS7a","SBS7b", "SBS10a","SBS10b", "SBS17a", "SBS17b", "SBS18", "SBS25", "SBS28", "SBS31","SBS35","SBSA","SBSB") 

signatures=t(ref[,gdsigs])
sample_list=paste0("N",c(0:1,3:6,8:12)) #N2, N7 and N13 are new polymerase signatures 10c,10d,10e - not decomposed
profiles=hdp_sigs[,sample_list]

signature_fraction = matrix(NA,nrow=nrow(signatures),ncol=length(sample_list))
rownames(signature_fraction) = rownames(signatures)
colnames(signature_fraction) = sample_list
maxiter <- 1000

for (j in 1:length(sample_list)) {
  freqs = profiles[,j]
  freqs[is.na(freqs)] = 0
  # EM algowith to estimate the signature contribution
  alpha = runif(nrow(signatures)); alpha=alpha/sum(alpha) # Random start (seems to give ~identical results)
  for (iter in 1:maxiter) {
    contr = t(array(alpha,dim=c(nrow(signatures),96))) * t(signatures)
    probs = contr/array(rowSums(contr),dim=dim(contr))
    probs = probs * freqs
    old_alpha = alpha
    alpha = colSums(probs)/sum(probs)
    if (sum(abs(alpha-old_alpha))<1e-5) {
      break
    }
  }
  # Saving the signature contributions for the sample
  print(j/length(sample_list))
  signature_fraction[,j] = alpha
}

#Plot initial deconvolution and save results
pdf("Deconvolution_hdp_sigs_R1.pdf", height=5, width=10)
dev.new(height=5, width=10)
color.palette = colorRampPalette(c("white", "orange", "purple"))
levelplot((signature_fraction[nrow(signature_fraction):1,]),col.regions=color.palette, aspect="fill", scales=list(x=list(rot=90)))
dev.off()

write.table(signature_fraction, "HDP_2020_03/hdp_known_sigs_broken_down_into_pcawg_gd_sigs.txt", sep="\t", col.names=T, row.names = T, quote=F)

sigs_deconv_R2=list()
for(n in 1:length(sample_list)){
  sigs_deconv_R2[[n]]=rownames(signature_fraction)[signature_fraction[,n]>0.15]
}
names(sigs_deconv_R2)=colnames(signature_fraction)

#Some signatures only have one major contributor - leaving those out
#N0 -> SBS5
#N8 -> SBSA
#N11 -> SBSB
sigs_to_deconv=names(sigs_deconv_R2)[unlist(lapply(sigs_deconv_R2,length))>1]

ref_sigs_R2=sort(unique(unlist(sigs_deconv_R2)))
signature_fractionR2=matrix(NA,ncol=length(sigs_to_deconv),nrow=length(ref_sigs_R2))
rownames(signature_fractionR2)=ref_sigs_R2
colnames(signature_fractionR2)=sigs_to_deconv
#repeat the deconvolution with the identified constitutive signatures
n=1
for(s in sigs_to_deconv){
  gdsigs <- sigs_deconv_R2[[s]]
  signatures <- t(ref[,gdsigs])
  
  signature_fraction = matrix(NA,nrow=nrow(signatures),ncol=length(sample_list))
  rownames(signature_fraction) = rownames(signatures)
  colnames(signature_fraction) = sample_list
  maxiter <- 1000
  
  freqs = profiles[,s]
  freqs[is.na(freqs)] = 0
  
  alpha = runif(nrow(signatures)); alpha=alpha/sum(alpha) # Random start (seems to give ~identical results)
  for (iter in 1:maxiter) {
    contr = t(array(alpha,dim=c(nrow(signatures),96))) * t(signatures)
    probs = contr/array(rowSums(contr),dim=dim(contr))
    probs = probs * freqs
    old_alpha = alpha
    alpha = colSums(probs)/sum(probs)
    if (sum(abs(alpha-old_alpha))<1e-5) {
      break
    }
  }
  # Saving the signature contributions for the sample
  signature_fractionR2[gdsigs,n] = alpha
  n=n+1
  reconsbs <- rep(0,96)
  for (g in gdsigs) {
    reconsbs=reconsbs+(ref[,g]*alpha[g])
  }
  cosine_reconst=cosine(x=reconsbs, y=hdp_sigs[,s])
  print(paste0(s,": ",cosine_reconst))
  pdf(paste0("HDP_2020_03/HDP_",s,"_reconstitution2.pdf"))
  par(mfrow=c(length(alpha)+2,1))
  par(mar=c(1,2,4,1))
  barplot(hdp_sigs[,s], col=mut.cols, main=paste0("HDP ",s),names.arg="")
  barplot(reconsbs, col=mut.cols, main=paste0("Reconstituted ",s," cosine similarity to original: ", round(cosine_reconst, digits=2)))
  for (g in gdsigs) {
    barplot(ref[,g], col=mut.cols, main=paste0("PCAWG ", g, " accounts for ", round(alpha[g], digits=2)))
  }
  dev.off()
}

sigs_deconv_R2$N2="SBS10C"
sigs_deconv_R2$N7="SBS10D"
sigs_deconv_R2$N13="SBS10E"
sigs_deconv_R2$N10="SBS35like"

saveRDS(sigs_deconv_R2,"hdp2refsigs.Rdata")


#Combine hdp signatures that did not get deconvolved and reference signatures into final table
final_sigs=cbind(sbs10b_prime,hdp_sigs[,c("N2","N7","N13","N10")],ref[,ref_sigs_R2])
#Rename the HDP components that didn't get deconvoluted
colnames(final_sigs)[1:4]=paste0("SBS10",LETTERS[2:5])
colnames(final_sigs)[5]="SBS35like"
#Order them - optional
final_sigs=final_sigs[,paste0("SBS",c(1,5,"7a","7b","10a","10B","10C","10D","10E","17b",28,"35like","A","B"))]
write.table(final_sigs,"final_sigs_2020_04_21.txt")
