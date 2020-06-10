#Polymerase paper, Figure 3

#---------
#Botseq
options(stringsAsFactors = F)
library(sigfit)
data("cosmic_signatures_v2")
final_sigs=t(read.table("~/Desktop/Phil_POL/final_sigs_2020_04_21.txt"))
counts=read.table("~/Desktop/Phil_POL/botseq_trinuc_counts.txt")
trinuc_factors=read.table("~/Desktop/Phil_POL/trinuc_factors_botseq.txt")
trinuc_freqs=read.table("~/Desktop/Phil_POL/trinuc_freqs_botseq.txt")

colnames(counts)=colnames(cosmic_signatures_v2)

trinucs_96=substr(colnames(cosmic_signatures_v2),1,3)
trinuc_freqs_all=trinuc_factors_all=matrix(nrow=4,ncol=96)
for(i in 1:4){
  trinuc_freqs_all[i,]=as.numeric(trinuc_freqs[trinucs_96,i])
  trinuc_factors_all[i,]=as.numeric(trinuc_factors[trinucs_96,i])
}
colnames(trinuc_freqs_all)=colnames(final_sigs)=colnames(cosmic_signatures_v2)
rownames(trinuc_freqs_all)=rownames(counts)
counts_adjusted=round(counts*trinuc_factors_all)

POLE_botseq_fit=fit_signatures(counts=counts_adjusted[c(1,3),], 
                               signatures = final_sigs[c("SBS1","SBS5","SBS10a","SBS10B","SBS28"),],
                               iter = 20000, 
                               warmup = 10000,
                               model="poisson",
                               chains = 2)
POLE_botseq_exposures <- retrieve_pars(POLE_botseq_fit, 
                                       par = "exposures", 
                                       hpd_prob = 0.95)
POLD_botseq_fit=fit_signatures(counts=counts[c(2,4),], 
                               signatures = final_sigs[c("SBS1","SBS5","SBS10C","SBS28"),],
                               iter = 20000, 
                               warmup = 10000,
                               model="poisson",
                               chains = 2)
POLD_botseq_exposures <- retrieve_pars(POLD_botseq_fit, 
                                       par = "exposures", 
                                       hpd_prob = 0.95)


#Endometrium
endo_counts=read.table("~/Desktop/Phil_POL/endo_sigs_input.txt")
colnames(endo_counts)=colnames(cosmic_signatures_v2)

endo_fit=fit_signatures(counts=endo_counts, 
                        signatures = final_sigs[c("SBS1","SBS5","SBS10a","SBS10B","SBS28"),],
                        iter = 20000, 
                        warmup = 10000,
                        model="poisson",
                        chains = 2)
endo_exposures=retrieve_pars(endo_fit, 
                             par = "exposures", 
                             hpd_prob = 0.95)

#Autopsy patient
patient="PD44594"
blocks=paste0("PD44594",c("c","f","h"))
PD44594_pars=list()
for (b in blocks){
  counts=read.table(paste0("~/Desktop/Phil_POL/",b,"_per_sample_sig_input.txt"))
  colnames(counts)=colnames(cosmic_signatures_v2)
  sigs_select=c("SBS1","SBS5","SBS10a","SBS10B","SBS10E","SBS28")
  if(b=="PD44594f")sigs_select=c(sigs_select,"SBS7a","SBS7b")
  
  fit=fit_signatures(counts=counts, 
                          signatures = final_sigs[sigs_select,],
                          iter = 20000, 
                          warmup = 10000,
                          model="poisson",
                          chains = 2)
  PD44594_pars[[b]]=retrieve_pars(fit, 
                               par = "exposures", 
                               hpd_prob = 0.95)
}
names(tissues)=blocks



#Fig 3A
pdf(width=10,height=6,"~/Desktop/Phil_POL/sigfit_all_tissues_restricted.pdf")
plot(0,type = "n", xaxt = "n", yaxt = "n", xlab = "", bty="n",
     ylab = "",xlim=c(-1,20),ylim=c(-0.55,1.1))
segments(x0=0,y0=0,y1=1)
segments(x0=0,x1=-0.1,y0=seq(0,1,length.out = 5))
text(x=-0.1,y=seq(0,1,length.out = 5),labels=seq(0,1,by=0.25),pos=2)
n=1

m=n
rect(xleft=0.5*(n+0.1),xright=0.5*(n+0.9),ybottom=0,ytop=sum(POLE_botseq_exposures$mean[1,c("SBS1","SBS5")]),col="red3")
rect(xleft=0.5*(n+0.1),xright=0.5*(n+0.9),ytop=1,ybottom=sum(POLE_botseq_exposures$mean[1,c("SBS1","SBS5")]),col="steelblue2")
n=n+1
rect(xleft=0.5*(n+0.1),xright=0.5*(n+0.9),ybottom=0,ytop=sum(POLD_botseq_exposures$mean[1,c("SBS1","SBS5")]),col="red3")
rect(xleft=0.5*(n+0.1),xright=0.5*(n+0.9),ytop=1,ybottom=sum(POLD_botseq_exposures$mean[1,c("SBS1","SBS5")]),col="green3")
text(x=(m+n)/4+0.5,y=-0.03,label="blood",pos=2,srt=90)
n=n+2

m=n
rect(xleft=0.5*(n+0.1),xright=0.5*(n+0.9),ybottom=0,ytop=sum(POLE_botseq_exposures$mean[2,c("SBS1","SBS5")]),col="red3")
rect(xleft=0.5*(n+0.1),xright=0.5*(n+0.9),ytop=1,ybottom=sum(POLE_botseq_exposures$mean[2,c("SBS1","SBS5")]),col="steelblue2")
n=n+1
rect(xleft=0.5*(n+0.1),xright=0.5*(n+0.9),ybottom=0,ytop=sum(POLD_botseq_exposures$mean[2,c("SBS1","SBS5")]),col="red3")
rect(xleft=0.5*(n+0.1),xright=0.5*(n+0.9),ytop=1,ybottom=sum(POLD_botseq_exposures$mean[2,c("SBS1","SBS5")]),col="green3")
text(x=(m+n)/4+0.5,y=-0.03,label="sperm",pos=2,srt=90)
n=n+2

m=n
for(i in 1:nrow(endo_counts)){
  rect(xleft=0.5*(n+0.1),xright=0.5*(n+0.9),ybottom=0,ytop=sum(endo_exposures$mean[i,c("SBS1","SBS5")]),col="red3")
  rect(xleft=0.5*(n+0.1),xright=0.5*(n+0.9),ytop=1,ybottom=sum(endo_exposures$mean[i,c("SBS1","SBS5")]),col="steelblue2")
  n=n+1
}
text(x=(m+n)/4+0.5,y=-0.03,label="endometrium",pos=2,srt=90)
n=n+1
for(b in blocks[blocks%in%paste0("PD44594",c("c","f","h"))]){
  pars=PD44594_pars[[b]]
  m=n
  for(s in 1:nrow(pars$mean)){
    start=0
    end=sum(pars$mean[s,c("SBS1","SBS5")])
    rect(xleft=0.5*(n+0.1),xright=0.5*(n+0.9),ybottom=start,ytop=end,col="red3")
    start=end
    if(b=="PD44594f"){
      end=start+sum(pars$mean[s,c("SBS7a","SBS7b")])
      rect(xleft=0.5*(n+0.1),xright=0.5*(n+0.9),ybottom=start,ytop=end,col="plum2")
      start=end
    }
    end=start+sum(pars$mean[s,c("SBS10a","SBS10B","SBS28")])
    rect(xleft=0.5*(n+0.1),xright=0.5*(n+0.9),ybottom=start,ytop=end,col="steelblue2")
    start=end
    end=1
    rect(xleft=0.5*(n+0.1),xright=0.5*(n+0.9),ybottom=start,ytop=end,col="navy")
    n=n+1
  }
  text(x=(m+n)/4+0.5,y=-0.03,label=tissues[b],pos=2,srt=90)
  n=n+1
}
dev.off()

#Fig 3B
#Paternal age effect numbers (Rahbari et al, 2016)
normal_rate=2.87
normal_rate_hi=3.64
normal_rate_lo=2.11
ymax=25
pdf("~/Desktop/Phil_POL/sperm_nanoseq_sigfit_rates.pdf",useDingbats = F,height=6,width=4)
plot(0,type = "n", xaxt = "n", yaxt = "n", xlab = "", bty="n",
     ylab = "",xlim=c(2.5,7),ylim=c(-2,ymax*1.2))
segments(x0=3,y0=0,y1=ymax)
segments(x0=3,y0=0,x1=6)
segments(x0=3,x1=2.9,y0=seq(0,ymax,by=5))
text(x=2.9,y=seq(0,ymax,by=5),labels=seq(0,ymax,by=5),pos=2)
points(x=4:6-0.5,y=c(POLE_sperm_rate,POLD1_sperm_rate,normal_rate),
       cex=1.3,pch=16)
segments(x0=5.5,y0=normal_rate_lo,y1=normal_rate_hi,lwd=2)
segments(x0=3.5,y0=sum(POLE_botseq_exposures$upper_95[2,c("SBS1","SBS5")])*POLE_sperm_rate,
         y1=sum(POLE_botseq_exposures$lower_95[2,c("SBS1","SBS5")])*POLE_sperm_rate,lwd=2,col='grey70')
segments(x0=4.5,y0=sum(POLD_botseq_exposures$upper_95[2,c("SBS1","SBS5")])*POLD1_sperm_rate,
         y1=sum(POLD_botseq_exposures$lower_95[2,c("SBS1","SBS5")])*POLD1_sperm_rate,lwd=2,col='grey70')

points(x=c(3.5,4.5),y=c(sum(POLE_botseq_exposures$mean[2,c("SBS1","SBS5")])*POLE_sperm_rate,
                        sum(POLD_botseq_exposures$mean[2,c("SBS1","SBS5")])*POLD1_sperm_rate),
       cex=1.3,pch=16,col="grey70")

text(y=0,x=4:6-0.5,labels=c("POLE L424V","POLD1 S478N","Normal"),pos=1,srt=90)
dev.off()


POLE_blood_rate=823.366012/3155133188*5.8E9/17
POLD1_blood_rate=1584.23248/2463758361*5.8E9/71
#Normal blood mutation rate (Osorio et al, 2018)
normal_rate=14.2
normal_rate_hi=22.4
normal_rate_lo=6.1
ymax=100

pdf("~/Desktop/Phil_POL/blood_nanoseq_sigfit_rates.pdf",useDingbats = F,height=6,width=4)
plot(0,type = "n", xaxt = "n", yaxt = "n", xlab = "", bty="n",
     ylab = "",xlim=c(2.5,7),ylim=c(-10,ymax*1.2))

segments(x0=3,y0=0,y1=ymax)
segments(x0=3,y0=0,x1=6)
segments(x0=3,x1=2.9,y0=seq(0,ymax,by=25))
text(x=2.9,y=seq(0,ymax,by=25),labels=seq(0,ymax,by=25),pos=2)
points(x=4:6-0.5,y=c(POLE_blood_rate,POLD1_blood_rate,normal_rate),
       cex=1.3,pch=16)
segments(x0=5.5,y0=normal_rate_lo,y1=normal_rate_hi,lwd=2)
segments(x0=3.5,y0=sum(POLE_botseq_exposures$upper_95[1,c("SBS1","SBS5")])*POLE_blood_rate,
         y1=sum(POLE_botseq_exposures$lower_95[1,c("SBS1","SBS5")])*POLE_blood_rate,lwd=2,col='grey70')
segments(x0=4.5,y0=sum(POLD_botseq_exposures$upper_95[1,c("SBS1","SBS5")])*POLD1_blood_rate,
         y1=sum(POLD_botseq_exposures$lower_95[1,c("SBS1","SBS5")])*POLD1_blood_rate,lwd=2,col='grey70')

points(x=c(3.5,4.5),y=c(sum(POLE_botseq_exposures$mean[1,c("SBS1","SBS5")])*POLE_blood_rate,
                        sum(POLD_botseq_exposures$mean[1,c("SBS1","SBS5")])*POLD1_blood_rate),
       cex=1.3,pch=16,col="grey70")

text(y=0,x=4:6-0.5,labels=c("POLE L424V","POLD1 S478N","Normal"),pos=1,srt=90)
dev.off()

ymax=300
pdf("~/Desktop/Phil_POL/endo_sigfit_rates.pdf",useDingbats = F,height=6,width=4)
plot(0,type = "n", xaxt = "n", yaxt = "n", xlab = "", bty="n",
     ylab = "",xlim=c(6.5,12.5),ylim=c(-15,ymax*1.1))
points(x=seq(8,10,length.out=11),y=rowSums(endo_counts)/47,pch=16,cex=1.3)
points(x=seq(8,10,length.out=11),y=rowSums(endo_counts)/47*rowSums(endo_exposures$mean[,c("SBS1","SBS5")]),pch=16,col='grey70',cex=1.3)
segments(x0=seq(8,10,length.out=11),y0=rowSums(endo_counts)/47*rowSums(endo_exposures$lower_95[,c("SBS1","SBS5")]),
         y1=rowSums(endo_counts)/47*rowSums(endo_exposures$upper_95[,c("SBS1","SBS5")]),lwd=1.3,col='grey70')

points(x=11,y=29,pch=16,cex=1.3)
segments(x0=11,y0=23,y1=34,lwd=2)

segments(x0=7.5,y0=0,y1=ymax)
segments(x0=7.5,y0=0,x1=11.5)
segments(x0=7.5,x1=7.3,y0=seq(0,ymax,by=50))
text(x=7.3,y=seq(0,ymax,by=50),labels=seq(0,ymax,by=50),pos=2)
text(y=0,x=c(8.5,10.5),labels=c("POLE\nL424V","Normal"),pos=1)
dev.off()


#POLE embryonic snvs - Fig 3D
# options(stringsAsFactors = F)
# sampleInfo=read.table("/lustre/scratch116/casm/cgp/users/pr10/POL/sampleinformation/POL_final_updated_20200111.txt",header=T,sep="\t")
# library(data.table)
# blood_samples=sampleInfo[sampleInfo$tissue=="blood","sample"]
# non_blood_samples=sampleInfo[sampleInfo$tissue!="blood","sample"]
# 
# blood_all_muts=c()
# patients=unique(substr(blood_samples,1,7))
# for (p in patients){
#   s=blood_samples[grepl(p,blood_samples)]
#   if(file.exists(paste0(p,"/snp_NR_filtered_all.txt"))){
#     NR_flt_2=read.table(paste0(p,"/snp_NR_filtered_all.txt"))
#     NV_flt_2=read.table(paste0(p,"/snp_NV_filtered_all.txt"))
#     NR_flt_2[NR_flt_2==0]=1
#     Genotype=NV_flt_2/NR_flt_2
#     if (length(s)>1){
#       present_in_blood=rowSums(NV_flt_2[,s]>1)==length(s)
#       other_samples=sampleInfo$sample[sampleInfo$patient==p&sampleInfo$tissue2!="blood"]
#       pregastrulation=rowSums(Genotype[present_in_blood,other_samples]>0.25)>0
#       muts=as.data.frame(matrix(ncol=4,unlist(strsplit(rownames(Genotype)[present_in_blood],split="_")),byrow = T))
#       vaf=rowMeans(Genotype[present_in_blood,s])
#       blood_all_muts=rbind(blood_all_muts,cbind(p,muts,vaf,pregastrulation))
#     }
#     if (length(s)==1){
#       present_in_blood=NV_flt_2[,s]>1
#       other_samples=sampleInfo$sample[sampleInfo$patient==p&sampleInfo$tissue2!="blood"]
#       pregastrulation=rowSums(Genotype[present_in_blood,other_samples]>0.25)>0
#       muts=as.data.frame(matrix(ncol=4,unlist(strsplit(rownames(Genotype)[present_in_blood],split="_")),byrow = T))
#       vaf=Genotype[present_in_blood,s]
#       blood_all_muts=rbind(blood_all_muts,cbind(p,muts,vaf,pregastrulation))
#     }
#   }
#   print(s)
# }
# table(blood_all_muts$p,blood_all_muts$pregastrulation)
# blood_all_muts=blood_all_muts[!(blood_all_muts$p=="PD44587"&blood_all_muts$V1=="X"),]
# POLE=paste0("PD445",c(80,86,87,89,92,93))
# bad=read.table("jbrowse_embryonic_bad.txt")
# blood_all_muts=blood_all_muts[!blood_all_muts$V2%in%(bad$V2+49),]
# 
# colnames(blood_all_muts)=c("sample","chr","pos","ref","alt","vaf","pregastrulation")
# blood_all_muts$chr=paste0("chr",blood_all_muts$chr)
# blood_all_muts$pos=as.numeric(blood_all_muts$pos)
# select=blood_all_muts$vaf>0.15
# library(deconstructSigs)
# sigs.input=mut.to.sigs.input(mut.ref=blood_all_muts[select,c("sample","chr","pos","ref","alt")],
#                              sample.id = "sample",
#                              chr = "chr",
#                              pos = "pos",
#                              ref = "ref",
#                              alt = "alt")
# write.table(sigs.input,"sigs_input_embryonic_onlyVAF015.txt")

library(sigfit)
final_sigs=t(read.table("~/Desktop/Phil_POL/final_sigs_2020_04_21.txt"))
embryonic_counts=read.table("~/Desktop/Phil_POL/sigs_input_embryonic_onlyVAF015.txt")
colnames(final_sigs)=colnames(embryonic_counts)=colnames(sigs_to_fit)=colnames(cosmic_signatures_v2)

POLE=paste0("PD445",c(80,86,87,89,91:93))
POLD=rownames(embryonic_counts)[!rownames(embryonic_counts)%in%POLE]
embryonic_fit_POLE=fit_signatures(counts=embryonic_counts[rownames(embryonic_counts)%in%POLE,], 
                             signatures = final_sigs[c("SBS1","SBS5","SBS10a","SBS10B","SBS28"),],
                             iter = 20000, 
                             warmup = 10000,
                             model="poisson",
                             chains = 2)
embryonic_exposures_POLE=retrieve_pars(embryonic_fit_POLE, 
                                  par = "exposures", 
                                  hpd_prob = 0.95)
embryonic_sig_counts_POLE=data.frame(normal=rowSums(embryonic_exposures_POLE$mean[POLE,c("SBS1","SBS5")])*rowSums(embryonic_counts[POLE,]),
                                pol=rowSums(embryonic_exposures_POLE$mean[POLE,c("SBS10a","SBS10B","SBS28")])*rowSums(embryonic_counts[POLE,]))
pdf("~/Desktop/Phil_POL/embryonic_POLE_muts.pdf",useDingbats = F,width=8,height=5)
barplot(t(as.matrix(embryonic_sig_counts_POLE)),las=2,col=c("red3","steelblue2"),ylim=c(0,150))
dev.off()

#POLD embryonic insertions - Fig 3E
options(stringsAsFactors = F)
sampleInfo=read.table("/lustre/scratch116/casm/cgp/users/pr10/POL/sampleinformation/POL_final_updated_20200111.txt",header=T,sep="\t")
library(data.table)
blood_samples=sampleInfo[sampleInfo$tissue=="blood","sample"]
blood_all_muts=c()
patients=unique(substr(blood_samples,1,7))
patients=patients[patients!="PD44583"]
read_threshold=2
for (p in patients){
  path="/lustre/scratch116/casm/cgp/users/pr10/POL/indel_new/count_matrices_submission"
  s=blood_samples[grepl(p,blood_samples)]
  
  #if(file.exists(paste0(path,p,"/indel_NR_filtered_all.txt"))){
    NR_flt_2=read.table(paste0(p,"/indel_NR_filtered_all.txt"))
    NV_flt_2=read.table(paste0(p,"/indel_NV_filtered_all.txt"))
    NR_flt_2[NR_flt_2==0]=1
    Genotype=NV_flt_2/NR_flt_2
    if (length(s)>1){
      present_in_blood=rowSums(NV_flt_2[,s]>read_threshold)==length(s)
      indels_on_branches=fread(paste0(p,"/indel_assigned_to_snv_branches.txt"),data.table = F)
      indels_on_branches$ID=paste(indels_on_branches$Chr,indels_on_branches$Pos,indels_on_branches$Ref,indels_on_branches$Alt,sep="_")
      pregastrulation=rownames(NV_flt_2[present_in_blood,])%in%indels_on_branches$ID
      
      muts=as.data.frame(matrix(ncol=4,unlist(strsplit(rownames(Genotype)[present_in_blood],split="_")),byrow = T))
      vaf=rowMeans(Genotype[present_in_blood,s])
      NR=unlist(apply(NR_flt_2[present_in_blood,s],1, function(x) paste(x,collapse=";")))
      NV=unlist(apply(NV_flt_2[present_in_blood,s],1, function(x) paste(x,collapse=";")))
      blood_all_muts=rbind(blood_all_muts,cbind(p,muts,NV,NR,vaf,pregastrulation))
    }
    if (length(s)==1){
      present_in_blood=NV_flt_2[,s]>read_threshold
      indels_on_branches=fread(paste0(p,"/indel_assigned_to_snv_branches.txt"),data.table = F)
      indels_on_branches$ID=paste(indels_on_branches$Chr,indels_on_branches$Pos,indels_on_branches$Ref,indels_on_branches$Alt,sep="_")
      pregastrulation=rownames(NV_flt_2[present_in_blood,])%in%indels_on_branches$ID
      muts=as.data.frame(matrix(ncol=4,unlist(strsplit(rownames(Genotype)[present_in_blood],split="_")),byrow = T))
      vaf=Genotype[present_in_blood,s]
      NR=NR_flt_2[present_in_blood,s]
      NV=NV_flt_2[present_in_blood,s]
      
      blood_all_muts=rbind(blood_all_muts,cbind(p,muts,NV,NR,vaf,pregastrulation))
    }
  print(s)
}
table(blood_all_muts$p,blood_all_muts$pregastrulation)

sb_insertions=nchar(blood_all_muts$V4)==2&nchar(blood_all_muts$V3)==1
insertions_T=substr(blood_all_muts$V4,nchar(blood_all_muts$V4),nchar(blood_all_muts$V4))%in%c("A","T")
sb_ins_T=sb_insertions&insertions_T&blood_all_muts$V3%in%c("A","T")
blood_all_muts_flt=blood_all_muts[sb_ins_T,]

library("GenomicRanges")
library("Rsamtools")
library("MASS")
genomeFile = "/nfs/cancer_ref01/Homo_sapiens/37/genome.fa"

blood_all_muts_flt$preceding = as.vector(scanFa(genomeFile, GRanges(blood_all_muts_flt$V1, IRanges(as.numeric(blood_all_muts_flt$V2)-5, 
                                                                                                   as.numeric(blood_all_muts_flt$V2)-1))))

blood_all_muts_flt$following = as.vector(scanFa(genomeFile, GRanges(blood_all_muts_flt$V1, IRanges(as.numeric(blood_all_muts_flt$V2)+1, 
                                                                                                   as.numeric(blood_all_muts_flt$V2)+5))))
homopolymer_runs=blood_all_muts_flt$preceding=="TTTTT"|blood_all_muts_flt$following=="AAAAA"|blood_all_muts_flt$preceding=="AAAAA"|blood_all_muts_flt$following=="TTTTT"
thresh=0.1
table(blood_all_muts_flt$p[homopolymer_runs&blood_all_muts_flt$vaf>thresh],blood_all_muts_flt$pregastrulation[homopolymer_runs&blood_all_muts_flt$vaf>thresh])

homopolymer_runs=(substr(blood_all_muts_flt$V4,nchar(blood_all_muts_flt$V4),nchar(blood_all_muts_flt$V4))=="A"&(blood_all_muts_flt$following=="AAAAA"|blood_all_muts_flt$preceding=="AAAAA"))|
  (substr(blood_all_muts_flt$V4,nchar(blood_all_muts_flt$V4),nchar(blood_all_muts_flt$V4))=="T"&(blood_all_muts_flt$following=="TTTTT"|blood_all_muts_flt$preceding=="TTTTT"))

pold_samples=c(paste0("PD4458",c(1,2,4,5,8)),"PD44590")
thresh=0.1
counts=table(blood_all_muts_flt$p[blood_all_muts_flt$vaf>thresh])

pdf("embryonic_POLD_ins_final.pdf",useDingbats = F,width=8,height=5)
barplot(t(as.matrix(counts)),las=2,col="grey50",ylim=c(0,10))
dev.off()

