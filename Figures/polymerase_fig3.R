#Polymerase paper, Figure 3

#---------
#Duplex sequencing blood and sperm
options(stringsAsFactors = F)
library(sigfit)
data("cosmic_signatures_v2")
final_sigs=t(read.table("final_sigs_2020_04_21.txt"))
counts=read.table("dupseq_trinuc_counts.txt")
trinuc_factors=read.table("trinuc_factors_dupseq.txt")
trinuc_freqs=read.table("trinuc_freqs_dupseq.txt")

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

POLE_dupseq_fit=fit_signatures(counts=counts_adjusted[c(1,3),], 
                               signatures = final_sigs[c("SBS1","SBS5","SBS10a","SBS10B","SBS28"),],
                               iter = 20000, 
                               warmup = 10000,
                               model="poisson",
                               chains = 2)
POLE_dupseq_exposures <- retrieve_pars(POLE_dupseq_fit, 
                                       par = "exposures", 
                                       hpd_prob = 0.95)
write.table(POLE_botseq_exposures$mean,"POLE_dupseq_exposures_mean.txt",sep="\t",quote=F)

POLD_dupseq_fit=fit_signatures(counts=counts[c(2,4),], 
                               signatures = final_sigs[c("SBS1","SBS5","SBS10C","SBS28"),],
                               iter = 20000, 
                               warmup = 10000,
                               model="poisson",
                               chains = 2)
POLD_dupseq_exposures <- retrieve_pars(POLD_dupseq_fit, 
                                       par = "exposures", 
                                       hpd_prob = 0.95)
write.table(POLD_dupseq_exposures$mean,"POLE_dupseq_exposures_mean.txt",sep="\t",quote=F)



#Endometrium
endo_counts=read.table("endo_sigs_input.txt")
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
  counts=read.table(paste0(b,"_per_sample_sig_input.txt"))
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
pdf(width=10,height=6,"sigfit_all_tissues_restricted.pdf")
plot(0,type = "n", xaxt = "n", yaxt = "n", xlab = "", bty="n",
     ylab = "",xlim=c(-1,20),ylim=c(-0.55,1.1))
segments(x0=0,y0=0,y1=1)
segments(x0=0,x1=-0.1,y0=seq(0,1,length.out = 5))
text(x=-0.1,y=seq(0,1,length.out = 5),labels=seq(0,1,by=0.25),pos=2)
n=1

m=n
rect(xleft=0.5*(n+0.1),xright=0.5*(n+0.9),ybottom=0,ytop=sum(POLE_dupseq_exposures$mean[1,c("SBS1","SBS5")]),col="red3")
rect(xleft=0.5*(n+0.1),xright=0.5*(n+0.9),ytop=1,ybottom=sum(POLE_dupseq_exposures$mean[1,c("SBS1","SBS5")]),col="steelblue2")
n=n+1
rect(xleft=0.5*(n+0.1),xright=0.5*(n+0.9),ybottom=0,ytop=sum(POLD_dupseq_exposures$mean[1,c("SBS1","SBS5")]),col="red3")
rect(xleft=0.5*(n+0.1),xright=0.5*(n+0.9),ytop=1,ybottom=sum(POLD_dupseq_exposures$mean[1,c("SBS1","SBS5")]),col="green3")
text(x=(m+n)/4+0.5,y=-0.03,label="blood",pos=2,srt=90)
n=n+2

m=n
rect(xleft=0.5*(n+0.1),xright=0.5*(n+0.9),ybottom=0,ytop=sum(POLE_dupseq_exposures$mean[2,c("SBS1","SBS5")]),col="red3")
rect(xleft=0.5*(n+0.1),xright=0.5*(n+0.9),ytop=1,ybottom=sum(POLE_dupseq_exposures$mean[2,c("SBS1","SBS5")]),col="steelblue2")
n=n+1
rect(xleft=0.5*(n+0.1),xright=0.5*(n+0.9),ybottom=0,ytop=sum(POLD_dupseq_exposures$mean[2,c("SBS1","SBS5")]),col="red3")
rect(xleft=0.5*(n+0.1),xright=0.5*(n+0.9),ytop=1,ybottom=sum(POLD_dupseq_exposures$mean[2,c("SBS1","SBS5")]),col="green3")
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
pdf("sperm_dupseq_sigfit_rates.pdf",useDingbats = F,height=6,width=4)
plot(0,type = "n", xaxt = "n", yaxt = "n", xlab = "", bty="n",
     ylab = "",xlim=c(2.5,7),ylim=c(-2,ymax*1.2))
segments(x0=3,y0=0,y1=ymax)
segments(x0=3,y0=0,x1=6)
segments(x0=3,x1=2.9,y0=seq(0,ymax,by=5))
text(x=2.9,y=seq(0,ymax,by=5),labels=seq(0,ymax,by=5),pos=2)
points(x=4:6-0.5,y=c(POLE_sperm_rate,POLD1_sperm_rate,normal_rate),
       cex=1.3,pch=16)
segments(x0=5.5,y0=normal_rate_lo,y1=normal_rate_hi,lwd=2)
segments(x0=3.5,y0=sum(POLE_dupseq_exposures$upper_95[2,c("SBS1","SBS5")])*POLE_sperm_rate,
         y1=sum(POLE_dupseq_exposures$lower_95[2,c("SBS1","SBS5")])*POLE_sperm_rate,lwd=2,col='grey70')
segments(x0=4.5,y0=sum(POLD_dupseq_exposures$upper_95[2,c("SBS1","SBS5")])*POLD1_sperm_rate,
         y1=sum(POLD_dupseq_exposures$lower_95[2,c("SBS1","SBS5")])*POLD1_sperm_rate,lwd=2,col='grey70')

points(x=c(3.5,4.5),y=c(sum(POLE_dupseq_exposures$mean[2,c("SBS1","SBS5")])*POLE_sperm_rate,
                        sum(POLD_dupseq_exposures$mean[2,c("SBS1","SBS5")])*POLD1_sperm_rate),
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

pdf("blood_dupseq_sigfit_rates.pdf",useDingbats = F,height=6,width=4)
plot(0,type = "n", xaxt = "n", yaxt = "n", xlab = "", bty="n",
     ylab = "",xlim=c(2.5,7),ylim=c(-10,ymax*1.2))

segments(x0=3,y0=0,y1=ymax)
segments(x0=3,y0=0,x1=6)
segments(x0=3,x1=2.9,y0=seq(0,ymax,by=25))
text(x=2.9,y=seq(0,ymax,by=25),labels=seq(0,ymax,by=25),pos=2)
points(x=4:6-0.5,y=c(POLE_blood_rate,POLD1_blood_rate,normal_rate),
       cex=1.3,pch=16)
segments(x0=5.5,y0=normal_rate_lo,y1=normal_rate_hi,lwd=2)
segments(x0=3.5,y0=sum(POLE_dupseq_exposures$upper_95[1,c("SBS1","SBS5")])*POLE_blood_rate,
         y1=sum(POLE_botseq_exposures$lower_95[1,c("SBS1","SBS5")])*POLE_blood_rate,lwd=2,col='grey70')
segments(x0=4.5,y0=sum(POLE_dupseq_exposures$upper_95[1,c("SBS1","SBS5")])*POLD1_blood_rate,
         y1=sum(POLD_dupseq_exposures$lower_95[1,c("SBS1","SBS5")])*POLD1_blood_rate,lwd=2,col='grey70')

points(x=c(3.5,4.5),y=c(sum(POLE_dupseq_exposures$mean[1,c("SBS1","SBS5")])*POLE_blood_rate,
                        sum(POLD_dupseq_exposures$mean[1,c("SBS1","SBS5")])*POLD1_blood_rate),
       cex=1.3,pch=16,col="grey70")

text(y=0,x=4:6-0.5,labels=c("POLE L424V","POLD1 S478N","Normal"),pos=1,srt=90)
dev.off()

ymax=300
pdf("endo_sigfit_rates.pdf",useDingbats = F,height=6,width=4)
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
library(sigfit)
final_sigs=t(read.table("final_sigs_2020_04_21.txt"))
embryonic_counts=read.table("sigs_input_embryonic_onlyVAF015.txt")
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
pdf("embryonic_POLE_muts.pdf",useDingbats = F,width=8,height=5)
barplot(t(as.matrix(embryonic_sig_counts_POLE)),las=2,col=c("red3","steelblue2"),ylim=c(0,150))
dev.off()

#POLD embryonic insertions - Fig 3E

counts=read.table("POLD_insertion_T_homopolymer_counts.txt")

pdf("embryonic_POLD_ins_final.pdf",useDingbats = F,width=8,height=5)
barplot(counts$V2,las=2,col="grey50",ylim=c(0,10),names.arg = counts$V1)
dev.off()
