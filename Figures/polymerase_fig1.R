#Figure 1

options(stringsAsFactors = F)
data=read.table("POL_final_updated_20200422.txt",header=T,sep="\t")
data=data[!data$patient%in%c("PD40097","PD44594"),]
data$sensitivity=NA
for(k in 1:nrow(data)) data$sensitivity[k]=pbinom(q=3,size=data$MEDIAN_COVERAGE[k],p=data$median_vaf[k],lower.tail=F)
data$sbs_burden_corrected=data$sbstotal/data$sensitivity
data$indel_burden_corrected=data$indeltotal/data$sensitivity

data$full_descr=paste0(data$patient," (",data$age," years)")
select=data$type=="normal"
data_select=data[select,]
data_select$germline_mutation=factor(data_select$germline_mutation,levels=unique(data$germline_mutation),ordered=T)
data_select$full_descr=factor(data_select$full_descr,levels=unique(data_select$full_descr[order(data_select$germline_mutation,data_select$age)]),ordered=T)
colvec=c(rep("dodgerblue",7),rep("red",4),"limegreen","magenta")
pt_col="grey60"

#Fig1A
pdf("snv_burden_boxplot_all_2020_05.pdf",width=14,useDingbats = F)
boxplot(sbs_burden_corrected~full_descr,data=data_select,col=colvec,outline=F,las=2,ylim=c(-10000,40000),
        frame.plot=F,xaxt='n',yaxt='n',xlim=c(-1,14.5))
stripchart(sbs_burden_corrected ~ full_descr, data = data_select, method = "jitter",
           jitter=0.3,add = TRUE,vertical=T,pch=21,bg=pt_col,cex=2)
segments(x0=c(7.5,11.5,12.5),y0=0,y1=40000,lwd=2)
segments(y0=0,x0=0,x1=13.5)
segments(y0=-750,y1=0,x0=seq(1,14,1),lwd=1.2)
segments(y0=0,x0=0,y1=40000)
segments(x0=0,x1=-0.2,y0=seq(0,40000,10000),lwd=1.2)
text(labels=seq(0,40000,10000),y=seq(0,40000,10000),x=-0.3)
text(labels=levels(unique(data_select$full_descr)),y=-1000,x=seq(1,14,1),srt=45,pos=1)
dev.off()

#Fig1C
pdf("indel_burden_boxplot_all_2020_05.pdf",width=14,useDingbats = F)
boxplot(indel_burden_corrected~full_descr,data=data_select,col=colvec,outline=F,las=2,ylim=c(-1000,5000),
        frame.plot=F,xaxt='n',yaxt='n',xlim=c(-1,13.5))
stripchart(indel_burden_corrected ~ full_descr, data = data_select, method = "jitter",
           jitter=0.3,add = TRUE,vertical=T,pch=21,bg=pt_col,cex=2)
segments(x0=c(7.5,11.5,12.5),y0=0,y1=8000,lwd=2)
segments(y0=0,x0=0,x1=13.5)
segments(y0=-150,y1=0,x0=seq(1,14,1),lwd=1.2)
segments(y0=0,x0=0,y1=8000)
segments(x0=0,x1=-0.2,y0=seq(0,5000,1000),lwd=1.2)
text(labels=seq(0,5000,1000),y=seq(0,5000,1000),x=-0.3)
text(labels=levels(unique(data_select$full_descr)),y=-200,x=seq(1,14,1),srt=45,pos=1)
dev.off()

mut_burden_mean=unique(data_select[,c("full_descr","age","patient")])
mut_burden_mean=mut_burden_mean[order(mut_burden_mean$full_descr),]
mut_burden_mean$snv_mean=mut_burden_mean$indel_mean=NA
for(n in 1:nrow(mut_burden_mean)){
  mut_burden_mean$snv_mean[n]=mean(data_select$sbs_burden_corrected[data_select$patient==mut_burden_mean$patient[n]])
  mut_burden_mean$indel_mean[n]=mean(data_select$indel_burden_corrected[data_select$patient==mut_burden_mean$patient[n]])
}
normal_age_vec=seq(0,max(as.numeric(data$age)),by=0.05)
mut_rate=43.6

#Fig1B
pdf("snv_regression.pdf",useDingbats = F,width=4,height=4)
plot(y=mut_burden_mean$snv_mean,x=as.numeric(mut_burden_mean$age),pch=21,bg="white",col="white",xlim=c(0,75),ylim=c(0,27500),
     xlab="",ylab="Mean SBS burden")
lines(x=normal_age_vec,y=mut_rate*normal_age_vec,lwd=2,lty='dashed')
POLE=lm(snv_mean~age+0,data=mut_burden_mean[1:7,])
lines(x=normal_age_vec,y=POLE$coefficients*normal_age_vec,col='dodgerblue',lwd=3)
POLD1=lm(snv_mean~age+0,data=mut_burden_mean[8:11,])
lines(x=normal_age_vec,y=POLD1$coefficients*normal_age_vec,col='red',lwd=3)
POLD1_L=lm(snv_mean~age+0,data=mut_burden_mean[12,])
lines(x=normal_age_vec,y=POLD1_L$coefficients*normal_age_vec,col='limegreen',lwd=3)
POLD1_D=lm(snv_mean~age+0,data=mut_burden_mean[13,])
lines(x=normal_age_vec,y=POLD1_D$coefficients*normal_age_vec,col='magenta',lwd=3)
points(y=mut_burden_mean$snv_mean,x=as.numeric(mut_burden_mean$age),pch=21,bg=colvec,cex=2)
dev.off()

#Fig1C
pdf("indel_regression.pdf",useDingbats = F,width=4,height=4)
plot(y=mut_burden_mean$indel_mean,x=as.numeric(mut_burden_mean$age),pch=21,bg="white",col="white",xlim=c(0,75),ylim=c(0,4000),
     xlab="Age",ylab="Mean ID burden")
lines(x=normal_age_vec,y=mut_rate/10*normal_age_vec,lwd=2,lty='dashed')
POLE=lm(indel_mean~age+0,data=mut_burden_mean[1:7,])
lines(x=normal_age_vec,y=POLE$coefficients*normal_age_vec,col='dodgerblue',lwd=3)
POLD1=lm(indel_mean~age+0,data=mut_burden_mean[8:11,])
lines(x=normal_age_vec,y=POLD1$coefficients*normal_age_vec,col='red',lwd=3)
POLD1_L=lm(indel_mean~age+0,data=mut_burden_mean[12,])
lines(x=normal_age_vec,y=POLD1_L$coefficients*normal_age_vec,col='limegreen',lwd=3)
POLD1_D=lm(indel_mean~age+0,data=mut_burden_mean[13,])
lines(x=normal_age_vec,y=POLD1_D$coefficients*normal_age_vec,col='magenta',lwd=3)

points(y=mut_burden_mean$indel_mean,x=as.numeric(mut_burden_mean$age),pch=21,bg=colvec,cex=2)
dev.off()


normal_age_vec=seq(0,max(as.numeric(data$age)),by=0.05)
mut_rate=43.6

#Fig1E
col=rep("dodgerblue",sum(select))
col[data[select,"germline_mutation"]=="POLD1 S478N"]="red"
col[data[select,"germline_mutation"]=="POLD1 L474P"]="limegreen"
col[data[select,"germline_mutation"]=="POLD1 D316N"]="magenta"

pdf("indels_snvs.pdf",useDingbats = F,height=4,width=4)
plot(y=data[select,"indel_burden_corrected"]/data[select,"age"],xlim=c(0,725),
     x=data[select,"sbs_burden_corrected"]/data[select,"age"],ylim=c(0,75),
     pch=21,bg=col,xlab="SBS rate (per year)",ylab="ID rate (per year)",cex=1.4)
points(x=43.6,y=4.4,cex=1.4,lwd=2,pch=4)
dev.off()


select2=grepl("_lo0",data$sample)&data$patient!="PD44594"&!grepl("normal",data$tissue)&
  !data$hist.comment%in%c("lymphocytes","connective tissue")
library(scales)

#Fig1F
pdf("Indels_vs_SNVs_polyps.pdf",useDingbats = F,width=4,height=4)
plot(x=data$sbs_burden_corrected[select2],y=data$indel_burden_corrected[select2],pch=21,bg="white",
     xlim=c(0,max(data$sbs_burden_corrected[select2],na.rm=T)),xlab="SBS burden",cex=1.5,col="white",
     ylim=c(0,max(data$indel_burden_corrected[select2],na.rm=T)),ylab="ID burden")
points(x=data$sbs_burden_corrected[select2],y=data$indel_burden_corrected[select2],cex=1.5,pch=21,bg=col_vec2)
rect(xleft=0,xright=max(data$sbs_burden_corrected[select]),ybottom=0,ytop=max(data$indel_burden_corrected[select]),col=alpha("grey50",0.5))
dev.off()
