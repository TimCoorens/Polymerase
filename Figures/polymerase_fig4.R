#Polymerase paper, Figure 4

#Skin
final_sigs=t(read.table("~/Desktop/Phil_POL/final_sigs_2020_04_21.txt"))
skin_counts=t(data.frame(all=as.numeric(read.table("PD44594f_all_sig_input.txt")[1,]),
                       gene=as.numeric(read.table("PD44594f_all_sig_input_genes.txt")[1,]))
select=c("SBS1","SBS5","SBS10a","SBS10B","SBS28","SBS7a","SBS7b")
colnames(skin_counts)=colnames(final_sigs)=colnames(cosmic_signatures_v2)

skin_counts=as.data.frame(skin_counts)
skin_fit=fit_signatures(counts=skin_counts["all",],signatures = final_sigs[select,],iter = 20000,warmup = 10000,model="poisson",chains = 2)
skin_gene_fit=fit_signatures(counts=skin_counts[2,],opportunities = "human-exome",signatures = final_sigs[select,],iter = 20000,warmup = 10000,model="poisson",chains = 2)

skin_exposures=retrieve_pars(skin_fit,par = "exposures",hpd_prob = 0.95)
skin_gene_exposures=retrieve_pars(skin_gene_fit,par = "exposures",hpd_prob = 0.95)

1/rowSums(skin_exposures$mean[,c("SBS1","SBS5","SBS7a","SBS7b")])
1/rowSums(skin_gene_exposures$mean[,c("SBS1","SBS5","SBS7a","SBS7b")])

1/rowSums(skin_exposures$lower_95[,c("SBS1","SBS5","SBS7a","SBS7b")])
1/rowSums(skin_gene_exposures$lower_95[,c("SBS1","SBS5","SBS7a","SBS7b")])

1/rowSums(skin_exposures$upper_95[,c("SBS1","SBS5","SBS7a","SBS7b")])
1/rowSums(skin_gene_exposures$upper_95[,c("SBS1","SBS5","SBS7a","SBS7b")])

#Endometrium
final_sigs=t(read.table("~/Desktop/Phil_POL/final_sigs_2020_04_21.txt"))
endo_counts=t(data.frame(all=as.numeric(read.table("~/Desktop/Phil_POL/endo_sigs_input_all.txt")[1,]),
                         gene=as.numeric(read.table("~/Desktop/Phil_POL/endo_sigs_input_genes_all.txt")[1,]))
select=c("SBS1","SBS5","SBS10a","SBS10B","SBS28")
colnames(endo_counts)=colnames(final_sigs)=colnames(cosmic_signatures_v2)

endo_counts=as.data.frame(endo_counts)
endo_fit=fit_signatures(counts=endo_counts["all",],signatures = final_sigs[select,],iter = 20000,warmup = 10000,model="poisson",chains = 2)
endo_gene_fit=fit_signatures(counts=endo_counts[2,],opportunities = "human-exome",signatures = final_sigs[select,],iter = 20000,warmup = 10000,model="poisson",chains = 2)

endo_exposures=retrieve_pars(endo_fit,par = "exposures",hpd_prob = 0.95)
endo_gene_exposures=retrieve_pars(endo_gene_fit,par = "exposures",hpd_prob = 0.95)

1/rowSums(endo_exposures$mean[,c("SBS1","SBS5")])
1/rowSums(endo_gene_exposures$mean[,c("SBS1","SBS5")])

1/rowSums(endo_exposures$lower_95[,c("SBS1","SBS5")])
1/rowSums(endo_gene_exposures$lower_95[,c("SBS1","SBS5")])

1/rowSums(endo_exposures$upper_95[,c("SBS1","SBS5")])
1/rowSums(endo_gene_exposures$upper_95[,c("SBS1","SBS5")])


y = as.numeric(endo_counts[3,]); maxy = max(y)
h = barplot(y, las=2, col=hist.cols, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr)
for (j in 1:length(sub_vec)) {
  xpos = h[c((j-1)*16+1,j*16)]
  rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=hist.cols[j*16])
  text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
} 

#Colon

final_sigs=t(read.table("final_sigs_2020_04_21.txt"))
colon_counts=rbind(read.table("colon_sigs_input_all.txt"),
                  read.table("colon_sigs_input_genes.txt"))
hdp_exposures=read.table("mean_assignment_hdp.txt")
colnames(hdp_exposures)=paste0("N",0:13)
hdp2ref=readRDS("hdp2refsigs.Rdata")

patients=rownames(colon_counts)[1:12]

ratios_colon=matrix(0,nrow=length(patients),ncol=6)
colnames(ratios_colon)=outer(c("all","genes"),c("_mean","_lower95","_upper95"), FUN = "paste0")[1:6]
rownames(ratios_colon)=patients
normal=c("SBS1","SBS5")
pol=c("SBS10a","SBS10B","SBS10C","SBS28")

for(patient in patients){
  hdp_sigs_present=colnames(hdp_exposures)[colSums(hdp_exposures[grepl(patient,rownames(hdp_exposures)),]>0.1)>0]
  ref_sigs_present=unique(c(unlist(hdp2ref[hdp_sigs_present]),"SBS1","SBS5"))
  ref_sigs_present[ref_sigs_present=="SBS10b"]="SBS10B"
  if(patient%in%POLE)ref_sigs_present=unique(c(ref_sigs_present,"SBS28"))
  ref_sigs_present=ref_sigs_present[ref_sigs_present!="SBS10D"]
  if(patient=="PD44588")ref_sigs_present=c(ref_sigs_present,"SBSA")
  if(patient=="PD44590")ref_sigs_present=c(ref_sigs_present,"SBS10C")
  
  fit=fit_signatures(counts=colon_counts[patient,],signatures = final_sigs[ref_sigs_present,],iter = 20000,warmup = 10000,model="poisson",chains = 2)
  colon_exposure=retrieve_pars(fit,par = "exposures",hpd_prob = 0.95)
  ratios_colon[patient,"all_mean"]=rowSums(colon_exposure$mean[patient,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$mean[patient,normal])
  ratios_colon[patient,"all_lower95"]=rowSums(colon_exposure$lower_95[patient,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$lower_95[patient,normal])
  ratios_colon[patient,"all_upper95"]=rowSums(colon_exposure$upper_95[patient,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$upper_95[patient,normal])
  
  fit=fit_signatures(counts=colon_counts[paste0(patient,'_gene'),],opportunities = "human-exome",signatures = final_sigs[ref_sigs_present,],iter = 20000,warmup = 10000,model="poisson",chains = 2)
  colon_exposure=retrieve_pars(fit,par = "exposures",hpd_prob = 0.95)
  ratios_colon[patient,"genes_mean"]=rowSums(colon_exposure$mean[1,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$mean[1,normal])
  ratios_colon[patient,"genes_lower95"]=rowSums(colon_exposure$lower_95[1,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$lower_95[1,normal])
  ratios_colon[patient,"genes_upper95"]=rowSums(colon_exposure$upper_95[1,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$upper_95[1,normal])
  
  print(patient)
}
POLD1_S478N=paste0("PD445",c(81,82,85))

ratios_colon_geno=matrix(0,nrow=4,ncol=9)
colnames(ratios_colon_geno)=outer(c("all","genes","cancergenes"),c("_mean","_lower95","_upper95"), FUN = "paste0")[1:9]
rownames(ratios_colon_geno)=c("POLE","POLD1_S478N","POLD1_L474P","POLD1_D316N")
POLE=c("PD44580","PD44586","PD44587","PD44589","PD44591","PD44592","PD44593")
              
all_counts_POLE=rbind(colon_counts[1,],colSums(colon_counts[POLE,]))
fit=fit_signatures(counts=all_counts_POLE[2,],signatures = final_sigs[c("SBS1","SBS5","SBS10a","SBS10B","SBSA","SBS35like","SBS28"),],iter = 20000,warmup = 10000,model="poisson",chains = 2)
colon_exposure=retrieve_pars(fit,par = "exposures",hpd_prob = 0.95)
ratios_colon_geno["POLE","all_mean"]=rowSums(colon_exposure$mean[1,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$mean[1,normal])
ratios_colon_geno["POLE","all_lower95"]=rowSums(colon_exposure$lower_95[1,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$lower_95[1,normal])
ratios_colon_geno["POLE","all_upper95"]=rowSums(colon_exposure$upper_95[1,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$upper_95[1,normal])
all_counts_POLE=rbind(colon_counts[1,],colSums(colon_counts[paste0(POLE,"_gene"),]))
fit=fit_signatures(counts=all_counts_POLE[2,],signatures = final_sigs[c("SBS1","SBS5","SBS10a","SBS10B","SBSA","SBS35like","SBS28"),],iter = 20000,warmup = 10000,model="poisson",chains = 2)
colon_exposure=retrieve_pars(fit,par = "exposures",hpd_prob = 0.95)
ratios_colon_geno["POLE","genes_mean"]=rowSums(colon_exposure$mean[1,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$mean[1,normal])
ratios_colon_geno["POLE","genes_lower95"]=rowSums(colon_exposure$lower_95[1,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$lower_95[1,normal])
ratios_colon_geno["POLE","genes_upper95"]=rowSums(colon_exposure$upper_95[1,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$upper_95[1,normal])

POLD1_S478N=paste0("PD445",c(81,82,85))
all_counts_POLD1=rbind(colon_counts[1,],colSums(colon_counts[POLD1_S478N,]))
fit=fit_signatures(counts=all_counts_POLD1[2,],signatures = final_sigs[c("SBS1","SBS5","SBS10C","SBSA","SBS28"),],iter = 20000,warmup = 10000,model="poisson",chains = 2)
colon_exposure=retrieve_pars(fit,par = "exposures",hpd_prob = 0.95)
ratios_colon_geno["POLD1_S478N","all_mean"]=rowSums(colon_exposure$mean[1,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$mean[1,normal])
ratios_colon_geno["POLD1_S478N","all_lower95"]=rowSums(colon_exposure$lower_95[1,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$lower_95[1,normal])
ratios_colon_geno["POLD1_S478N","all_upper95"]=rowSums(colon_exposure$upper_95[1,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$upper_95[1,normal])
all_counts_POLD1=rbind(colon_counts[1,],colSums(colon_counts[paste0(POLD1_S478N,"_gene"),]))
fit=fit_signatures(counts=all_counts_POLD1[2,],signatures = final_sigs[c("SBS1","SBS5","SBS10C","SBSA","SBS28"),],opportunities = "human-exome",iter = 20000,warmup = 10000,model="poisson",chains = 2)
colon_exposure=retrieve_pars(fit,par = "exposures",hpd_prob = 0.95)
ratios_colon_geno["POLD1_S478N","genes_mean"]=rowSums(colon_exposure$mean[1,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$mean[1,normal])
ratios_colon_geno["POLD1_S478N","genes_lower95"]=rowSums(colon_exposure$lower_95[1,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$lower_95[1,normal])
ratios_colon_geno["POLD1_S478N","genes_upper95"]=rowSums(colon_exposure$upper_95[1,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$upper_95[1,normal])

patient="PD44588"
ref_sigs_present=c("SBS1","SBS5","SBS10C","SBSA","SBSB","SBS28")
fit=fit_signatures(counts=colon_counts[patient,],signatures = final_sigs[ref_sigs_present,],iter = 20000,warmup = 10000,model="poisson",chains = 2)
colon_exposure=retrieve_pars(fit,par = "exposures",hpd_prob = 0.95)
ratios_colon_geno["POLD1_L474P","all_mean"]=rowSums(colon_exposure$mean[patient,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$mean[patient,normal])
ratios_colon_geno["POLD1_L474P","all_lower95"]=rowSums(colon_exposure$lower_95[patient,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$lower_95[patient,normal])
ratios_colon_geno["POLD1_L474P","all_upper95"]=rowSums(colon_exposure$upper_95[patient,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$upper_95[patient,normal])
fit=fit_signatures(counts=colon_counts[paste0(patient,'_gene'),],opportunities = "human-exome",signatures = final_sigs[ref_sigs_present,],iter = 20000,warmup = 10000,model="poisson",chains = 2)
colon_exposure=retrieve_pars(fit,par = "exposures",hpd_prob = 0.95)
ratios_colon_geno["POLD1_L474P","genes_mean"]=rowSums(colon_exposure$mean[1,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$mean[1,normal])
ratios_colon_geno["POLD1_L474P","genes_lower95"]=rowSums(colon_exposure$lower_95[1,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$lower_95[1,normal])
ratios_colon_geno["POLD1_L474P","genes_upper95"]=rowSums(colon_exposure$upper_95[1,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$upper_95[1,normal])

patient="PD44590"
ref_sigs_present=c("SBS1","SBS5","SBS10C","SBSA","SBS28")
fit=fit_signatures(counts=colon_counts[patient,],signatures = final_sigs[ref_sigs_present,],iter = 20000,warmup = 10000,model="poisson",chains = 2)
colon_exposure=retrieve_pars(fit,par = "exposures",hpd_prob = 0.95)
ratios_colon_geno["POLD1_D316N","all_mean"]=rowSums(colon_exposure$mean[patient,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$mean[patient,normal])
ratios_colon_geno["POLD1_D316N","all_lower95"]=rowSums(colon_exposure$lower_95[patient,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$lower_95[patient,normal])
ratios_colon_geno["POLD1_D316N","all_upper95"]=rowSums(colon_exposure$upper_95[patient,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$upper_95[patient,normal])
fit=fit_signatures(counts=colon_counts[paste0(patient,'_gene'),],opportunities = "human-exome",signatures = final_sigs[ref_sigs_present,],iter = 20000,warmup = 10000,model="poisson",chains = 2)
colon_exposure=retrieve_pars(fit,par = "exposures",hpd_prob = 0.95)
ratios_colon_geno["POLD1_D316N","genes_mean"]=rowSums(colon_exposure$mean[1,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$mean[1,normal])
ratios_colon_geno["POLD1_D316N","genes_lower95"]=rowSums(colon_exposure$lower_95[1,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$lower_95[1,normal])
ratios_colon_geno["POLD1_D316N","genes_upper95"]=rowSums(colon_exposure$upper_95[1,ref_sigs_present%in%c(normal,pol)])/rowSums(colon_exposure$upper_95[1,normal])

POLE_nanoseq_exposures_mean=read.table("POLE_nanoseq_exposures_mean.txt",sep="\t")
POLD_nanoseq_exposures_mean=read.table("POLD_nanoseq_exposures_mean.txt",sep="\t")

wg_means=c(ratios_colon_geno["POLE","all_mean"], 
           1/rowSums(endo_exposures$mean[1,c("SBS1","SBS5")]),
           1/rowSums(POLE_nanoseq_exposures_mean$mean[2:1,c("SBS1","SBS5")]),
           1/rowSums(skin_exposures$mean[1,c("SBS1","SBS5","SBS7a","SBS7b")]),
           1/rowSums(POLD_nanoseq_exposures_mean[2,c("SBS1","SBS5")]),
           ratios_colon_geno["POLD1_S478N","all_mean"],
           1/rowSums(POLD_nanoseq_exposures_mean[1,c("SBS1","SBS5")]),
           ratios_colon_geno[c("POLD1_L474P","POLD1_D316N"),"all_mean"])

#Fig4A
pdf("wg_increase.pdf",useDingbats = F)
plot(wg_means,ylab="Mutation rate relative to SBS1&5",cex=2,pch=16,x=c(1:5,7:9,11,13), bty='l',ylim=c(1,12))
abline(v=c(6,10,12),lty='dashed')
abline(h=1,lty='dashed')
dev.off()

#Fig4B
coding_means=c(ratios_colon_geno["POLE","genes_mean"], 
           1/rowSums(endo_gene_exposures$mean[1,c("SBS1","SBS5")]),
           1/rowSums(skin_gene_exposures$mean[1,c("SBS1","SBS5","SBS7a","SBS7b")]),
           ratios_colon_geno["POLD1_S478N","genes_mean"],
           ratios_colon_geno[c("POLD1_L474P","POLD1_D316N"),"genes_mean"])

1/rowSums(skin_gene_exposures$lower_95[1,c("SBS1","SBS5","SBS7a","SBS7b")])
1/rowSums(skin_gene_exposures$upper_95[1,c("SBS1","SBS5","SBS7a","SBS7b")])

pdf("coding_increase.pdf",useDingbats = F,height=6,width=5)
plot(coding_means,ylab="Mutation rate relative to SBS1&5",cex=2,pch=16,x=c(1:3,5,7,9), bty='l',ylim=c(1,2),xlim=c(0,10))
abline(v=c(4,6,8,10),lty='dashed')
dev.off()

