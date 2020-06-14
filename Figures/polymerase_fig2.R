#Fit signatures to observed counts
options(stringsAsFactors = F)
library(sigfit)
library(RColorBrewer)
library(ape)
library(ggtree)
data("cosmic_signatures_v2")
final_sigs=t(read.table("final_sigs_2020_04_21.txt"))
all_counts=read.table("trinuc_mut_mat_2020_03_17.txt")
hdp_exposures=read.table("mean_assignment_hdp.txt")
colnames(hdp_exposures)=paste0("N",0:13)
hdp2ref=readRDS("hdp2refsigs.Rdata")
colnames(all_counts)=colnames(final_sigs)=colnames(cosmic_signatures_v2)
POLE=c("PD44580","PD44586","PD44587","PD44588","PD44589","PD44592","PD44593","PD44594")

all_cols=rep(NA,nrow(final_sigs))
names(all_cols)=rownames(final_sigs)
all_cols[c("SBS10a","SBS10B","SBS28","SBS10E")]=c("lightblue1","deepskyblue","dodgerblue3","navy")
all_cols[c("SBS10C","SBS10D")]=brewer.pal("Greens",n = 3)[2:3]
all_cols[c("SBS1","SBS5")]=c("salmon","firebrick")
all_cols[c("SBSA","SBSB")]=c("grey60","grey20")
all_cols[c("SBS7a","SBS7b")]=c("plum2","mediumorchid3")

all_cols[c("SBS35like","SBS17b")]=c("orange","yellow2")

sig_order=1:16
names(sig_order)=c("SBS10a","SBS10B","SBS10C","SBS10D","SBS28","SBS10E","SBS1","SBS5","SBSA","SBSB","SBS7a","SBS7b","SBS35like","SBS17b")

patients=unique(substr(rownames(all_counts),1,7))
for(patient in patients[patients!="PD44584"]){
hdp_sigs_present=colnames(hdp_exposures)[colSums(hdp_exposures[grepl(patient,rownames(hdp_exposures)),]>0.1)>0]
ref_sigs_present=unique(c(unlist(hdp2ref[hdp_sigs_present]),"SBS1","SBS5"))
ref_sigs_present[ref_sigs_present=="SBS10b"]="SBS10B"
if(patient%in%POLE)ref_sigs_present=unique(c(ref_sigs_present,"SBS28"))
 ref_sigs_present=ref_sigs_present[ref_sigs_present!="SBS10D"]
if(patient=="PD44588")ref_sigs_present=c(ref_sigs_present,"SBSA")
 if(patient=="PD44590")ref_sigs_present=c(ref_sigs_present,"SBS10C")
 
counts_patient=all_counts[grepl(patient,rownames(all_counts)),]
colnames(counts_patient)=colnames(cosmic_signatures_v2)
ref_sigs_present=names(sort(sig_order[ref_sigs_present]))

fit=fit_signatures(counts=counts_patient,signatures = final_sigs[ref_sigs_present,],
                   iter = 20000,warmup = 10000,model="poisson",chains = 2)
pars <- retrieve_pars(fit,par = "exposures",hpd_prob = 0.95)
saveRDS(pars,paste0(patient,"_restricted.rda"))
library(ape)
library(ggtree)
tree = read.tree(paste0(patient,".snp_tree_with_branch_length.tree"))

tree_df=fortify(tree)

cols=all_cols[ref_sigs_present]

exposures=t(pars$mean[,ref_sigs_present])
samples=colnames(exposures)[grepl(patient,colnames(exposures))]
branches=substr(samples,9,nchar(samples))

pdf(paste0(patient,"_tree_with_ref_signatures_restricted.pdf"))
plot(tree,label.offset=0.01*max(tree_df$x),show.tip.label=F)
for (k in 1:length(samples)){
  n=as.numeric(branches[k])
  x_end=tree_df$x[n]
  x_start=tree_df$x[tree_df$parent[n]]
  x_intv=x_end-x_start
  y=node.height(tree)[n]
  tipnum=sum(tree_df$isTip)
  for (s in ref_sigs_present){
    x_end=x_start+exposures[s,samples[k]]*x_intv
    rect(ybottom=y-min(0.02*tipnum,0.4),ytop=y+min(0.02*tipnum,0.4),xleft=x_start,xright=x_end,col=cols[s],lwd=0.5)
    #rect(ybottom=y-0.4,ytop=y+0.4,xleft=x_start,xright=x_end,col=cols[s],lwd=0.5)
    
    x_start=x_end
  }
}
axisPhylo(side = 1,backward=F)
legend("topright",title="Signatures", legend=ref_sigs_present, 
      fill=cols, bty="n",cex=0.8, ncol=1, xjust=0.5)
dev.off()

tree_collapsed=tree
tree_collapsed$edge.length=rep(1,nrow(tree_collapsed$edge))
tree_collapsed$edge.length[tree$edge.length==0]=0
tree_collapsed_df=fortify(tree_collapsed)

pdf(paste0(patient,"_tree_collapsed_with_ref_signatures.pdf"))
plot(tree_collapsed,label.offset=0.01*max(tree_collapsed_df$x))
for (k in 1:length(samples)){
  n=as.numeric(branches[k])
  x_end=tree_collapsed_df$x[n]
  x_start=tree_collapsed_df$x[tree_collapsed_df$parent[n]]
  x_intv=x_end-x_start
  y=node.height(tree_collapsed)[n]
  tipnum=sum(tree_collapsed_df$isTip)

  for (s in ref_sigs_present){
    x_end=x_start+exposures[s,samples[k]]*x_intv
    rect(ybottom=y-min(0.02*tipnum,0.4),ytop=y+min(0.02*tipnum,0.4),xleft=x_start,xright=x_end,col=cols[s],lwd=0.5)
    x_start=x_end
  }
}
axisPhylo(side = 1,backward=F)
legend("topright",title="Signatures", legend=ref_sigs_present,
       fill=cols, bty="n",cex=0.8, ncol=1, xjust=0.5)
dev.off()
print(patient)
}

#------------------------------
#PD44584
patient="PD44584"
hdp_sigs_present=colnames(hdp_exposures)[colSums(hdp_exposures[grepl(patient,rownames(hdp_exposures)),]>0.1)>0]
ref_sigs_present=unique(c(unlist(hdp2ref[hdp_sigs_present]),"SBS1","SBS5"))
ref_sigs_present=ref_sigs_present[ref_sigs_present!="SBS10D"]

tree = read.tree(paste0(patient,".snp_tree_with_branch_length.tree"))
tree_df=fortify(tree)

find_parent=function(sample,tree_df){
  child=which(tree_df$label==sample)
  parent=tree_df$parent[child]
  all_nodes=c(child,parent)
  while(child!=parent){
    child=parent
    parent=tree_df$parent[child]
    all_nodes=c(all_nodes,parent)
  }
  return(unique(all_nodes))
}
polyp_samples=paste0("PD44584c_lo000",1:3)
polyp_branches=unique(unlist(lapply(polyp_samples,find_parent,tree_df)))
polyp=rownames(counts_patient)%in%paste0("PD44584_",polyp_branches)

counts_patient=all_counts[grepl(patient,rownames(all_counts)),]
ref_sigs_present=names(sort(sig_order[ref_sigs_present]))

fit=fit_signatures(counts=counts_patient[!polyp,],signatures = final_sigs[ref_sigs_present,],
                     iter = 20000,warmup = 10000,model="poisson",chains = 2)
fit_polyp=fit_signatures(counts=counts_patient[polyp,],signatures = final_sigs[c(ref_sigs_present,"SBS10D"),],
                         iter = 20000,warmup = 10000,model="poisson",chains = 2)
pars_normal <- retrieve_pars(fit,par = "exposures",hpd_prob = 0.95)
pars_polyp <- retrieve_pars(fit_polyp,par = "exposures",hpd_prob = 0.95)
pars=list(mean=rbind(pars_polyp$mean,cbind(pars_normal$mean,data.frame(SBS10D=0))),
          lower_95=rbind(pars_polyp$lower_95,cbind(pars_normal$lower_95,data.frame(SBS10D=0))),
          upper_95=rbind(pars_polyp$upper_95,cbind(pars_normal$upper_95,data.frame(SBS10D=0))))
saveRDS(pars,paste0(patient,".rda"))

ref_sigs_present=c("SBS10C","SBS10D","SBS1","SBS5","SBSA")
cols=all_cols[ref_sigs_present]
  
exposures=t(pars$mean[,ref_sigs_present])
samples=colnames(exposures)[grepl(patient,colnames(exposures))]
branches=substr(samples,9,nchar(samples))
  
pdf(paste0(patient,"_tree_with_ref_signatures_figure.pdf"),width=25)
plot(tree,label.offset=0.01*max(tree_df$x))
for (k in 1:length(samples)){
  n=as.numeric(branches[k])
  x_end=tree_df$x[n]
  x_start=tree_df$x[tree_df$parent[n]]
  x_intv=x_end-x_start
  y=node.height(tree)[n]
  tipnum=sum(tree_df$isTip)
  for (s in ref_sigs_present){
    x_end=x_start+exposures[s,samples[k]]*x_intv
    rect(ybottom=y-min(0.02*tipnum,0.4),ytop=y+min(0.02*tipnum,0.4),xleft=x_start,xright=x_end,col=cols[s],lwd=0.5)
    x_start=x_end
  }
}
axisPhylo(side = 1,backward=F)
legend("right",title="Signatures", legend=ref_sigs_present, 
       fill=cols, bty="n",cex=0.8, ncol=1, xjust=0.5)
dev.off()
  
  tree_collapsed=tree
  tree_collapsed$edge.length=rep(1,nrow(tree_collapsed$edge))
  tree_collapsed$edge.length[tree$edge.length==0]=0
  tree_collapsed_df=fortify(tree_collapsed)
  
  pdf(paste0(patient,"_tree_collapsed_with_ref_signatures.pdf"))
  plot(tree_collapsed,label.offset=0.01*max(tree_collapsed_df$x))
  for (k in 1:length(samples)){
    n=as.numeric(branches[k])
    x_end=tree_collapsed_df$x[n]
    x_start=tree_collapsed_df$x[tree_collapsed_df$parent[n]]
    x_intv=x_end-x_start
    y=node.height(tree_collapsed)[n]
    tipnum=sum(tree_collapsed_df$isTip)
    
    for (s in ref_sigs_present){
      x_end=x_start+exposures[s,samples[k]]*x_intv
      rect(ybottom=y-min(0.02*tipnum,0.4),ytop=y+min(0.02*tipnum,0.4),xleft=x_start,xright=x_end,col=cols[s],lwd=0.5)
      x_start=x_end
    }
  }
  axisPhylo(side = 1,backward=F)
  legend("topright",title="Signatures", legend=ref_sigs_present, 
         fill=cols, bty="n",cex=0.8, ncol=1, xjust=0.5)
  dev.off()
  
