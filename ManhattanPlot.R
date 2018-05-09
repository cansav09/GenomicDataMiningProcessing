# Manhattan plots
# Author: Candace Savonen 
# Last Update: 11/10/17

library(qqman)
library(beepr)
home="/Users/cansav091/Desktop/ManhattanPlot"
setwd(home)
#home="/volumes/transsci/Bernstein_Lab/Projects/P1/E5/Differential_Methylation_Analysis/Beta_regression/"
fdr.files=grep("FDR.txt",dir(home),value=TRUE)
# Load the already made concatenated files with Beta regression p values
load("/volumes/transsci/Bernstein_Lab/Projects/P1/E5/Differential_Methylation_Analysis/Beta_regression/fivehmc_glmmtmb_out.RData")
load("/volumes/transsci/Bernstein_Lab/Projects/P1/E5/Differential_Methylation_Analysis/Beta_regression/fivemc_glmmtmb_out.RData")

# Load the annotation for these so I can retreieve the chr locations of the cpg's.
load("/volumes/transsci/Bernstein_Lab/Projects/P1/E5/5mC_5hmc_Estimation/Manifest/EPIC.manifest.rda")

# The list of datasets we need plots from. Need a 2 plots from each of these datasets. Used a loop so its all consistent. No copy pasting.
datasets=c("fivehmc.glmmtmb.all.cing","fivehmc.glmmtmb.all.pari","fivemc.glmmtmb.cingulate","fivemc.glmmtmb.parietal")

mod=c("5hmc","5hmc","5mc","5mc") # For labeling purposes
tissues=rep(c("Cingulate","Parietal"),2) # For labeling purposes

if(dir.exists(paste0(home,"/ManhattanPlots"))==FALSE){ # For storing the graphs
  dir.create(paste0(home,"/ManhattanPlots"))
}
tissues.f=c("Cingulate","Cingulate","Parietal","Parietal")
stage=c("limbic","neocortical","limbic","neocortical")

fdr.p=data.frame(t(c(rep(0,4))))
colnames(fdr.p)=c("FC","HolP","tissues","stage")
fdr.p=fdr.p[-1,]

for(ii in 1:length(fdr.files)){
xx=read.table(paste0(home,"/",fdr.files[ii]),row.names=1,skip=1)
xx=cbind(xx,rep(tissues.f[ii],nrow(xx)),rep(stage[ii],nrow(xx)))
colnames(xx)=c("FC","HolP","tissues","stage")
fdr.p=rbind(fdr.p,xx)
}


for(ii in 1:length(datasets)){
data=eval(parse(text=datasets[ii]))# Set the next dataset in our list for this iteration of the loop.

probes=rownames(data) 

# Match the probe names to the manifest annotation dataset. 
if(ii<3){
yy=fdr.p[which(fdr.p[,3]==tissues[ii]),]
yy=yy[match(probes,rownames(yy)),]
yy=yy[which(!is.na(yy[,1])),]

yy.l=yy[which(yy[,4]=="limbic"),]
yy.n=yy[which(yy[,4]=="neocortical"),]

  probes.l=rownames(yy.l)
}else{
  probes.l=rownames(data) 
}

xx.l=match(probes.l,EPIC.manifest@ranges@NAMES)

chrs.l=EPIC.manifest@elementMetadata@listData$chrmA[xx.l] 
#manhattan function requires chromosomes to be noted as numeric vector with X, Y, and MT chrs being 23:25 respectively. 
chrs.l=gsub("chr","",chrs.l) 
chrs.l=gsub("X",23,chrs.l)
chrs.l=gsub("Y",24,chrs.l)
chrs.l=gsub("MT",25,chrs.l)

if(ii<3){
probes.n=rownames(yy.n)
xx.n=match(probes.n,EPIC.manifest@ranges@NAMES)

chrs.n=EPIC.manifest@elementMetadata@listData$chrmA[xx.n] 
#manhattan function requires chromosomes to be noted as numeric vector with X, Y, and MT chrs being 23:25 respectively. 
chrs.n=gsub("chr","",chrs.n) 
chrs.n=gsub("X",23,chrs.n)
chrs.n=gsub("Y",24,chrs.n)
chrs.n=gsub("MT",25,chrs.n)

# Make a dataframe with the CPG's, chromsomes, bp start position, and FDR q values. 
manh.data.l=data.frame(probes.l,         
           as.numeric(chrs.l),
           EPIC.manifest@ranges@start[xx.l],
           -log10(yy.l[,2]),
           stringsAsFactors = FALSE)

manh.data.n=data.frame(probes.n,         
           as.numeric(chrs.n),
           EPIC.manifest@ranges@start[xx.n],
           -log10(yy.n[,2]),
           stringsAsFactors = FALSE)
colnames(manh.data.l)=colnames(manh.data.n)=c("CPG","CHR","BP","P")

}else{
  manh.data.l=data.frame(probes.l,         
                         as.numeric(chrs.l),
                         EPIC.manifest@ranges@start[xx.l],
                         -log10(data$LimbicVSNoneFDR),
                         -log10(data$NeocorticalVSNoneFDR),
                         stringsAsFactors = FALSE)
  colnames(manh.data.l)=c("CPG","CHR","BP","P","NP")
  manh.data.l=manh.data.l[-which(manh.data.l$NP==Inf),]
  
  xx=unique(which(is.nan(manh.data.l$NP)),which(is.nan(manh.data.l$P)))
  if(length(xx)>0){
  manh.data.l=manh.data.l[-xx,]
  }
}

# Label them as such. 

manh.data.l=manh.data.l[-which(manh.data.l$P==Inf),]
manh.data.l$CPG=as.character(manh.data.l$CPG) # CPG's need to be a character vector
manh.data.l$CHR=as.numeric(manh.data.l$CHR) # Chromsomal locations need to be numeric
manh.data.l$BP= as.numeric(manh.data.l$BP)
manh.data.l=manh.data.l[which(!is.na(manh.data.l$CHR)),]

if(ii<3){
manh.data.n=manh.data.n[-which(manh.data.n$P==Inf),]
manh.data.n$CPG=as.character(manh.data.n$CPG) # CPG's need to be a character vector
manh.data.n$CHR=as.numeric(manh.data.n$CHR) # Chromsomal locations need to be numeric
manh.data.n$BP= as.numeric(manh.data.n$BP)
manh.data.n=manh.data.n[which(!is.na(manh.data.n$CHR)),]

}
# Adjust the 0's to be the minimum value that R seems to retain. 
# This way we don't get Infinity as our -log10(0). We can't plot Inf.
#manh.data$LP[which(manh.data$LP==0)]= min(manh.data$LP[-which(manh.data$LP==0)])
#manh.data$NP[which(manh.data$NP==0)]= min(manh.data$NP[-which(manh.data$NP==0)])

# Manhattan Plot for Limbic data
jpeg(paste0(home,"/ManhattanPlots/",datasets[ii],"MidManhattan2.jpeg"))
if(ii<3){
sig=manh.data.l$CPG[which(manh.data.l$P>3)]
manhattan(manh.data.l,chr="CHR",bp="BP",p="P",snp="CPG",logp=F,ylim=c(0,round(max(manh.data.l$P))+10),highlight=sig,chrlabs=c(1:22, "X", "Y"),suggestiveline=FALSE,genomewideline=TRUE,main=paste( mod[ii],tissues[ii],"Mid Stage Disease"))
}else{
sig=manh.data.l$CPG[which(manh.data.l$P>3)]
manhattan(manh.data.l,chr="CHR",bp="BP",p="P",snp="CPG",ylim=c(0,round(max(manh.data.l$P[!is.na(manh.data.l$P)]))),logp=F,highlight=sig,chrlabs=c(1:22, "X", "Y"),suggestiveline=FALSE,genomewideline=TRUE,main=paste( mod[ii],tissues[ii],"Mid Stage Disease"))
}
# this vector will allow us to highlight probes that exceed our cutoff. 
abline(h=-log10(.001),col="red") # plot cutoff line. 
dev.off()

# Manhattan Plot for Neocortical data
jpeg(paste0(home,"/ManhattanPlots/",datasets[ii],"LateManhattan2.jpeg"))
if(ii<3){
  sig=manh.data.n$CPG[which(manh.data.n$P>3)]
  manhattan(manh.data.n,chr="CHR",bp="BP",p="P",snp="CPG",ylim=c(0,round(max(manh.data.l$P))+10),logp=F,highlight=sig,chrlabs=c(1:22, "X", "Y"),suggestiveline=FALSE,genomewideline=TRUE,main=paste(mod[ii],tissues[ii],"Late Stage Disease"))
}else{
  sig=manh.data.l$CPG[which(manh.data.l$NP>3)]
  manhattan(manh.data.l,chr="CHR",bp="BP",p="NP",snp="CPG",ylim=c(0,round(max(manh.data.l$NP))),logp=F,highlight=sig,chrlabs=c(1:22, "X", "Y"),suggestiveline=FALSE,genomewideline=TRUE,main=paste(mod[ii],tissues[ii],"Late Stage Disease"))
}
  abline(h=-log10(.001),col="red")
  dev.off()


}
# Just a nifty way to signal that your graphs are finished being made. 
beep(sound=2)


