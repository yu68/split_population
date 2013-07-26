#!/usr/bin/env Rscript
#generate relative intensities of histone modification in single nucleosome from absolute intensities.
# Version: 0.1
setwd("./")

list.of.packages <- c("argparse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos="http://R-Forge.R-project.org")
suppressPackageStartupMessages(require(argparse))

parser <- ArgumentParser(description="Generate relative intensities of histone modification in single nucleosome from absolute intensities",epilog="Need 'argparse' library")
parser$add_argument("input",help="input files contain absolute intensities of histone modification in mononucleosome level")
parser$add_argument("-p","--pdf",default="Distribution_of_Max_Intensity.pdf",help="pdf file to store the distribution of maximum intensities for each histone modification. [default \"%(default)s\"]")
parser$add_argument("-o","--output",default="Epi_Intensity_nucleosome_relative.txt",help="output txt file storing the relative intensities. [default \"%(default)s\"]")
parser$add_argument("-s","--size",type="integer",default=5000,help="size of the sampled sublist for distribution of max values. [default: %(default)s]")
parser$add_argument("-i","--iteration",type="integer",default=1000,help="iteration number to get a pool of max values for its distribution. [default: %(default)s]")
parser$add_argument("--pvalue",type="double",default=0.01,help="pvalue cutoff of fuzziness and summit value to select high-quality nucleosomes. [default: %(default)s]")

args <- parser$parse_args()

input=args$input
pdf=args$pdf   
output=args$output   
size=args$size
iter=args$iteration
p=args$pvalue

Epi_Intensity=read.table(input,sep='\t',header=T)
n_hmark=(ncol(Epi_Intensity)-8)/2 ## number of histone modifications
Epi_Intensity=Epi_Intensity[rowSums(Epi_Intensity[,9:(8+n_hmark)])!=0,]
Epi_Intensity=Epi_Intensity[(Epi_Intensity$fuzziness_pval<p)&(Epi_Intensity$smt_pval<p),]



find_max<-function(mark_data,size,iter,name)
{
  max_pool=NULL
  for (i in 1:iter){
    max_i=max(mark_data[sample(length(mark_data))<=size])
    max_pool=c(max_pool,max_i)
  }
  box=boxplot(max_pool,range=3,plot=F)
  layout(matrix(c(1,2),ncol=1),heights=c(0.4,0.1))
  min=as.integer(box$stats[1,1])-5
  max=as.integer(box$stats[5,1])+10
  par(mar=c(0,4,3,1))
  a=hist(max_pool,breaks=c(0:max,1E8),xlim=c(min,max),xaxt='n',main=name,col='grey')
  
  cutoff=box$stats[5,1]
  num_outlier=sum(max_pool>=cutoff)
  legend("topright",paste(num_outlier,"/",iter," points outliers",sep=''))
  par(mar=c(4,4,0,1))
  boxplot(max_pool,range=3,ylim=c(min,max),horizontal=T,frame=F,col='grey')
  mtext("Max intensity",1,line=2)
  
  mean_max=mean(max_pool[max_pool<cutoff])
  sd_max=sd(max_pool[max_pool<cutoff])
  relative_data=mark_data/mean_max
  cat("mean and S.D. of sampled max values\t",c(mean_max,sd_max),'\n')
  return(relative_data)
}


#find relative intensities
Epi_Intensity_relative=Epi_Intensity
selection=Epi_Intensity$H3K4me3>=0
pdf(pdf)
cat("Distribution of sampled max values:\n")
for (i in 9:(8+n_hmark)) {
  cat("  ",colnames(Epi_Intensity)[i],'\t')
  Epi_Intensity_relative[,i]=find_max(Epi_Intensity[,i],size,iter,colnames(Epi_Intensity)[i])
  Epi_Intensity_relative[Epi_Intensity_relative[,i]>1,i]=1
}
dev.off()
#choose only ones with relative intensity less than one
write.table(Epi_Intensity_relative[,c(1:(8+n_hmark))],output,sep='\t',col.names=T,row.names=F,quote=F)


