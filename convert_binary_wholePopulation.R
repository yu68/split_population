setwd("~/split_population/program")
#!/bin/Rscript
args <- commandArgs(TRUE)
if (length(args)<4) {
  stop("Rscript convert_binary_wholePopulation.R <input_raw_intensity1> <input_raw_intensity2> <pvalue> <output>")
}
input1=args[1]
input2=args[2]
p=eval(parse(text=args[3]))
output=args[4]
Epi_Intensity=read.table(input1,sep='\t',header=T)
a=read.table(input2,sep='\t',header=T)
Epi_Intensity=cbind(Epi_Intensity,a[,9:12])
rm(a)
Epi_Intensity[,9:16]=apply(Epi_Intensity[,9:16],2,function(x) (x>=qpois((1-p),mean(x)))*1)
write.table(Epi_Intensity,output,sep='\t',col.names=T,row.names=F,quote=F)


names=gsub("H3","",colnames(Epi_Intensity)[c(9:13,15:16)])

barplot_count <- function(data,title){
counts=NULL
for (i in 0:7) {
  index=rowSums(data)==i
  print (paste("# of Nucleosome with ",i," marks: ",sum(index),sep=''))
  if ((i==0)|(i>=4)) { next; }
  comb=apply(data[index,],1,function(x) 
                            do.call(paste, c(as.list(names[x==1]), sep="/")) )
  count=as.data.frame(table(comb))
  count=count[order(-count$Freq),]
  if(i!=1) { count=count[count$Freq>0.04*sum(index),] }
  counts=rbind(counts,count)
  print (count)
  
}
colors=c('#DDECEF','#BFD9DA','#87AAAE')
Num = nchar(gsub("[^/]","",counts$comb))+1 # histone modification numbers to count in counts$comb variable

par(mar=c(8,4,3,0),cex.axis=0.8)
barplot(counts$Freq,names.arg=counts$comb,las=2,col=colors[Num],main=title,
        legend=1:3, args.legend = list(fill=colors,title='# of marks:',ncol=3,x='topright'))
}

pdf("epi-patterns_in_single_nucleosome.pdf",width=6,height=3.5)
barplot_count(data = Epi_Intensity[,c(9:13,15:16)],"Whole population")
dev.off()

#for sub-populations
tmp=read.table("Clut-3_Hist-7_Nucleo-974874_beta-500_matrix.txt",sep='\t',header=T)
pdf("epi-patterns_in_single_nucleosome_subP.pdf",width=6,height=10)
par(mfrow=c(3,1))
barplot_count(data=tmp[,9:15],"Sub-population 1")
barplot_count(data=tmp[,16:22],"Sub-population 2")
barplot_count(data=tmp[,23:29],"Sub-population 3")


dev.off()
