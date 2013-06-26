#!/bin/Rscript
args <- commandArgs(TRUE)
if (length(args)<3) {
stop("plot_colocalization.R <input_file> <cluster_n> <output.pdf>\nPlease specify all three parameters")}
input = args[1]
cluster_n = eval (parse(text=args[2]))
output = args[3]

library("fields")
library("reshape")
library("scales")
library("ggplot2")
library('grid')

data = read.table(input,sep='\t',header=T,nrow=50000)
if ((ncol(data)-8)%%(cluster_n+1) != 0){
stop("cluster_n is inconsistent with input data")}
hist_n = (ncol(data)-8)/(cluster_n+1)
hist_name = colnames(data)[(9+hist_n*cluster_n):ncol(data)]
print(hist_name)
co.matrix <- function(cluster){
  # find the matrix for co-loc log-ratios amoung hist-marks in given subpopulation
  tmp_data = data[(9+hist_n*(cluster-1)):(8+hist_n*cluster)]
  co_matrix = apply(tmp_data,2,
                  function(y) 
                       apply(tmp_data,2,function(x) sum(y*x)/sum(x))/sum(y)
                 )*nrow(data)
  for (i in 1:hist_n) {co_matrix[i,i] = 1}
  co_matrix = log(co_matrix)
  colnames(co_matrix)=rownames(co_matrix)=hist_name
  min = min(co_matrix[co_matrix!=-Inf])
  max = max(co_matrix[co_matrix!=-Inf])
  print (c(min,max))
  b=c(seq(min,0,length.out=50),seq(0,max,length.out=50))
  my.colors = colorRampPalette(c("blue",'white','red'))(100)
  po.nopanel <- list(opts(panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank()))
  df_matrix=melt(co_matrix)
  df_matrix$X1 <- factor(df_matrix$X1, levels=rev(unique(as.character(df_matrix$X2)))) # reorder one dimension
  plot = ggplot(melt(co_matrix))+geom_tile(aes(x=X1,y=X2,fill=value),colour='black')+scale_fill_gradientn(values=rescale(b),colours=my.colors,name='log-ratio')+po.nopanel
  plot = plot + theme(axis.text.x = element_text(size=5,colour='black'), 
                      axis.text.y = element_text(size=5,colour='black'),
                      plot.title = element_text(size=10))
  plot = plot + xlab('') + ylab('') + ggtitle(paste("Co-localization in sub-population",cluster,sep=' '))
  #mtext(text=1:7,side=2,line=0.3,at=seq(0,(hist_n-1))/(hist_n-1),las=1)
  #mtext(text=hist_name,side=1,line=0.3,at=seq(0,(hist_n-1))/(hist_n-1),las=1)
  return(list(plot=plot,co_matrix=co_matrix))
}


pdf(output,width=5,height=6)
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 1)))
for (i in 1:cluster_n) {
result=co.matrix(i)
print(result$plot, vp = vplayout(i,1))
}
dev.off()
