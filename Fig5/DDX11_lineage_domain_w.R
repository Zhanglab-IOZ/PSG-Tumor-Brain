libray(pheatmap)
library(RColorBrewer)
#sum1.tab comtain the W ratio of different lineage in differnt domain
sum1<-read.table("sum1.tab",header=T,row.names=1)
pheatmap(sum1,scale="none",color=colorRampPalette(colors=c("white","firebrick1"))(100),cluster_rows=FALSE,cluster_cols=FALSE,display_numbers=T,number_color="black",filename="domainheat2.pdf",height=3,width=7)