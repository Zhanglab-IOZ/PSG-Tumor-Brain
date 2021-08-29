library(dplyr)
library(ggplot2)
library(reshape2)

load("integrate_msigdb_all_result.rds")
load("integrate_signal_gene_result.rds")

PSG_up_gene<-read.table("PSG_up.txt",stringsAsFactors=F)$V1

integrate_msigdb_result<-integrate_msigdb_result%>%filter(gene%in%PSG_up_gene)
msigdb_name<-read.table("msigdb_hallmark_name.txt",stringsAsFactors=F)
msigdb_name$V1<-paste("HALLMARK_",msigdb_name$V1,sep="")


integrate_msigdb_result$pathway<-factor(integrate_msigdb_result$pathway,levels=msigdb_name$V1,labels=msigdb_name$V2)

pathway_order<-levels(integrate_msigdb_result$pathway)

result<-integrate_signal_gene_result%>%filter(gene%in%PSG_up_gene)


#hallmark_result<-data.frame(PSG_up_gene);names(hallmark_result)<-"gene"
get_info<-function(target_gene){
  temp_info<-integrate_msigdb_result%>%filter(gene==target_gene)
  temp<-data.frame(data.frame(target_gene));names(temp)<-"gene"
  for(j in pathway_order){
    haha<-temp_info%>%filter(pathway==j)%>%select(gene,signal_row,ratio,pvalue)
    names(haha)[2:4]<-c("Overlap_gene_number","Odds_ratio","P_value")
    names(haha)[2:4]<-paste(j,names(haha)[2:4],sep="_")
    temp<-merge(temp,haha,by="gene")
  }
  return(temp)
}

hallmark_result<-do.call(rbind,lapply(PSG_up_gene,get_info))

result<-merge(result,hallmark_result,by="gene")

write.table(result,file="PSG_all_msigdb_hallmark.txt",row.names=F,col.names=T,sep="\t",quote=F)
