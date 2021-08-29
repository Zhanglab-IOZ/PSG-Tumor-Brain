library(ggplot2);library(dplyr)
Args<-commandArgs()

target_cancer_type<-Args[6]

dir.create(target_cancer_type)
setwd(target_cancer_type)


load(paste(".",target_cancer_type,"/mydata_TSS_tumor.RData",sep=""))
load(paste(".",target_cancer_type,"/mydata_TSS_normal.RData",sep=""))


rownames(mydata_TSS_tumor)<-mydata_TSS_tumor$gene
rownames(mydata_TSS_normal)<-mydata_TSS_normal$gene
mydata_TSS_tumor<-mydata_TSS_tumor[-1]
mydata_TSS_normal<-mydata_TSS_normal[-1]

get_info<-function(genename,tumor_data,normal_data){
  tumor_methylation<-as.numeric(tumor_data[genename,])
  normal_methylation<-as.numeric(normal_data[genename,])
  normal_methylation<-normal_methylation[!is.na(normal_methylation)]
  tumor_methylation<-tumor_methylation[!is.na(tumor_methylation)]
  normal_num<-length(normal_methylation);
  normal_up_cutoff<-median(normal_methylation)+2*sd(normal_methylation)
  normal_down_cutoff<-median(normal_methylation)-2*sd(normal_methylation)
  hyper_ratio<-round(sum(tumor_methylation>normal_up_cutoff)/length(tumor_methylation),4)
  hypo_ratio<-round(sum(tumor_methylation<normal_down_cutoff)/length(tumor_methylation),4)
  return(list(genename,normal_num,hyper_ratio,hypo_ratio))
}

gene_list<-rownames(mydata_TSS_tumor)
gene_list<-intersect(gene_list,read.table("expressed_gene_pancancer.txt",stringsAsFactors=F)$V1)

integrate_info<-lapply(gene_list,get_info,mydata_TSS_tumor,mydata_TSS_normal)
integrate_info<-data.frame(do.call(rbind,integrate_info))
names(integrate_info)<-c("gene","normal_num","hyper_ratio","hypo_ratio")
#temp$gene<-as.character(temp$gene)
integrate_info[,2:4]<-lapply(2:4,function(x){as.numeric(as.character(integrate_info[,x]))})
integrate_info[,1]<-lapply(1,function(x){as.character(integrate_info[,x])})
integrate_info$cancer_type<-target_cancer_type

load("tau_SPM_age_chr_info_hybrid_result.rds")
integrate_info<-merge(integrate_info,tau_SPM_age_chr_info_hybrid_result,by="gene")

save(integrate_info,file="integrate_info.rds")
