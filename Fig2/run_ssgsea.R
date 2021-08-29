source("http://bioconductor.org/biocLite.R")
library("GSVA")
library("reshape2")
library("psych")
library("dplyr")

Args<-commandArgs()
cancer_type<-Args[6]

dir.create(cancer_type)
setwd(cancer_type)


##ssgsea gene list
load("tau_SPM_age_chr_info_hybrid_result.rds")
tau_SPM_age_chr_info_result<-tau_SPM_age_chr_info_hybrid_result
target_gene<-read.table("expressed_gene_pancancer.txt",stringsAsFactors=F)$V1
tau_SPM_age_chr_info_result<-tau_SPM_age_chr_info_result%>%filter(gene%in%target_gene)
target_gene<-tau_SPM_age_chr_info_result$gene


gene_list<-lapply(levels(tau_SPM_age_chr_info_result$phylo_age),function(i){tau_SPM_age_chr_info_result%>%filter(phylo_age==i)%>%.$gene})
names(gene_list)<-levels(tau_SPM_age_chr_info_result$phylo_age)




#get_tumor_normal_ssgsea<-function(cancer_type){
load(paste("./",cancer_type,"/tumor_log_result.rds",sep=""))
tumor_expression<-tumor_final[-1];rm(tumor_final)
tumor_expression<-tumor_expression[target_gene,]
load(paste("./",cancer_type,"/purity_result.rds",sep=""))
tumor_expression<-tumor_expression[intersect(purity_result%>%filter(purity_value>=0.7)%>%.$sampleID,names(tumor_expression))]
#run gsva package
tumor_result<-gsva(as.matrix(tumor_expression), gene_list, method="ssgsea", verbose=FALSE, parallel.sz=5,kcdf="Gaussian",abs.ranking=F,mx.diff=F,ssgsea.norm=F)
tumor_result<-data.frame(melt(t(tumor_result)))
names(tumor_result)<-c("sample","age_4_group","ssgsea_value")
tumor_result$type<-"tumor"

load(paste("./",cancer_type,"/normal_log_result.rds",sep=""))
normal_expression<-normal_expression[-1]
normal_expression<-normal_expression[target_gene,]
normal_result<-gsva(as.matrix(normal_expression), gene_list, method="ssgsea", verbose=FALSE, parallel.sz=5,kcdf="Gaussian",abs.ranking=F,mx.diff=T,ssgsea.norm=F)
normal_result<-data.frame(melt(t(normal_result)))
names(normal_result)<-c("sample","age_4_group","ssgsea_value")
normal_result$type<-"normal"
integrate_result<-rbind(tumor_result,normal_result)
integrate_result$cancer_type<-cancer_type
ssgsea_result<-integrate_result
save(ssgsea_result,file="ssgsea_result.rds")
