library(dplyr);library(ggplot2)


Args<-commandArgs()

target_cancer_type<-Args[6]

dir.create(target_cancer_type)
setwd(target_cancer_type)

all_data<-read.table(paste("./",target_cancer_type,"/all_data_by_genes.txt",sep=""),stringsAsFactors=F,sep="\t",h=T)
all_data<-all_data[-c(2:3)]
names(all_data)[1]<-"gene_symbol"
names(all_data)[2:ncol(all_data)]<-strtrim(names(all_data)[2:ncol(all_data)],12)
names(all_data)[2:ncol(all_data)]<-gsub("\\.","_",names(all_data)[2:ncol(all_data)])

load("tau_SPM_age_chr_info_hybrid_result.rds")
tau_SPM_age_chr_info_hybrid_result<-tau_SPM_age_chr_info_hybrid_result%>%filter(gene%in%read.table("~/TCGA/kallisto_gencode_18_resource/TCGA/tumor/exp_summary/filter_version/expressed_gene_pancancer.txt",stringsAsFactors=F)$V1)
load(paste("./",target_cancer_type,"/tumor_log_result.rds",sep=""))

all_data<-cbind(all_data[1],all_data[intersect(names(tumor_final),names(all_data))])
rownames(all_data)<-all_data$gene_symbol
all_data<-all_data[-1]

temp<-data.frame(t(apply(all_data,1,function(x){del_prop<-round(sum(x<=-0.2)/length(x),4);amp_prop<-round(sum(x>=0.2)/length(x),4);return(c(del_prop,amp_prop))})))
names(temp)<-c("del_prop","amp_prop")
temp$gene_symbol<-rownames(temp)

integrate_info<-merge(temp,tau_SPM_age_chr_info_hybrid_result,by="gene_symbol")
integrate_info$cancer_type<-target_cancer_type

save(integrate_info,file="integrate_info.rds")
