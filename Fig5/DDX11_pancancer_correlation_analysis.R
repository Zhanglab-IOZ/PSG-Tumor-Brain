library(dplyr)
library(ggplot2)
Args<-commandArgs()
suppressMessages(library(ppcor))


## get correlation information
A_gene<-Args[6]

load("tau_SPM_age_chr_info_hybrid_result.rds")
A_gene_symbol<-as.character(tau_SPM_age_chr_info_hybrid_result%>%filter(gene==A_gene)%>%.$gene_symbol)
all_cancer_type<-read.table("cancer_type",stringsAsFactors=F)$V1


get_cor<-function(target_cancer_type){
  load(paste("./",target_cancer_type,"/tumor_log_result.rds",sep=""))
  tumor_exp<-tumor_final[-1];rm(tumor_final)
  load("integrate_exp_per_count.rds")
  target_gene<-intersect(read.table("expressed_gene_pancancer.txt",stringsAsFactors=F)$V1,integrate_exp_per_count%>%filter(cancer_type==target_cancer_type,zero_prop<=0.2)%>%.$gene)
  tumor_exp<-tumor_exp[target_gene,]
  target_exp<-as.numeric(tumor_exp[A_gene,])
  tumor_exp<-tumor_exp[-which(rownames(tumor_exp)==A_gene),]
  load(paste("./",target_cancer_type,"/purity_result.rds",sep=""))
  comm_sample<-intersect(purity_result$sampleID,names(tumor_exp))
  purity_result<-purity_result[match(comm_sample,purity_result$sampleID),]
  tumor_exp<-tumor_exp[match(comm_sample,names(tumor_exp))]
  
  spearman_cor<-apply(tumor_exp,1,function(x){as.numeric(spcor.test(as.numeric(x),target_exp,purity_result$purity_value,method="spearman")[1])})
  #spearman_cor<-cor(t(tumor_exp),t(tumor_exp[A_gene,]),method="spearman")
  spearman_cor<-data.frame(spearman_cor);spearman_cor$gene<-rownames(spearman_cor)
  
  pearson_cor<-apply(tumor_exp,1,function(x){as.numeric(spcor.test(as.numeric(x),target_exp,purity_result$purity_value,method="pearson")[1])})
  #pearson_cor<-cor(t(tumor_exp),t(tumor_exp[A_gene,]),method="pearson")
  pearson_cor<-data.frame(pearson_cor);pearson_cor$gene<-rownames(pearson_cor)
  
  integrate_cor<-merge(pearson_cor,spearman_cor,by="gene")
  integrate_cor$cancer_type<-target_cancer_type
  integrate_cor<-merge(integrate_cor,tau_SPM_age_chr_info_hybrid_result%>%dplyr::select("gene","gene_symbol"))
}

integrate_result<-data.frame(do.call(rbind,lapply(all_cancer_type,get_cor)))

integrate_result$target_gene<-A_gene
integrate_result$target_gene_symbol<-A_gene_symbol


save(integrate_result,file="integrate_cor_result.rds")



## parse DDX11 pan-cancer correlation results
load("integrate_cor_update_result.rds")

integrate_info<-data.frame(integrate_result%>%group_by(gene,gene_symbol)%>%summarize(num=n(),spearman_value=median(spearman_cor)))

final_info<-integrate_info%>%filter(gene_symbol=="TIMELESS"|grepl("^E2F",gene_symbol))
final_info<-integrate_info%>%filter(gene_symbol=="TIMELESS"|grepl("^E2F",gene_symbol))

final_info<-final_info%>%mutate(class=ifelse(gene_symbol=="TIMELESS","yes","no"))
final_info$rank_info<-as.numeric(unlist(lapply(final_info$spearman_value,function(x){ecdf(integrate_info$spearman_value)(x)})))

pdf("ddx11_timeless_E2Fs.pdf",width=4,height=3)
ggplot()+theme(legend.position="none",axis.text=element_text(size=10),axis.title=element_text(size=12),panel.background=element_blank(),panel.border=element_rect(fill=NA,colour="black",size=1))+stat_ecdf(data=integrate_info,aes(x=spearman_value),size=1.2,colour="grey80")+geom_point(data=final_info,aes(x=spearman_value,y=rank_info,shape=class,colour=class,size=class))+labs(x="Correlation coefficient with DDX11",y="Cumulative proportion")+scale_y_continuous(breaks=c(0,0.5,1))+geom_hline(yintercept=.9,linetype="dashed",colour="black")+scale_size_manual(values=c(1,1))+scale_colour_manual(values=c("#0C0CEA","#FF3E96"))
dev.off()
