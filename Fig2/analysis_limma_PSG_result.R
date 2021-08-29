library(dplyr)
library(ggplot2)

load("integrate_limma_age_result.rds")

load("tau_SPM_age_chr_info_hybrid_result.rds")
target_gene<-read.table("expressed_gene_pancancer.txt",stringsAsFactors=F)$V1

integrate_limma_result<-merge(integrate_limma_age_result%>%select("gene","logFC","cancer_type"),tau_SPM_age_chr_info_hybrid_result,by="gene");rm(integrate_limma_age_result)


integrate_result<-data.frame(integrate_limma_result%>%filter(gene%in%target_gene,gene%in%c(tau_SPM_age_chr_info_hybrid_result%>%filter(phylo_age=="Primate")%>%.$gene))%>%group_by(cancer_type)%>%do(logFC=median(.$logFC),pvalue=wilcox.test(.$logFC,mu=0)$p.value)%>%summarize(cancer_type,logFC,pvalue))%>%arrange(desc(logFC))

integrate_result<-integrate_result%>%mutate(class=ifelse(logFC>0&pvalue<=0.05,"up",ifelse(logFC<0&pvalue<=0.05,"down","NA")))
integrate_result$class<-factor(integrate_result$class,levels=c("up","NA","down"))

integrate_result$cancer_type<-factor(integrate_result$cancer_type,levels=integrate_result$cancer_type)

pdf("limma_PSG_pancancer.pdf",width=5.5,height=5)
ggplot(data=integrate_result,aes(x=cancer_type,y=logFC,fill=class))+geom_bar(stat="identity",size=.7)+theme_bw()+theme(legend.background=element_blank(),legend.title=element_blank(),legend.text=element_text(size=15),legend.position=c(0.6,0.9),axis.text.x=element_text(size=15,angle=30,hjust=1,vjust=1),axis.title.x=element_blank(),axis.text.y=element_text(size=15),axis.title.y=element_text(size=16))+scale_fill_manual(values=c("hotpink1","grey50","cornflowerblue"),labels=c("Upregulated in tumor (P<0.05)","Not significant","Downregulated in tumor (P<0.05)"))+scale_y_continuous(limits=c(-0.3,0.45),breaks=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4),name="Log Fold Change\n(tumor vs normal)")
dev.off()
