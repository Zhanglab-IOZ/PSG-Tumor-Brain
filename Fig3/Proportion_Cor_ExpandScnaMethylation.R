library(dplyr)
library(reshape2)
library(ggplot2)
Args<-commandArgs(trailingOnly = T)

load("integrate_DE_result.rds")


p_cutoff<-Args[0]
n_cutoff<-Args[1]


## SCNA data
load("SCNA_cor_result.rds")
integrate_cor_result<-merge(integrate_cor_result,integrate_DE_result%>%select(gene,cancer_type,DE_class),by=c("gene","cancer_type"))
#integrate_cor_result<-data.frame(integrate_cor_result%>%group_by(cancer_type)%>%mutate(pearson_fdr=p.adjust(pearson_pvalue,method="fdr")))
#integrate_cor_result<-data.frame(integrate_cor_result%>%group_by(cancer_type)%>%mutate(rank_info=length(scna_rs)+1-rank(scna_rs)))


PSG_broad<-integrate_cor_result%>%filter(DE_class=="up",grepl("PSG",phylo_age_class),tau_tissue=="broad")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(scna_rs>=p_cutoff),prop=num_signal/num)%>%mutate(class="B-PSGs")
PSG_tissue<-integrate_cor_result%>%filter(DE_class=="up",grepl("PSG",phylo_age_class),tau_tissue!="broad")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(scna_rs>=p_cutoff),prop=num_signal/num)%>%mutate(class="T-PSGs")
UC<-integrate_cor_result%>%filter(DE_class=="up",grepl("UC",phylo_age_class))%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(scna_rs>=p_cutoff),prop=num_signal/num)%>%mutate(class="UC")
EM<-integrate_cor_result%>%filter(DE_class=="down",grepl("EM",phylo_age_class))%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(scna_rs>=p_cutoff),prop=num_signal/num)%>%mutate(class="EM")
scna_info<-rbind(UC,EM,PSG_broad,PSG_tissue)%>%mutate(group="scna")


## methylation
load("methylation_cor_result.rds")
integrate_cor_result<-merge(integrate_cor_result,integrate_DE_result%>%select(gene,cancer_type,DE_class),by=c("gene","cancer_type"))
#integrate_cor_result<-data.frame(integrate_cor_result%>%group_by(cancer_type)%>%mutate(pearson_fdr=p.adjust(pearson_pvalue,method="fdr")))
#integrate_cor_result<-data.frame(integrate_cor_result%>%group_by(cancer_type)%>%mutate(rank_info=rank(methylation_rs)))

PSG_broad<-integrate_cor_result%>%filter(DE_class=="up",grepl("PSG",phylo_age_class),tau_tissue=="broad")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(methylation_r<=n_cutoff),prop=num_signal/num)%>%mutate(class="B-PSGs")
PSG_tissue<-integrate_cor_result%>%filter(DE_class=="up",grepl("PSG",phylo_age_class),tau_tissue!="broad")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(methylation_r<=n_cutoff),prop=num_signal/num)%>%mutate(class="T-PSGs")
UC<-integrate_cor_result%>%filter(DE_class=="up",grepl("UC",phylo_age_class))%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(methylation_r<=n_cutoff),prop=num_signal/num)%>%mutate(class="UC")
EM<-integrate_cor_result%>%filter(DE_class=="down",grepl("EM",phylo_age_class))%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(methylation_r<=n_cutoff),prop=num_signal/num)%>%mutate(class="EM")
methylation_info<-rbind(UC,EM,PSG_broad,PSG_tissue)%>%mutate(group="methylation")

integrate_info<-rbind(scna_info,methylation_info)


integrate_info$class<-factor(integrate_info$class,levels=c("UC","EM","B-PSGs","T-PSGs"))
integrate_info$group<-factor(integrate_info$group,levels=c("scna","methylation"),labels=c("SCNA","Promoter methylation"))
integrate_info<-integrate_info%>%mutate(name=paste(class,group,sep="_"))

integrate_info$name<-factor(integrate_info$name,levels=c("UC_SCNA","UC_Promoter methylation","EM_SCNA","EM_Promoter methylation","B-PSGs_SCNA","B-PSGs_Promoter methylation","T-PSGs_SCNA","T-PSGs_Promoter methylation"))


pdf("Fig3.pdf",width=7.5,height=4)
ggplot()+theme_classic()+theme(axis.ticks.x=element_blank(),legend.key=element_blank(),strip.text=element_text(size=15),axis.title=element_text(size=15),axis.text.x=element_blank(),axis.text.y=element_text(size=15),legend.title=element_blank(),legend.text=element_text(size=15),legend.position="top")+geom_violin(data=integrate_info,aes(x=name,y=prop,fill=group,colour=group))+geom_boxplot(data=integrate_info,aes(x=name,y=prop),width=.05,colour="black",fill="black",outlier.colour=NA)+stat_summary(data=integrate_info,aes(x=name,y=prop),fun=median,geom="point",fill="white",shape=21,size=2)+labs(x="",y="Proportion")+scale_fill_manual(values=c("royalblue1","tomato2"))+scale_colour_manual(values=c("royalblue1","tomato2"))
dev.off()

temp<-dcast(data=integrate_info%>%select(cancer_type,prop,class,group),class+cancer_type~group,value.var="prop")
names(temp)[3]<-"methylation"

temp%>%group_by(class)%>%summarize(value=wilcox.test(methylation,SCNA,paired=T)$p.value)

 