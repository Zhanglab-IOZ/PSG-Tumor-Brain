library(dplyr)
library(ggplot2)
library(reshape2)

load("integrate_info.rds")
load("tau_SPM_age_chr_info_hybrid_result.rds")
integrate_info<-merge(integrate_info,tau_SPM_age_chr_info_hybrid_result%>%select(gene,tau_tissue,phylo_age,phylo_age_class),by="gene")

PSG_data<-data.frame(rbind(integrate_info%>%filter(up_sample_num>=20,!is.na(amp_up_num),phylo_age=="Primate",tau_tissue!="broad")%>%group_by(cancer_type)%>%summarize(value=median(amp_up_num/up_sample_num))%>%mutate(class="Amplification",tau_tissue="Tissue biased",group="Up"),integrate_info%>%filter(other_sample_num>=20,!is.na(amp_other_num),phylo_age=="Primate",tau_tissue!="broad")%>%group_by(cancer_type)%>%summarize(value=median(amp_other_num/other_sample_num))%>%mutate(class="Amplification",tau_tissue="Tissue biased",group="No_up"),integrate_info%>%filter(up_sample_num>=20,!is.na(amp_up_num),phylo_age=="Primate",tau_tissue=="broad")%>%group_by(cancer_type)%>%summarize(value=median(amp_up_num/up_sample_num))%>%mutate(class="Amplification",tau_tissue="Broadly expressed",group="Up"),integrate_info%>%filter(other_sample_num>=20,!is.na(amp_other_num),phylo_age=="Primate",tau_tissue=="broad")%>%group_by(cancer_type)%>%summarize(value=median(amp_other_num/other_sample_num))%>%mutate(class="Amplification",tau_tissue="Broadly expressed",group="No_up"),integrate_info%>%filter(normal_num>=10,up_sample_num>=20,!is.na(hypo_up_num),phylo_age=="Primate",tau_tissue!="broad")%>%group_by(cancer_type)%>%summarize(value=median(hypo_up_num/up_sample_num))%>%mutate(class="Hypomethylation",tau_tissue="Tissue biased",group="Up"),integrate_info%>%filter(normal_num>=10,other_sample_num>=20,!is.na(hypo_other_num),phylo_age=="Primate",tau_tissue!="broad")%>%group_by(cancer_type)%>%summarize(value=median(hypo_other_num/other_sample_num))%>%mutate(class="Hypomethylation",tau_tissue="Tissue biased",group="No_up"),integrate_info%>%filter(normal_num>=10,up_sample_num>=20,!is.na(hypo_up_num),phylo_age=="Primate",tau_tissue=="broad")%>%group_by(cancer_type)%>%summarize(value=median(hypo_up_num/up_sample_num))%>%mutate(class="Hypomethylation",tau_tissue="Broadly expressed",group="Up"),integrate_info%>%filter(normal_num>=10,other_sample_num>=20,!is.na(hypo_other_num),phylo_age=="Primate",tau_tissue=="broad")%>%group_by(cancer_type)%>%summarize(value=median(hypo_other_num/other_sample_num))%>%mutate(class="Hypomethylation",tau_tissue="Broadly expressed",group="No_up")))

PSG_data$class<-factor(PSG_data$class,levels=c("Amplification","Hypomethylation"))

PSG_data$tau_tissue<-factor(PSG_data$tau_tissue,levels=c("Broadly expressed","Tissue biased"),labels=c("B-PSGs","T-PSGs"))

PSG_data$class<-paste(PSG_data$class,PSG_data$group,sep="_")
PSG_data$class<-factor(PSG_data$class,levels=c("Amplification_Up","Amplification_No_up","Hypomethylation_Up","Hypomethylation_No_up"))

PSG_data$group<-factor(PSG_data$group,levels=c("Up","No_up"),labels=c("Upregulated samples","Non-upregulated samples"))


pdf("PSG_tissue_broad_hypo_amp_compare.pdf",width=4.5,height=3.5)
ggplot(data=PSG_data,aes(x=class,y=value))+theme(panel.grid=element_blank(),strip.background=element_blank(),axis.ticks.x=element_blank(),legend.key=element_blank(),panel.border=element_rect(fill=NA,size=.5,colour="grey30"),panel.background=element_blank(),strip.text=element_text(size=15),axis.title=element_text(size=15),axis.text.x=element_blank(),axis.text.y=element_text(size=15),legend.title=element_blank(),legend.text=element_text(size=15),legend.position="top")+geom_violin(aes(fill=group))+geom_boxplot(width=.1,colour="black",fill="black",outlier.colour=NA)+stat_summary(fun=median,geom="point",fill="white",shape=21,size=2)+facet_grid(tau_tissue~.)+scale_fill_manual(values=c("#FFB5C5","#87CEFF"))+scale_y_continuous(limits=c(0,0.6),breaks=c(0.1,0.3,0.5))+labs(x="",y="Proportion")
#ggplot(data=PSG_data,aes(x=class,y=value))+theme(strip.background=element_blank(),axis.ticks.x=element_blank(),legend.key=element_blank(),panel.border=element_rect(fill=NA,size=.5,colour="grey30"),panel.background=element_blank(),strip.text=element_text(size=15),axis.title=element_text(size=15),axis.text.x=element_blank(),axis.text.y=element_text(size=15),legend.title=element_blank(),legend.text=element_text(size=15),legend.position="top")+geom_violin(aes(fill=group))+geom_boxplot(width=.1,colour="black",fill="black",outlier.colour=NA)+stat_summary(fun=median,geom="point",fill="white",shape=21,size=2)+facet_grid(tau_tissue~.)+scale_fill_manual(values=c("#87CEFF","#FFB5C5"))+scale_y_continuous(limits=c(0,0.6),breaks=c(0.1,0.3,0.5))+labs(x="",y="Proportion")
dev.off()


PSG_data<-dcast(PSG_data,cancer_type+tau_tissue~class,value.var="value")

data.frame(PSG_data%>%group_by(tau_tissue)%>%summarize(amp_pvalue=wilcox.test(Amplification_Up,Amplification_No_up,paired=T)$p.value,hypo_pvalue=wilcox.test(Hypomethylation_Up,Hypomethylation_No_up,paired=T)$p.value))

wilcox.test(PSG_data%>%filter(tau_tissue=="B-PSGs")%>%.$Amplification_Up,PSG_data%>%filter(tau_tissue=="B-PSGs")%>%.$Hypomethylation_Up)
wilcox.test(PSG_data%>%filter(tau_tissue=="T-PSGs")%>%.$Amplification_Up,PSG_data%>%filter(tau_tissue=="T-PSGs")%>%.$Hypomethylation_Up)

