library(dplyr)
library(reshape2)
library(ggplot2)


get_distance<-function(methylation_cor,scna_cor){
  methylation_cor<-as.numeric(methylation_cor)
  scna_cor<-as.numeric(scna_cor)
  temp_model<-fBasics::ks2Test(methylation_cor,scna_cor)
  D_statistic<-as.numeric(temp_model@test$statistic)
  D_value<-ifelse(D_statistic[1]==D_statistic[2],D_statistic[1],-D_statistic[1])
  return(c(D_value))
}

get_median_distance<-function(methylation_cor,scna_cor){
  methylation_cor<-as.numeric(methylation_cor)
  scna_cor<-as.numeric(scna_cor)
  D_value<-median(methylation_cor)-median(scna_cor)
  return(c(D_value))
}

load("integrate_info.rds")
load("tau_SPM_age_chr_info_hybrid_result.rds")

integrate_info<-integrate_info%>%mutate(tau_tissue_class=ifelse(tau_tissue=="broad","broad","others"))
integrate_info$phylo_age_class<-as.character(integrate_info$phylo_age_class)
integrate_info$phylo_age_class<-gsub("_early|_late","",integrate_info$phylo_age_class)
integrate_info$phylo_age_class<-integrate_info$phylo_age_class<-factor(integrate_info$phylo_age_class,levels=c("UC","EM","Vertebrate","no_info","MM","PSG"),labels=c("UC","EM","Vertebrate","no_info","MM","PSGs"))
integrate_info<-integrate_info%>%filter(normal_num>=10)
integrate_info<-data.frame(integrate_info%>%group_by(cancer_type)%>%mutate(hyper_rank=1-rank(hyper_ratio)/length(hyper_ratio),hypo_rank=1-rank(hypo_ratio)/length(hypo_ratio)))

integrate_extend_info<-integrate_info%>%filter(normal_num>=10)%>%select(gene,hyper_ratio,hypo_ratio,phylo_age,phylo_age_class,cancer_type,tau_tissue_class)
integrate_extend_info<-melt(integrate_extend_info)
integrate_extend_info$variable<-factor(integrate_extend_info$variable,labels=c("Hypermethylation","Hypomethylation"))



cutoff<-0.3
cutoff<-as.numeric(cutoff)

hypo_info<-data.frame(rbind(integrate_info%>%filter(normal_num>=10)%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hypo_ratio>=cutoff),prop=num_signal/num)%>%mutate(freq_class="High",class="GB"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="broad")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hypo_ratio>=cutoff),prop=num_signal/num)%>%mutate(freq_class="High",class="B-PSGs"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="others")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hypo_ratio>=cutoff),prop=num_signal/num)%>%mutate(freq_class="High",class="T-PSGs"),integrate_info%>%filter(normal_num>=10)%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hypo_ratio==0),prop=num_signal/num)%>%mutate(freq_class="No",class="GB"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="broad")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hypo_ratio==0),prop=num_signal/num)%>%mutate(freq_class="No",class="B-PSGs"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="others")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hypo_ratio==0),prop=num_signal/num)%>%mutate(freq_class="No",class="T-PSGs"),integrate_info%>%filter(normal_num>=10)%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hypo_ratio>0&hypo_ratio<cutoff),prop=num_signal/num)%>%mutate(freq_class="Low",class="GB"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="broad")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hypo_ratio>0&hypo_ratio<cutoff),prop=num_signal/num)%>%mutate(freq_class="Low",class="B-PSGs"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="others")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hypo_ratio>0&hypo_ratio<cutoff),prop=num_signal/num)%>%mutate(freq_class="Low",class="T-PSGs")))
hypo_info$methylation_class<-"Hypomethylation"

hyper_info<-data.frame(rbind(integrate_info%>%filter(normal_num>=10)%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hyper_ratio>=cutoff),prop=num_signal/num)%>%mutate(freq_class="High",class="GB"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="broad")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hyper_ratio>=cutoff),prop=num_signal/num)%>%mutate(freq_class="High",class="B-PSGs"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="others")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hyper_ratio>=cutoff),prop=num_signal/num)%>%mutate(freq_class="High",class="T-PSGs"),integrate_info%>%filter(normal_num>=10)%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hyper_ratio==0),prop=num_signal/num)%>%mutate(freq_class="No",class="GB"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="broad")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hyper_ratio==0),prop=num_signal/num)%>%mutate(freq_class="No",class="B-PSGs"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="others")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hyper_ratio==0),prop=num_signal/num)%>%mutate(freq_class="No",class="T-PSGs"),integrate_info%>%filter(normal_num>=10)%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hyper_ratio>0&hyper_ratio<cutoff),prop=num_signal/num)%>%mutate(freq_class="Low",class="GB"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="broad")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hyper_ratio>0&hyper_ratio<cutoff),prop=num_signal/num)%>%mutate(freq_class="Low",class="B-PSGs"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="others")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hyper_ratio>0&hyper_ratio<cutoff),prop=num_signal/num)%>%mutate(freq_class="Low",class="T-PSGs")))
hyper_info$methylation_class<-"Hypermethylation"





integrate_result<-rbind(hypo_info,hyper_info)

integrate_result$methylation_class<-factor(integrate_result$methylation_class,levels=c("Hypermethylation","Hypomethylation"))



integrate_reshape_result<-dcast(data=integrate_result,class+methylation_class+cancer_type~freq_class,value.var="prop")
integrate_reshape_result<-integrate_reshape_result%>%mutate(ratio=High/(High+Low))
integrate_reshape_result$class<-factor(integrate_reshape_result$class,levels=c("B-PSGs","GB","T-PSGs"),labels=c("B-PSGs","Genome background","T-PSGs"))
integrate_reshape_result$jitter_num<-jitter(as.numeric(integrate_reshape_result$class),factor=.5)


##Figure 3D
pdf("BG_PSG_tissue_broad_recurrent_hypo.pdf",width=3.8,height=3.6)
ggplot()+theme(axis.ticks.x=element_blank(),legend.key=element_blank(),panel.border=element_rect(fill=NA,size=.5,colour="grey30"),panel.background=element_blank(),legend.text=element_text(size=15),legend.title=element_blank(),legend.position="top",axis.text.x=element_blank(),axis.text.y=element_text(size=15),axis.title=element_text(size=15),strip.text=element_text(size=15))+geom_violin(data=integrate_reshape_result%>%filter(methylation_class=="Hypomethylation"),aes(x=class,y=ratio,fill=class,colour=class))+geom_boxplot(data=integrate_reshape_result%>%filter(methylation_class=="Hypomethylation"),aes(x=class,y=ratio),width=.07,outlier.shape=NA,colour="black",fill="black")+stat_summary(data=integrate_reshape_result%>%filter(methylation_class=="Hypomethylation"),aes(x=class,y=ratio),fun=median,geom="point",fill="white",shape=21,size=2)+labs(x="",y="Recurrent\nhypomethylation ratio")+scale_fill_manual(values=c("#440154","gray70","#21908C"))+scale_colour_manual(values=c("#440154","gray70","#21908C"))
dev.off()


pvalue_result<-dcast(integrate_reshape_result%>%filter(methylation_class=="Hypomethylation")%>%select(class,cancer_type,ratio),cancer_type~class,value.var="ratio")
names(pvalue_result)[c(2,3,4)]<-c("B_PSGs","GB","T_PSGs")

wilcox.test(pvalue_result$B_PSGs,pvalue_result$T_PSGs,paired=T,a="g")
wilcox.test(pvalue_result$T_PSGs,pvalue_result$GB,paired=T,a="g")

## rank_version
cutoff<-0.1
cutoff<-as.numeric(cutoff)

hypo_info<-data.frame(rbind(integrate_info%>%filter(normal_num>=10)%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hypo_rank<=cutoff),prop=num_signal/num)%>%mutate(freq_class="High",class="GB"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="broad")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hypo_rank<=cutoff),prop=num_signal/num)%>%mutate(freq_class="High",class="B-PSGs"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="others")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hypo_rank<=cutoff),prop=num_signal/num)%>%mutate(freq_class="High",class="T-PSGs"),integrate_info%>%filter(normal_num>=10)%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hypo_ratio==0),prop=num_signal/num)%>%mutate(freq_class="No",class="GB"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="broad")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hypo_ratio==0),prop=num_signal/num)%>%mutate(freq_class="No",class="B-PSGs"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="others")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hypo_ratio==0),prop=num_signal/num)%>%mutate(freq_class="No",class="T-PSGs"),integrate_info%>%filter(normal_num>=10)%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hypo_ratio>0&hypo_rank>cutoff),prop=num_signal/num)%>%mutate(freq_class="Low",class="GB"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="broad")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hypo_ratio>0&hypo_rank>cutoff),prop=num_signal/num)%>%mutate(freq_class="Low",class="B-PSGs"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="others")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hypo_ratio>0&hypo_rank>cutoff),prop=num_signal/num)%>%mutate(freq_class="Low",class="T-PSGs")))
hypo_info$methylation_class<-"Hypomethylation"

hyper_info<-data.frame(rbind(integrate_info%>%filter(normal_num>=10)%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hyper_rank<=cutoff),prop=num_signal/num)%>%mutate(freq_class="High",class="GB"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="broad")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hyper_rank<=cutoff),prop=num_signal/num)%>%mutate(freq_class="High",class="B-PSGs"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="others")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hyper_rank<=cutoff),prop=num_signal/num)%>%mutate(freq_class="High",class="T-PSGs"),integrate_info%>%filter(normal_num>=10)%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hyper_ratio==0),prop=num_signal/num)%>%mutate(freq_class="No",class="GB"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="broad")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hyper_ratio==0),prop=num_signal/num)%>%mutate(freq_class="No",class="B-PSGs"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="others")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hyper_ratio==0),prop=num_signal/num)%>%mutate(freq_class="No",class="T-PSGs"),integrate_info%>%filter(normal_num>=10)%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hyper_ratio>0&hyper_rank>cutoff),prop=num_signal/num)%>%mutate(freq_class="Low",class="GB"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="broad")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hyper_ratio>0&hyper_rank>cutoff),prop=num_signal/num)%>%mutate(freq_class="Low",class="B-PSGs"),integrate_info%>%filter(normal_num>=10,phylo_age=="Primate",tau_tissue_class=="others")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(hyper_ratio>0&hyper_rank>cutoff),prop=num_signal/num)%>%mutate(freq_class="Low",class="T-PSGs")))
hyper_info$methylation_class<-"Hypermethylation"

integrate_result<-rbind(hypo_info,hyper_info)
integrate_result$methylation_class<-factor(integrate_result$methylation_class,levels=c("Hypermethylation","Hypomethylation"))
integrate_reshape_result<-dcast(data=integrate_result,class+methylation_class+cancer_type~freq_class,value.var="prop")
integrate_reshape_result<-integrate_reshape_result%>%mutate(ratio=High/(High+Low))
integrate_reshape_result$class<-factor(integrate_reshape_result$class,levels=c("B-PSGs","GB","T-PSGs"),labels=c("B-PSGs","Genome background","T-PSGs"))
integrate_reshape_result$jitter_num<-jitter(as.numeric(integrate_reshape_result$class),factor=.5)


##Figure.S3D
pdf("BG_PSG_tissue_broad_recurrent_hypo_top10.pdf",width=3.8,height=3.6)

ggplot()+theme(axis.ticks.x=element_blank(),legend.key=element_blank(),panel.border=element_rect(fill=NA,size=.5,colour="grey30"),panel.background=element_blank(),legend.text=element_text(size=15),legend.title=element_blank(),legend.position="top",axis.text.x=element_blank(),axis.text.y=element_text(size=15),axis.title=element_text(size=15),strip.text=element_text(size=15))+geom_violin(data=integrate_reshape_result%>%filter(methylation_class=="Hypomethylation"),aes(x=class,y=ratio,fill=class,colour=class))+geom_boxplot(data=integrate_reshape_result%>%filter(methylation_class=="Hypomethylation"),aes(x=class,y=ratio),width=.07,outlier.shape=NA,colour="black",fill="black")+stat_summary(data=integrate_reshape_result%>%filter(methylation_class=="Hypomethylation"),aes(x=class,y=ratio),fun=median,geom="point",fill="white",shape=21,size=2)+labs(x="",y="Recurrent\nhypomethylation ratio")+scale_fill_manual(values=c("#440154","gray70","#21908C"))+scale_colour_manual(values=c("#440154","gray70","#21908C"))
dev.off()

pvalue_result<-dcast(integrate_reshape_result%>%filter(methylation_class=="Hypomethylation")%>%select(class,cancer_type,ratio),cancer_type~class,value.var="ratio")
names(pvalue_result)[c(2,3,4)]<-c("B_PSGs","GB","T_PSGs")

wilcox.test(pvalue_result$T_PSGs,pvalue_result$B_PSGs,paired=T,a="g")
wilcox.test(pvalue_result$T_PSGs,pvalue_result$GB,paired=T,a="g")


integrate_reshape_result<-dcast(data=integrate_reshape_result,cancer_type+methylation_class~class,value.var="ratio")
names(integrate_reshape_result)<-gsub("-","_",names(integrate_reshape_result))
names(integrate_reshape_result)[4]<-"GB"

wilcox.test(integrate_reshape_result%>%filter(methylation_class=="Hypomethylation")%>%.$T_PSGs,integrate_reshape_result%>%filter(methylation_class=="Hypomethylation")%>%.$B_PSGs,paired=T)

wilcox.test(integrate_reshape_result%>%filter(methylation_class=="Hypomethylation")%>%.$T_PSGs,integrate_reshape_result%>%filter(methylation_class=="Hypomethylation")%>%.$GB,paired=T)
wilcox.test(integrate_reshape_result%>%filter(methylation_class=="Hypomethylation")%>%.$B_PSGs,integrate_reshape_result%>%filter(methylation_class=="Hypomethylation")%>%.$GB,paired=T)


PSG_tissue_info<-data.frame(integrate_extend_info%>%filter(phylo_age=="Primate")%>%group_by(tau_tissue_class,variable,cancer_type)%>%summarize(value=median(value)))
PSG_tissue_info$tau_tissue_class<-factor(PSG_tissue_info$tau_tissue_class,levels=c("broad","others"),labels=c("B-PSGs","T-PSGs"))
BG_info<-data.frame(integrate_extend_info%>%filter()%>%group_by(variable,cancer_type)%>%summarize(value=median(value)))
BG_info$tau_tissue_class<-"Genome background"
BG_PSG_tissue_info<-rbind(BG_info,PSG_tissue_info)
BG_PSG_tissue_info$tau_tissue_class<-factor(BG_PSG_tissue_info$tau_tissue_class,levels=c("B-PSGs","Genome background","T-PSGs"),labels=c("B-PSGs",paste("Genome","background",sep="\n"),"T-PSGs"))

BG_PSG_tissue_info$class<-paste(as.character(BG_PSG_tissue_info$tau_tissue_class),BG_PSG_tissue_info$variable,sep="/")
BG_PSG_tissue_info$class<-factor(BG_PSG_tissue_info$class,levels=c("B-PSGs/Hypermethylation","B-PSGs/Hypomethylation","Genome\nbackground/Hypermethylation","Genome\nbackground/Hypomethylation","T-PSGs/Hypermethylation","T-PSGs/Hypomethylation"))



##Figure.3B
pdf("BG_PSG_tissue_hypo_hyper.pdf",width=4.5,height=3.5)
ggplot()+theme_classic()+theme(axis.ticks.x=element_blank(),legend.key=element_blank(),strip.text=element_text(size=15),axis.title=element_text(size=15),axis.text.x=element_blank(),axis.text.y=element_text(size=15),legend.title=element_blank(),legend.text=element_text(size=15),legend.position="top")+geom_violin(data=BG_PSG_tissue_info,aes(x=class,y=value,fill=variable,colour=variable))+geom_boxplot(data=BG_PSG_tissue_info,aes(x=class,y=value),width=.07,colour="black",fill="black",outlier.colour=NA)+stat_summary(data=BG_PSG_tissue_info,aes(x=class,y=value),fun=median,geom="point",fill="white",shape=21,size=2)+labs(x="",y="Proportion")+scale_fill_manual(values=c("#CD3278","#FB9507"))+scale_colour_manual(values=c("#CD3278","#FB9507"))
dev.off()

BG_PSG_tissue_reshape_info<-dcast(BG_PSG_tissue_info,tau_tissue_class+cancer_type~variable,value.var="value")
data.frame(BG_PSG_tissue_reshape_info%>%group_by(tau_tissue_class)%>%summarize(pvalue=wilcox.test(Hypermethylation,Hypomethylation,paired=T,a="g")$p.value))
data.frame(BG_PSG_tissue_reshape_info%>%group_by(tau_tissue_class)%>%summarize(pvalue=wilcox.test(Hypermethylation,Hypomethylation,paired=T,a="l")$p.value))

wilcox.test(BG_PSG_tissue_info%>%filter(variable=="Hypomethylation",tau_tissue_class=="T-PSGs")%>%.$value,BG_PSG_tissue_info%>%filter(variable=="Hypomethylation",tau_tissue_class=="Genome\nbackground")%>%.$value,paired=T)
wilcox.test(BG_PSG_tissue_info%>%filter(variable=="Hypomethylation",tau_tissue_class=="T-PSGs")%>%.$value,BG_PSG_tissue_info%>%filter(variable=="Hypomethylation",tau_tissue_class=="B-PSGs")%>%.$value,paired=T)
