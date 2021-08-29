library(dplyr)
library(ggplot2)
library(reshape2)


load("integrate_info.rds")
load("tau_SPM_age_chr_info_hybrid_result.rds")
integrate_info<-merge(integrate_info,tau_SPM_age_chr_info_hybrid_result%>%select(gene,phylo_age_class),by="gene")
integrate_info<-integrate_info%>%mutate(tau_tissue_class=ifelse(tau_tissue=="broad","broad","others"))
integrate_info$tau_tissue_class<-factor(integrate_info$tau_tissue_class,levels=c("broad","others"),labels=c("Broadly expressed","Tissue biased"))
integrate_info<-data.frame(integrate_info%>%group_by(cancer_type)%>%mutate(amp_rank=1-rank(amp_prop)/length(amp_prop),del_rank=1-rank(del_prop)/length(del_prop)))

PSG_data<-data.frame(integrate_info%>%filter(phylo_age=="Primate")%>%group_by(tau_tissue_class,cancer_type)%>%summarize(amp_value=median(amp_prop),del_value=median(del_prop)))
PSG_data$tau_tissue_class<-factor(PSG_data$tau_tissue_class,labels=c("B-PSGs","T-PSGs"))
GB_info<-data.frame(integrate_info%>%filter()%>%group_by(cancer_type)%>%summarize(amp_value=median(amp_prop),del_value=median(del_prop)))%>%mutate(tau_tissue_class="GB")
PSG_data<-rbind(PSG_data,GB_info)
PSG_data$tau_tissue_class<-factor(PSG_data$tau_tissue_class,levels=c("GB","B-PSGs","T-PSGs"),label=c(paste("Genome","background",sep="\n"),"B-PSGs","T-PSGs"))
data.frame(PSG_data%>%group_by(tau_tissue_class)%>%summarize(pvalue=wilcox.test(amp_value,del_value,paired=T,a="l")$p.value))
data.frame(PSG_data%>%group_by(tau_tissue_class)%>%summarize(pvalue=wilcox.test(amp_value,del_value,paired=T,a="g")$p.value))

PSG_data<-melt(PSG_data)

PSG_data$variable<-factor(PSG_data$variable,labels=c("Amplification","Deletion"))
PSG_data$class<-paste(as.character(PSG_data$tau_tissue_class),as.character(PSG_data$variable),sep="/")
PSG_data$class<-factor(PSG_data$class,levels=c("B-PSGs/Amplification","B-PSGs/Deletion","Genome\nbackground/Amplification","Genome\nbackground/Deletion","T-PSGs/Amplification","T-PSGs/Deletion"))

##Figure.3A

pdf("BG_PSG_tissue_broad_amp_del.pdf",width=6.3,height=3.5)
ggplot()+theme_classic()+theme(axis.ticks.x=element_blank(),legend.key=element_blank(),strip.text=element_text(size=15),axis.title=element_text(size=15),axis.text.x=element_blank(),axis.text.y=element_text(size=15),legend.title=element_blank(),legend.text=element_text(size=15),legend.position="top")+geom_violin(data=PSG_data,aes(x=class,y=value,fill=variable,colour=variable))+geom_boxplot(data=PSG_data,aes(x=class,y=value),width=.07,colour="black",fill="black",outlier.colour=NA)+stat_summary(data=PSG_data,aes(x=class,y=value),fun=median,geom="point",fill="white",shape=21,size=2)+scale_y_continuous(expand=c(0,0),limits=c(0,0.25),breaks=c(0,0.05,0.10,0.15,0.20))+labs(x="",y="Proportion")+scale_fill_manual(values=c("tomato2","royalblue1"))+scale_colour_manual(values=c("tomato2","royalblue1"))
dev.off()


wilcox.test(PSG_data%>%filter(tau_tissue_class=="B-PSGs",variable=="Amplification")%>%.$value,PSG_data%>%filter(tau_tissue_class=="T-PSGs",variable=="Amplification")%>%.$value,paired=T,a="g")
wilcox.test(PSG_data%>%filter(tau_tissue_class=="B-PSGs",variable=="Amplification")%>%.$value,PSG_data%>%filter(tau_tissue_class=="Genome\nbackground",variable=="Amplification")%>%.$value,paired=T,a="g")


cutoff<-0.25
cutoff<-as.numeric(cutoff)
integrate_info$tau_tissue_class<-factor(integrate_info$tau_tissue_class,labels=c("broad","others"))
amp_info<-data.frame(rbind(integrate_info%>%filter()%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(amp_prop>=cutoff),prop=num_signal/num)%>%mutate(freq_class="High",class="GB"),integrate_info%>%filter(phylo_age=="Primate",tau_tissue_class=="broad")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(amp_prop>=cutoff),prop=num_signal/num)%>%mutate(freq_class="High",class="B-PSGs"),integrate_info%>%filter(phylo_age=="Primate",tau_tissue_class=="others")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(amp_prop>=cutoff),prop=num_signal/num)%>%mutate(freq_class="High",class="T-PSGs"),integrate_info%>%filter()%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(amp_prop==0),prop=num_signal/num)%>%mutate(freq_class="No",class="GB"),integrate_info%>%filter(phylo_age=="Primate",tau_tissue_class=="broad")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(amp_prop==0),prop=num_signal/num)%>%mutate(freq_class="No",class="B-PSGs"),integrate_info%>%filter(phylo_age=="Primate",tau_tissue_class=="others")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(amp_prop==0),prop=num_signal/num)%>%mutate(freq_class="No",class="T-PSGs"),integrate_info%>%filter()%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(amp_prop>0&amp_prop<cutoff),prop=num_signal/num)%>%mutate(freq_class="Low",class="GB"),integrate_info%>%filter(phylo_age=="Primate",tau_tissue_class=="broad")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(amp_prop>0&amp_prop<cutoff),prop=num_signal/num)%>%mutate(freq_class="Low",class="B-PSGs"),integrate_info%>%filter(phylo_age=="Primate",tau_tissue_class=="others")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(amp_prop>0&amp_prop<cutoff),prop=num_signal/num)%>%mutate(freq_class="Low",class="T-PSGs")))

integrate_reshape_result<-dcast(data=amp_info,class+cancer_type~freq_class,value.var="prop")
integrate_reshape_result<-integrate_reshape_result%>%mutate(ratio=High/(High+Low))
integrate_reshape_result$class<-factor(integrate_reshape_result$class,levels=c("B-PSGs","GB","T-PSGs"),labels=c("B-PSGs","Genome background","T-PSGs"))
integrate_reshape_result$jitter_num<-jitter(as.numeric(integrate_reshape_result$class),factor=.5)

##Figure 3C
pdf("BG_PSG_tissue_broad_recurrent_amp.pdf",width=3.8,height=3.6)
ggplot()+theme(axis.ticks.x=element_blank(),legend.key=element_blank(),panel.border=element_rect(fill=NA,size=.5,colour="grey30"),panel.background=element_blank(),legend.text=element_text(size=15),legend.title=element_blank(),legend.position="top",axis.text.x=element_blank(),axis.text.y=element_text(size=15),axis.title=element_text(size=15),strip.text=element_text(size=15))+geom_violin(data=integrate_reshape_result,aes(x=class,y=ratio,fill=class,colour=class))+geom_boxplot(data=integrate_reshape_result,aes(x=class,y=ratio),width=.07,outlier.shape=NA,colour="black",fill="black")+stat_summary(data=integrate_reshape_result,aes(x=class,y=ratio),fun=median,geom="point",fill="white",shape=21,size=2)+labs(x="",y="Recurrent\namplification ratio")+scale_fill_manual(values=c("#440154","gray70","#21908C"))+scale_colour_manual(values=c("#440154","gray70","#21908C"))+scale_y_continuous(limits=c(0,0.45))
dev.off()

wilcox.test(integrate_reshape_result%>%filter(class=="B-PSGs")%>%.$ratio,integrate_reshape_result%>%filter(class=="Genome background")%>%.$ratio,paired=T,a="g")
wilcox.test(integrate_reshape_result%>%filter(class=="B-PSGs")%>%.$ratio,integrate_reshape_result%>%filter(class=="T-PSGs")%>%.$ratio,paired=T,a="g")



##based on rank
cutoff<-0.1
amp_info<-data.frame(rbind(integrate_info%>%filter()%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(amp_rank<=cutoff),prop=num_signal/num)%>%mutate(freq_class="High",class="GB"),integrate_info%>%filter(phylo_age=="Primate",tau_tissue_class=="broad")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(amp_rank<=cutoff),prop=num_signal/num)%>%mutate(freq_class="High",class="B-PSGs"),integrate_info%>%filter(phylo_age=="Primate",tau_tissue_class=="others")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(amp_rank<=cutoff),prop=num_signal/num)%>%mutate(freq_class="High",class="T-PSGs"),integrate_info%>%filter()%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(amp_prop==0),prop=num_signal/num)%>%mutate(freq_class="No",class="GB"),integrate_info%>%filter(phylo_age=="Primate",tau_tissue_class=="broad")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(amp_prop==0),prop=num_signal/num)%>%mutate(freq_class="No",class="B-PSGs"),integrate_info%>%filter(phylo_age=="Primate",tau_tissue_class=="others")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(amp_prop==0),prop=num_signal/num)%>%mutate(freq_class="No",class="T-PSGs"),integrate_info%>%filter()%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(amp_prop>0&amp_rank>cutoff),prop=num_signal/num)%>%mutate(freq_class="Low",class="GB"),integrate_info%>%filter(phylo_age=="Primate",tau_tissue_class=="broad")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(amp_prop>0&amp_rank>cutoff),prop=num_signal/num)%>%mutate(freq_class="Low",class="B-PSGs"),integrate_info%>%filter(phylo_age=="Primate",tau_tissue_class=="others")%>%group_by(cancer_type)%>%summarize(num=n(),num_signal=sum(amp_prop>0&amp_rank>cutoff),prop=num_signal/num)%>%mutate(freq_class="Low",class="T-PSGs")))


integrate_reshape_result<-dcast(data=amp_info,class+cancer_type~freq_class,value.var="prop")
integrate_reshape_result<-integrate_reshape_result%>%mutate(ratio=High/(High+Low))
integrate_reshape_result$class<-factor(integrate_reshape_result$class,levels=c("B-PSGs","GB","T-PSGs"),labels=c("B-PSGs","Genome background","T-PSGs"))
integrate_reshape_result$jitter_num<-jitter(as.numeric(integrate_reshape_result$class),factor=.5)


##Figure S3C
pdf("BG_PSG_tissue_broad_recurrent_amp_top_10.pdf",width=3.8,height=3.6)
ggplot()+theme(axis.ticks.x=element_blank(),legend.key=element_blank(),panel.border=element_rect(fill=NA,size=.5,colour="grey30"),panel.background=element_blank(),legend.text=element_text(size=15),legend.title=element_blank(),legend.position="top",axis.text.x=element_blank(),axis.text.y=element_text(size=15),axis.title=element_text(size=15),strip.text=element_text(size=15))+geom_violin(data=integrate_reshape_result,aes(x=class,y=ratio,fill=class,colour=class))+geom_boxplot(data=integrate_reshape_result,aes(x=class,y=ratio),width=.07,outlier.shape=NA,colour="black",fill="black")+stat_summary(data=integrate_reshape_result,aes(x=class,y=ratio),fun=median,geom="point",fill="white",shape=21,size=2)+labs(x="",y="Recurrent\namplification ratio")+scale_fill_manual(values=c("#440154","gray70","#21908C"))+scale_colour_manual(values=c("#440154","gray70","#21908C"))
dev.off()

wilcox.test(integrate_reshape_result%>%filter(class=="B-PSGs")%>%.$ratio,integrate_reshape_result%>%filter(class=="Genome background")%>%.$ratio,paired=T)
wilcox.test(integrate_reshape_result%>%filter(class=="B-PSGs")%>%.$ratio,integrate_reshape_result%>%filter(class=="T-PSGs")%>%.$ratio,paired=T,a="g")


