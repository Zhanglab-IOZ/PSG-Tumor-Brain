library(dplyr)
library(ggplot2)
library(egg)

load("tau_SPM_age_chr_info_hybrid_result.rds")

tau_SPM_age_chr_info_hybrid_result<-tau_SPM_age_chr_info_hybrid_result%>%mutate(tau_tissue_class=ifelse(tau_tissue=="broad","broad","others"))

tau_SPM_age_chr_info_hybrid_result$tau_tissue_class<-factor(tau_SPM_age_chr_info_hybrid_result$tau_tissue_class,levels=c("broad","others"),labels=c("Broadly expressed","Tissue biased expressed"))

tau_SPM_age_chr_info_hybrid_result$phylo_age_class<-gsub("_early|_late","",as.character(tau_SPM_age_chr_info_hybrid_result$phylo_age_class))
tau_SPM_age_chr_info_hybrid_result$phylo_age_class<-factor(tau_SPM_age_chr_info_hybrid_result$phylo_age_class,levels=c("UC","EM","Vertebrate","no_info","MM","PSG"))

data.frame(tau_SPM_age_chr_info_hybrid_result%>%filter(!phylo_age%in%c("no_info","Vertebrate"))%>%group_by(phylo_age_class)%>%summarize(num=n()))
data.frame(tau_SPM_age_chr_info_hybrid_result%>%filter(!phylo_age%in%c("no_info","Vertebrate"))%>%group_by(phylo_age_class)%>%summarize(num=n(),sum_broad=sum(tau_tissue=="broad"),sum_tissue=sum(tau_tissue!="broad"),broad_prop=sum_broad/num,tissue_prop=sum_tissue/num))


PSG_info<-tau_SPM_age_chr_info_hybrid_result%>%filter(phylo_age=="Primate")%>%group_by(tau_tissue_class)%>%summarize(num=n())
PSG_info<-PSG_info%>%mutate(prop=num/sum(num))
PSG_plot<-ggplot()+geom_bar(data=PSG_info,aes(x="",y=prop,fill=tau_tissue_class),stat="identity", alpha=.6,size=1,width=1, color="white")+coord_polar("y",start=0)+theme_void()+scale_fill_manual(values=c("#440154","#21908C"))+theme(legend.text=element_text(size=15),legend.title=element_blank(),legend.position="none")

BG_info<-tau_SPM_age_chr_info_hybrid_result%>%filter()%>%group_by(tau_tissue_class)%>%summarize(num=n())
BG_info<-BG_info%>%mutate(prop=num/sum(num))
BG_plot<-ggplot()+geom_bar(data=BG_info,aes(x="",y=prop,fill=tau_tissue_class),stat="identity", alpha=.6,size=1,width=1, color="white")+coord_polar("y",start=0)+theme_void()+scale_fill_manual(values=c("#440154","#21908C"))+theme(legend.text=element_text(size=15),legend.title=element_blank(),legend.position="none")

ggsave(ggarrange(PSG_plot,BG_plot,nrow=1,ncol=2,draw=F),file="PSG_BG_broad_tissue_pie.pdf",width=6,height=3)

UC_info<-tau_SPM_age_chr_info_hybrid_result%>%filter(phylo_age_class=="UC")%>%group_by(tau_tissue_class)%>%summarize(num=n())
UC_info<-UC_info%>%mutate(prop=num/sum(num))
UC_plot<-ggplot()+geom_bar(data=UC_info,aes(x="",y=prop,fill=tau_tissue_class),stat="identity", alpha=.6,size=1,width=1, color="white")+coord_polar("y",start=0)+theme_void()+scale_fill_manual(values=c("#440154","#21908C"))+theme(legend.text=element_text(size=15),legend.title=element_blank(),legend.position="none")


EM_info<-tau_SPM_age_chr_info_hybrid_result%>%filter(phylo_age_class=="EM")%>%group_by(tau_tissue_class)%>%summarize(num=n())
EM_info<-EM_info%>%mutate(prop=num/sum(num))
EM_plot<-ggplot()+geom_bar(data=EM_info,aes(x="",y=prop,fill=tau_tissue_class),stat="identity", alpha=.6,size=1,width=1, color="white")+coord_polar("y",start=0)+theme_void()+scale_fill_manual(values=c("#440154","#21908C"))+theme(legend.text=element_text(size=15),legend.title=element_blank(),legend.position="none")

MM_info<-tau_SPM_age_chr_info_hybrid_result%>%filter(phylo_age_class=="MM")%>%group_by(tau_tissue_class)%>%summarize(num=n())
MM_info<-MM_info%>%mutate(prop=num/sum(num))
MM_plot<-ggplot()+geom_bar(data=MM_info,aes(x="",y=prop,fill=tau_tissue_class),stat="identity", alpha=.6,size=1,width=1, color="white")+coord_polar("y",start=0)+theme_void()+scale_fill_manual(values=c("#440154","#21908C"))+theme(legend.text=element_text(size=15),legend.title=element_blank(),legend.position="none")


ggsave(ggarrange(UC_plot,EM_plot,MM_plot,PSG_plot,nrow=2,ncol=2,draw=F),file="phylo_age_broad_tissue_pie.pdf",width=3,height=3)






#ggplot()+geom_bar(data=PSG_info,aes(x="",y=prop,fill=tau_tissue_class),stat="identity", size=2,width=1, color="white",alpha=.5)+coord_polar("y",start=0)+theme_void()+scale_fill_manual(values=c("violetred1","grey30","grey75","palegreen3"))+theme(legend.text=element_text(size=15),legend.title=element_blank(),legend.position="none")


#pdf("PSG_tissue_pie.pdf",height=3.5,width=5)
#ggplot()+geom_bar(data=PSG_info,aes(x="",y=prop,fill=tau_tissue_class),stat="identity", alpha=.6,size=2,width=1, color="white")+coord_polar("y",start=0)+theme_void()+scale_fill_manual(values=c("#440154","#21908C","#FDE725"))+theme(legend.text=element_text(size=15),legend.title=element_blank(),legend.position="right")


#ggplot()+geom_bar(data=PSG_info,aes(x="",y=prop,fill=tau_tissue_class),stat="identity", size=2,width=1, color="white")+coord_polar("y",start=0)+theme_void()+scale_fill_viridis_d()+theme(legend.text=element_text(size=15),legend.title=element_blank(),legend.position="right")
#dev.off()
