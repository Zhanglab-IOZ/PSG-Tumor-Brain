library(dplyr)
library(ggplot2)
library(egg)

load("integrate_DE_info_result_update.rds")
#load("integrate_DE_info_result.rds")

info<-data.frame(integrate_DE_info_result%>%group_by(phylo_age)%>%summarize(num_all=sum(cancer_type_num>=3),num_signal=sum(cancer_type_positive>=3&cancer_type_positive>3*cancer_type_negative),prop=num_signal/num_all),stringsAsFactors=F)

BG_info<-data.frame(phylo_age="BG",num_all=sum(info$num_all),num_signal=sum(info$num_signal),prop=sum(info$num_signal)/sum(info$num_all))
PSG_info<-info%>%filter(phylo_age=="Primate")
EM_info<-data.frame(phylo_age="EM",num_all=sum(info%>%filter(phylo_age%in%c("Metazoa","Eumetazoa","Bilateria","Chordata","Euteleostomi","Tetrapoda","Amniota"))%>%.$num_all),num_signal=sum(info%>%filter(phylo_age%in%c("Metazoa","Eumetazoa","Bilateria","Chordata","Euteleostomi","Tetrapoda","Amniota"))%>%.$num_signal),prop=sum(info%>%filter(phylo_age%in%c("Metazoa","Eumetazoa","Bilateria","Chordata","Euteleostomi","Tetrapoda","Amniota"))%>%.$num_signal)/sum(info%>%filter(phylo_age%in%c("Metazoa","Eumetazoa","Bilateria","Chordata","Euteleostomi","Tetrapoda","Amniota"))%>%.$num_all))

up_info<-data.frame(rbind(PSG_info,BG_info,EM_info))
up_info$phylo_age<-factor(up_info$phylo_age,levels=c("Primate","BG","EM"))
up_pvalue_PSG_BG<-as.numeric(binom.test(up_info[1,3],up_info[1,2],up_info[2,3]/up_info[2,2])$p.value)
up_pvalue_EM_BG<-as.numeric(binom.test(up_info[3,3],up_info[3,2],up_info[2,3]/up_info[2,2])$p.value)

info<-data.frame(integrate_DE_info_result%>%group_by(phylo_age)%>%summarize(num_all=sum(cancer_type_num>=3),num_signal=sum(cancer_type_negative>=3&cancer_type_negative>3*cancer_type_positive),prop=num_signal/num_all),stringsAsFactors=F)
BG_info<-data.frame(phylo_age="BG",num_all=sum(info$num_all),num_signal=sum(info$num_signal),prop=sum(info$num_signal)/sum(info$num_all))
EM_info<-data.frame(phylo_age="EM",num_all=sum(info%>%filter(phylo_age%in%c("Metazoa","Eumetazoa","Bilateria","Chordata","Euteleostomi","Tetrapoda","Amniota"))%>%.$num_all),num_signal=sum(info%>%filter(phylo_age%in%c("Metazoa","Eumetazoa","Bilateria","Chordata","Euteleostomi","Tetrapoda","Amniota"))%>%.$num_signal),prop=sum(info%>%filter(phylo_age%in%c("Metazoa","Eumetazoa","Bilateria","Chordata","Euteleostomi","Tetrapoda","Amniota"))%>%.$num_signal)/sum(info%>%filter(phylo_age%in%c("Metazoa","Eumetazoa","Bilateria","Chordata","Euteleostomi","Tetrapoda","Amniota"))%>%.$num_all))
PSG_info<-info%>%filter(phylo_age=="Primate")

down_info<-data.frame(rbind(EM_info,BG_info,PSG_info))
down_info$phylo_age<-factor(down_info$phylo_age,levels=c("Primate","BG","EM"))
down_pvalue_EM_BG<-as.numeric(binom.test(down_info[1,3],down_info[1,2],down_info[2,3]/down_info[2,2])$p.value)
down_pvalue_PSG_BG<-as.numeric(binom.test(down_info[3,3],down_info[3,2],down_info[2,3]/down_info[2,2])$p.value)



down_plot<-ggplot()+theme_bw()+theme(legend.text=element_text(size=15),legend.title=element_blank(),legend.position="top",axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=15),axis.title=element_text(size=15))+geom_bar(data=down_info,aes(x=phylo_age,y=prop,fill=phylo_age),width=.5,stat="identity")+scale_y_continuous(expand=c(0,0),position="right",breaks=c(0,0.1,0.2,0.3),limit=c(0,0.4))+labs(x="",y="Proportion of pan-cancer\ndownregulated genes")+scale_fill_manual(values=c("violetred3","grey20","palegreen4"))

up_plot<-ggplot()+theme_bw()+theme(legend.text=element_text(size=15),legend.title=element_blank(),legend.position="top",axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=15),axis.title=element_text(size=15))+geom_bar(data=up_info,aes(x=phylo_age,y=prop,fill=phylo_age),width=.5,stat="identity")+scale_y_continuous(expand=c(0,0),breaks=c(0,0.1,0.2,0.3),limit=c(0,0.4))+labs(x="",y="Proportion of pan-cancer\nupregulated genes")+scale_fill_manual(values=c("violetred3","grey20","palegreen4"))

ggsave(ggarrange(up_plot,down_plot,nrow=1,ncol=2,draw=F),file="PSG_BG_EM_all_gene.pdf",width=5.5,height=4)






