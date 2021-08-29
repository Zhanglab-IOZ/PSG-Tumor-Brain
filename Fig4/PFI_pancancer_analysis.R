library(dplyr)
library(ggplot2)
library(egg)

load("integrate_DE_info_result.rds")
load("integrate_DE_result.rds")
load("integrate_cox_result.rds")

integrate_cox_result<-integrate_cox_result%>%select(gene,multi_p_value,multi_HR,multi_fdr,cancer_type)
integrate_cox_result<-integrate_cox_result%>%mutate(cox_class=ifelse(multi_HR>0&multi_fdr<=0.05,"up",ifelse(multi_HR<0&multi_fdr<=0.05,"down","no")))

load("tau_SPM_age_chr_info_hybrid_result.rds")

## pie chart
integrate_cox_result_number<-integrate_cox_result%>%group_by(gene)%>%summarize(cox_num=length(unique(cancer_type[which(cox_class!="no")])),cox_up=length(unique(cancer_type[which(cox_class=="up")])),cox_down=length(unique(cancer_type[which(cox_class=="down")])),cox_net=cox_up-cox_down)

integrate_cox_result_number<-merge(integrate_cox_result_number,tau_SPM_age_chr_info_hybrid_result%>%select(gene,gene_symbol,phylo_age,phylo_age_class),by="gene")
integrate_cox_result_number<-merge(integrate_cox_result_number,integrate_DE_info_result%>%select(gene,cancer_type_num,cancer_type_positive,cancer_type_negative),by="gene")
integrate_cox_result_number<-integrate_cox_result_number%>%mutate(DE_class=ifelse(cancer_type_positive>=3&cancer_type_positive>3*cancer_type_negative,"up",ifelse(cancer_type_negative>=3&cancer_type_negative>3*cancer_type_positive,"down","no")))
integrate_cox_result_number<-merge(integrate_cox_result_number,tau_SPM_age_chr_info_hybrid_result%>%select(gene,tau_tissue),by="gene")

PSG_up<-integrate_cox_result_number%>%filter(phylo_age=="Primate",DE_class=="up")%>%.$gene
EM_down<-integrate_cox_result_number%>%filter(phylo_age_class=="EM",DE_class=="down")%>%.$gene


UC_up<-integrate_cox_result_number%>%filter(as.numeric(phylo_age)<=3,DE_class=="up")%>%.$gene

pie_data<-rbind(data.frame(signal_num=c(nrow(integrate_cox_result_number%>%filter(gene%in%UC_up,cox_net>0)),nrow(integrate_cox_result_number%>%filter(gene%in%UC_up,cox_net<0)),nrow(integrate_cox_result_number%>%filter(gene%in%UC_up,cox_num>0,cox_net==0)),nrow(integrate_cox_result_number%>%filter(gene%in%UC_up,cox_num==0))),all_num=nrow(integrate_cox_result_number%>%filter(gene%in%UC_up)),class="UC"),data.frame(signal_num=c(nrow(integrate_cox_result_number%>%filter(gene%in%EM_down,cox_net>0)),nrow(integrate_cox_result_number%>%filter(gene%in%EM_down,cox_net<0)),nrow(integrate_cox_result_number%>%filter(gene%in%EM_down,cox_num>0,cox_net==0)),nrow(integrate_cox_result_number%>%filter(gene%in%EM_down,cox_num==0))),all_num=nrow(integrate_cox_result_number%>%filter(gene%in%EM_down)),class="EM"),data.frame(signal_num=c(nrow(integrate_cox_result_number%>%filter(cox_net>0)),nrow(integrate_cox_result_number%>%filter(cox_net<0)),nrow(integrate_cox_result_number%>%filter(cox_num>0,cox_net==0)),nrow(integrate_cox_result_number%>%filter(cox_num==0))),all_num=nrow(integrate_cox_result_number),class="BG"),data.frame(signal_num=c(nrow(integrate_cox_result_number%>%filter(gene%in%PSG_up,cox_net>0)),nrow(integrate_cox_result_number%>%filter(gene%in%PSG_up,cox_net<0)),nrow(integrate_cox_result_number%>%filter(gene%in%PSG_up,cox_num>0,cox_net==0)),nrow(integrate_cox_result_number%>%filter(gene%in%PSG_up,cox_num==0))),all_num=nrow(integrate_cox_result_number%>%filter(gene%in%PSG_up)),class="PSG"))
pie_data<-data.frame(pie_data%>%mutate(prop=signal_num/all_num))
pie_data$group<-rep(c("up","down","intermediate","no"),4)
pie_data$group<-factor(pie_data$group,levels=c("up","no","intermediate","down"))

pie_data
pie_EM<-ggplot()+geom_bar(data=pie_data%>%filter(class=="EM"),aes(x="",y=prop,fill=group),stat="identity", size=2,width=1, color="white",alpha=.5)+coord_polar("y",start=0)+theme_void()+scale_fill_manual(values=c("violetred1","grey30","grey75","palegreen3"))+theme(legend.text=element_text(size=15),legend.title=element_blank(),legend.position="none")
pie_PSG<-ggplot()+geom_bar(data=pie_data%>%filter(class=="PSG"),aes(x="",y=prop,fill=group),stat="identity", size=2,width=1, color="white",alpha=.5)+coord_polar("y",start=0)+theme_void()+scale_fill_manual(values=c("violetred1","grey30","grey75","palegreen3"))+theme(legend.text=element_text(size=15),legend.title=element_blank(),legend.position="none")
pie_BG<-ggplot()+geom_bar(data=pie_data%>%filter(class=="BG"),aes(x="",y=prop,fill=group),stat="identity", size=2,width=1, color="white",alpha=.5)+coord_polar("y",start=0)+theme_void()+scale_fill_manual(values=c("violetred1","grey30","grey75","palegreen3"))+theme(legend.text=element_text(size=15),legend.title=element_blank(),legend.position="none")
pie_UC<-ggplot()+geom_bar(data=pie_data%>%filter(class=="UC"),aes(x="",y=prop,fill=group),stat="identity", size=2,width=1, color="white",alpha=.5)+coord_polar("y",start=0)+theme_void()+scale_fill_manual(values=c("violetred1","grey30","grey75","palegreen3"))+theme(legend.text=element_text(size=15),legend.title=element_blank(),legend.position="none")
pie_blank<-ggplot()
#ggplot2::ggsave(file="PFI_fdr_0.05_UC.pdf",ggarrange(pie_UC,pie_BG,pie_PSG,pie_EM,nrow=2,ncol=2,widths=c(1,1),heights=c(1,1),draw=F),width=5.5,height=5.5,device="pdf")

ggplot2::ggsave(file="PFI_pancancer.pdf",ggarrange(pie_blank,pie_EM,pie_UC,pie_PSG,nrow=2,ncol=2,widths=c(1,1),heights=c(1,1),draw=F),width=4,height=4,device="pdf")


chi_pvalue_PSG_BG<-as.numeric(chisq.test(pie_data%>%filter(class=="PSG")%>%.$signal_num,p=c(pie_data%>%filter(class=="BG")%>%.$signal_num),rescale.p=T)%>%.$p.value)
chi_pvalue_PSG_EM<-as.numeric(chisq.test(pie_data%>%filter(class=="PSG")%>%.$signal_num,p=c(pie_data%>%filter(class=="EM")%>%.$signal_num),rescale.p=T)%>%.$p.value)
chi_pvalue_EM_BG<-as.numeric(chisq.test(pie_data%>%filter(class=="EM")%>%.$signal_num,p=c(pie_data%>%filter(class=="BG")%>%.$signal_num),rescale.p=T)%>%.$p.value)
chi_pvalue_UC_EM<-as.numeric(chisq.test(pie_data%>%filter(class=="UC")%>%.$signal_num,p=c(pie_data%>%filter(class=="EM")%>%.$signal_num),rescale.p=T)%>%.$p.value)
chi_pvalue_UC_PSG<-as.numeric(chisq.test(pie_data%>%filter(class=="PSG")%>%.$signal_num,p=c(pie_data%>%filter(class=="UC")%>%.$signal_num),rescale.p=T)%>%.$p.value)

