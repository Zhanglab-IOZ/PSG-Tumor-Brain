library(dplyr)
library(ggplot2)

load("R.Data")

load("tau_SPM_age_chr_info_hybrid_result_18195.rds")
load("integrate_DE_info_result_update.rds")
target_gene<-tau_SPM_age_chr_info_hybrid_result$gene

cell_cycle_gene<-tau_SPM_age_chr_info_hybrid_result%>%filter(gene_symbol%in%data$HGNC_gene_symbol)%>%.$gene


uPSG_cell_cycle<-read.table("cell_cycle_uPSG_46.txt",stringsAsFactors=F)$V1

all_up_gene<-integrate_DE_info_result%>%filter(up_class=="up")%>%.$gene

data.frame(all_gene=length(target_gene),signal_gene=length(intersect(target_gene,cell_cycle_gene)),class="Genomic background")
data.frame(all_gene=length(uPSG_cell_cycle),signal_gene=length(intersect(uPSG_cell_cycle,cell_cycle_gene)),class="uPSG_cell_cycle")

BG_info<-rbind(data.frame(all_gene=length(target_gene),signal_gene=length(intersect(target_gene,all_up_gene)),class="All genes",group="Background"),data.frame(all_gene=length(cell_cycle_gene),signal_gene=length(intersect(cell_cycle_gene,all_up_gene)),class="All genes",group="Cell cycle genes"))
PSG_info<-rbind(data.frame(all_gene=length(tau_SPM_age_chr_info_hybrid_result%>%filter(phylo_age=="Primate")%>%.$gene),signal_gene=length(intersect(tau_SPM_age_chr_info_hybrid_result%>%filter(phylo_age=="Primate")%>%.$gene,all_up_gene)),class="PSGs",group="Background"),data.frame(all_gene=length(intersect(cell_cycle_gene,tau_SPM_age_chr_info_hybrid_result%>%filter(phylo_age=="Primate")%>%.$gene)),signal_gene=length(intersect(all_up_gene,intersect(cell_cycle_gene,tau_SPM_age_chr_info_hybrid_result%>%filter(phylo_age=="Primate")%>%.$gene))),class="PSGs",group="Cell cycle genes"))
info<-rbind(BG_info,PSG_info);info<-info%>%mutate(prop=signal_gene/all_gene)
info<-info%>%mutate(category=paste(class,group,sep=" "))
info$category<-factor(info$category,levels=c("All genes Background","All genes Cell cycle genes","PSGs Background","PSGs Cell cycle genes"))


pdf("cell_cycle_enrich.pdf",width=5,height=4.2)
ggplot()+theme_classic()+theme(legend.position="top",legend.title=element_blank(),legend.text=element_text(size=15),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=15,colour="black"),axis.title.x=element_blank(),axis.title.y=element_text(size=15,colour="black"))+geom_bar(data=info,aes(x=category,y=prop,fill=group),width=.8,alpha=.6,stat="identity")+scale_y_continuous(expand=c(0,0),limits=c(0,0.82),breaks=c(0.0,0.2,0.4,0.6))+labs(x="",y="Proportion")+scale_fill_manual(values=c("#333333","#CD3278"))
dev.off()

