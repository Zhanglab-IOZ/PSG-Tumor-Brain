library(dplyr)
library(ggplot2)
library(reshape2)
library(egg)

load("integrate_stage_result.rds")

## stage info
PSG_info<-data.frame(table(integrate_stage_result%>%filter(median_exp>=0.2,phylo_age=="Primate")%>%group_by(gene)%>%slice_max(exp)%>%filter(zscore>=1.2)%>%.$stage))
names(PSG_info)<-c("stage","PSG_num")
PSG_num<-length(unique(integrate_stage_result%>%filter(median_exp>=0.2,phylo_age=="Primate")%>%.$gene))

BG_info<-data.frame(table(integrate_stage_result%>%filter(median_exp>=0.2)%>%group_by(gene)%>%slice_max(exp)%>%filter(zscore>=1.2)%>%.$stage))
names(BG_info)<-c("stage","BG_num")
BG_num<-length(unique(integrate_stage_result%>%filter(median_exp>=0.2)%>%.$gene))


load("GR_data.rds")
load("tau_SPM_age_chr_info_hybrid_result_18195.rds")


## GR positive gene
positive_gene<-GR_data%>%filter(OmegaGC12>=-2)%>%.$gene_symbol;positive_gene<-tau_SPM_age_chr_info_hybrid_result%>%filter(gene_symbol%in%positive_gene)%>%.$gene
positive_info<-data.frame(table(integrate_stage_result%>%filter(median_exp>=0.2,gene%in%positive_gene)%>%group_by(gene)%>%slice_max(exp)%>%filter(zscore>=1.2)%>%.$stage))
names(positive_info)<-c("stage","positive_num")
positive_num<-length(unique(integrate_stage_result%>%filter(median_exp>=0.2,gene%in%positive_gene)%>%.$gene))


## NG promoter positive gene
data<-read.table("/rd/macy/brainspan/zhangyf_data/gencode_v18/exp_results/period/modify_13pcw/disease/human_pass.txt",stringsAsFactors=F,h=T,sep="\t")
gene_symbol<-data%>%filter(boot.p.vals.med<=0.05)%>%.$HGNC.s.
select_gene_symbol<-gene_symbol[grep("[a-zA-Z]",gene_symbol)]

promoter_info<-data.frame(table(integrate_stage_result%>%filter(median_exp>=0.2,gene_symbol%in%select_gene_symbol)%>%group_by(gene)%>%slice_max(exp)%>%filter(zscore>=1.2)%>%.$stage))
names(promoter_info)<-c("stage","promoter_num")
promoter_num<-length(unique(integrate_stage_result%>%filter(median_exp>=0.2,gene_symbol%in%select_gene_symbol)%>%.$gene))


info<-merge(merge(merge(PSG_info,BG_info,by="stage"),positive_info,by="stage"),promoter_info,by="stage")
info$PSG_count<-PSG_num
info$BG_count<-BG_num
info$positive_count<-positive_num
info$promoter_count<-promoter_num

#info<-info%>%mutate(ratio=(PSG_num/PSG_count)/(BG_num/BG_count))%>%arrange(ratio)
info<-info%>%mutate(PSG_ratio=(PSG_num/PSG_count)/(BG_num/BG_count),positive_ratio=(positive_num/positive_count)/(BG_num/BG_count),promoter_ratio=(promoter_num/promoter_count)/(BG_num/BG_count))

info$PSG_pvalue<-apply(info,1,function(x){binom.test(as.numeric(x[2]),as.numeric(x[6]),as.numeric(x[3])/as.numeric(x[7]))$p.value})
info$positive_pvalue<-apply(info,1,function(x){binom.test(as.numeric(x[4]),as.numeric(x[8]),as.numeric(x[3])/as.numeric(x[7]))$p.value})
info$promoter_pvalue<-apply(info,1,function(x){binom.test(as.numeric(x[5]),as.numeric(x[9]),as.numeric(x[3])/as.numeric(x[7]))$p.value})

time_order<-c("P1","P2","P3","P4","P5","P6","P8","P10","P11","P12","P13","P14")
align_info<-read.table("./period_time_align",stringsAsFactors=F,sep="\t")
align_info<-align_info%>%filter(V1%in%time_order)


info$stage<-factor(info$stage,levels=time_order)
info$stage<-factor(info$stage,levels=time_order,labels=align_info$V2)
info<-info%>%mutate(PSG_prop=PSG_num/PSG_count,BG_prop=BG_num/BG_count,positive_prop=positive_num/positive_count,promoter_prop=promoter_num/promoter_count)
info<-melt(info%>%select(stage,BG_prop,PSG_prop,positive_prop,promoter_prop))
info$variable<-factor(info$variable,levels=c("BG_prop","PSG_prop","positive_prop","promoter_prop"),labels=c("Genomic background","PSGs","Positively selected genes\n(Coding sequence)","Positively selected genes\n(Promoter)"))

info<-rbind(info,info%>%filter(variable=="Genomic background")%>%mutate(variable="Blank"))
info$variable<-factor(info$variable,levels=c("Genomic background","Blank","PSGs","Positively selected genes\n(Coding sequence)","Positively selected genes\n(Promoter)"))


pdf("brain_stage_exp_selection_disease.pdf",width=13,height=6)
ggplot(data=info,aes(x=stage,y=value))+theme_classic()+theme(axis.ticks.x=element_blank(),legend.position=c(0.8,0.8),legend.title=element_blank(),legend.text=element_text(size=15),axis.text.x=element_text(size=15,angle=45,vjust=.5,hjust=.5),axis.text.y=element_text(size=15),axis.title=element_text(size=15))+geom_bar(aes(fill=variable),alpha=.7,stat="identity",position="dodge")+scale_y_continuous(expand=c(0,0))+labs(x="",y="Proportion")+scale_fill_manual(values=c("#999999","#999999","#CDAD00","#FFC125","#8B6914"))
dev.off()


