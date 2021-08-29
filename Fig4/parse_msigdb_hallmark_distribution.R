library(dplyr)
library(ggplot2)
library(reshape2)

load("integrate_msigdb_all_result.rds")
load("tau_SPM_age_chr_info_hybrid_result.rds")
#load("~/TCGA/kallisto_gencode_18_resource/TCGA/tumor/DE_united/filter_version/all_gene/integrate_DE_info_result.rds")
load("integrate_DE_info_result_update.rds")
select_hallmark<-read.table("select_hallmark",stringsAsFactors=F)
names(select_hallmark)<-c("pathway","hallmark")



integrate_msigdb_result<-merge(integrate_msigdb_result,tau_SPM_age_chr_info_hybrid_result%>%select(gene,phylo_age,tau_tissue,phylo_age_class),by="gene")

integrate_DE_info_result<-integrate_DE_info_result%>%mutate(DE_class=ifelse(cancer_type_positive>=3&cancer_type_positive>3*cancer_type_negative,"up",ifelse(cancer_type_negative>=3&cancer_type_negative>3*cancer_type_positive,"down","no")))

integrate_msigdb_result<-merge(integrate_msigdb_result,integrate_DE_info_result%>%select(gene,DE_class),by="gene")

integrate_msigdb_result[,c(3:9)]<-lapply(c(3:9),function(x){as.numeric(as.character(integrate_msigdb_result[,x]))})
integrate_msigdb_result$pathway<-as.character(integrate_msigdb_result$pathway)

integrate_msigdb_result<-integrate_msigdb_result%>%mutate(tau_tissue_class=ifelse(tau_tissue=="broad","broad",ifelse(tau_tissue=="testis","testis","others")))


gene_info<-integrate_msigdb_result%>%group_by(gene,gene_symbol,tau_tissue,tau_tissue_class,DE_class,phylo_age)%>%summarize(num=unique(BG_signal))


cell_cycle_term<-c("HALLMARK_G2M_CHECKPOINT","HALLMARK_E2F_TARGETS","HALLMARK_MITOTIC_SPINDLE")

## ratio>=1&&pvalue<=0.05;
BG_prop<-data.frame(table(integrate_msigdb_result%>%filter(ratio>=1,pvalue<=0.05)%>%.$pathway))
names(BG_prop)<-c("pathway","BG_signal_num")
BG_prop$BG_all_num<-sum(BG_prop$BG_signal_num)
BG_prop$BG_prop<-BG_prop$BG_signal_num/BG_prop$BG_all_num

PSG_up_prop<-data.frame(table(integrate_msigdb_result%>%filter(DE_class=="up",phylo_age=="Primate",ratio>=1,pvalue<=0.05)%>%.$pathway))
names(PSG_up_prop)<-c("pathway","PSG_up_signal_num")
PSG_up_prop$PSG_up_all_num<-sum(PSG_up_prop$PSG_up_signal_num)
PSG_up_prop$PSG_up_prop<-PSG_up_prop$PSG_up_signal_num/PSG_up_prop$PSG_up_all_num


all_up_prop<-data.frame(table(integrate_msigdb_result%>%filter(DE_class=="up",ratio>=1,pvalue<=0.05)%>%.$pathway))
names(all_up_prop)<-c("pathway","all_up_signal_num")
all_up_prop$all_up_all_num<-sum(all_up_prop$all_up_signal_num)
all_up_prop$all_up_prop<-all_up_prop$all_up_signal_num/all_up_prop$all_up_all_num


all_down_prop<-data.frame(table(integrate_msigdb_result%>%filter(DE_class=="down",ratio>=1,pvalue<=0.05)%>%.$pathway))
names(all_down_prop)<-c("pathway","all_down_signal_num")
all_down_prop$all_down_all_num<-sum(all_down_prop$all_down_signal_num)
all_down_prop$all_down_prop<-all_down_prop$all_down_signal_num/all_down_prop$all_down_all_num


#info<-merge(PSG_up_prop,BG_prop,by="pathway",all=T)
info<-merge(merge(PSG_up_prop,all_up_prop,by="pathway",all=T),BG_prop,by="pathway",all=T)
info$PSG_up_signal_num[is.na(info$PSG_up_signal_num)]<-0
info$PSG_up_prop[is.na(info$PSG_up_prop)]<-0
info$PSG_up_all_num[is.na(info$PSG_up_all_num)]<-unique(info$PSG_up_all_num[!is.na(info$PSG_up_all_num)])
#info<-info%>%mutate(ratio=PSG_up_prop/BG_prop)
info<-info%>%mutate(PSG_up_ratio=PSG_up_prop/BG_prop,all_up_ratio=all_up_prop/BG_prop)


info$PSG_pvalue_up<-as.numeric(apply(info,1,function(x){binom.test(as.numeric(x[2]),as.numeric(x[3]),as.numeric(x[8])/as.numeric(x[9]),a="g")$p.value}))
info$PSG_pvalue_down<-as.numeric(apply(info,1,function(x){binom.test(as.numeric(x[2]),as.numeric(x[3]),as.numeric(x[8])/as.numeric(x[9]),a="l")$p.value}))
info$PSG_pvalue<-as.numeric(apply(info,1,function(x){binom.test(as.numeric(x[2]),as.numeric(x[3]),as.numeric(x[8])/as.numeric(x[9]))$p.value}))


info$all_pvalue_up<-as.numeric(apply(info,1,function(x){binom.test(as.numeric(x[5]),as.numeric(x[6]),as.numeric(x[8])/as.numeric(x[9]),a="g")$p.value}))
info$all_pvalue_down<-as.numeric(apply(info,1,function(x){binom.test(as.numeric(x[5]),as.numeric(x[6]),as.numeric(x[8])/as.numeric(x[9]),a="l")$p.value}))
info$all_pvalue<-as.numeric(apply(info,1,function(x){binom.test(as.numeric(x[5]),as.numeric(x[6]),as.numeric(x[8])/as.numeric(x[9]))$p.value}))

info$pathway<-factor(info$pathway,levels=as.character(info%>%arrange(desc(all_up_prop),desc(BG_prop))%>%.$pathway))
raw_info<-info
info<-melt(info%>%select(pathway,PSG_up_prop,all_up_prop,BG_prop))


empty_bar<-c("add_1","add_2","add_3","add_4")
to_add<-data.frame(pathway=rep(empty_bar,each=3),variable=rep(levels(info$variable)),value=rep(NA,3*length(empty_bar)))
colnames(to_add) <- colnames(info)
add_info<-rbind(info,to_add)

add_info$pathway<-factor(add_info$pathway,levels=c(levels(info$pathway),empty_bar))

add_info$id<-as.numeric(add_info$pathway)
add_info$variable<-factor(add_info$variable,levels=c("BG_prop","all_up_prop","PSG_up_prop"))

## add the proliferation pathway info
proliferation_pathway<-c("HALLMARK_MITOTIC_SPINDLE","HALLMARK_MYC_TARGETS_V1","HALLMARK_G2M_CHECKPOINT","HALLMARK_E2F_TARGETS","HALLMARK_MYC_TARGETS_V2","HALLMARK_P53_PATHWAY")
proliferation_num<-data.frame(num=match(proliferation_pathway,levels(add_info$pathway)))


## add all up pathway
all_up_gene_pathway<-as.character(raw_info%>%mutate(fdr=p.adjust(all_pvalue_up,method="fdr"))%>%filter(fdr<=0.05)%>%.$pathway)

## add PSG up pathway
PSG_up_gene_pathway<-as.character(raw_info%>%mutate(fdr=p.adjust(PSG_pvalue_up,method="fdr"))%>%filter(fdr<=0.05)%>%.$pathway)


##Fig.4B
pdf("PSG_BG_50_pathway_enrich.pdf",width=5,height=5)
ggplot()+theme_minimal()+scale_y_reverse()+theme(legend.title=element_blank(),axis.text.x=element_text(size=15,angle=45,vjust=.5,hjust=.5),axis.text.y = element_blank(),axis.title=element_blank(),legend.position="none",panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank(),panel.grid.major.y=element_line(size=.35,colour="grey75",linetype="dashed"))+geom_bar(data=add_info,aes(x=as.numeric(id), y=value, fill=variable),alpha=0.6,position="dodge",stat="identity", width=0.8)+coord_polar()+scale_fill_manual(values=c("#333333","#009ACD","#CD3278"))+geom_segment(data=proliferation_num,aes(x=num-0.3,y=-0.008,xend=num+0.3,yend=-0.008),size=1,colour="#0C0CEA")+geom_point(data=add_info%>%filter(pathway%in%all_up_gene_pathway),aes(x=as.numeric(id),y=-0.02),size=3,colour="#009ACD",alpha=.6)+geom_point(data=add_info%>%filter(pathway%in%PSG_up_gene_pathway),aes(x=as.numeric(id),y=-0.04),size=3,colour="#CD3278",alpha=.6)
dev.off()


