library(dplyr)
library(ggplot2)
library(reshape2)
library(egg)


load("tau_SPM_age_chr_info_hybrid_result_18195.rds")

microcephaly_gene<-read.table("microcephaly.txt",stringsAsFactors=F,h=T)%>%.$GENE_SYMBOL
microcephaly_gene<-tau_SPM_age_chr_info_hybrid_result%>%filter(gene_symbol%in%microcephaly_gene)%>%.$gene
microcephaly_info<-data.frame(table(integrate_stage_result%>%filter(median_exp>=0.2,gene%in%microcephaly_gene)%>%group_by(gene)%>%slice_max(exp)%>%filter(zscore>=1.2)%>%.$stage))
names(microcephaly_info)<-c("stage","microcephaly_num")
microcephaly_num<-length(unique(integrate_stage_result%>%filter(median_exp>=0.2,gene%in%microcephaly_gene)%>%.$gene))



macrocephaly_gene<-read.table("macrocephaly.txt",stringsAsFactors=F,h=T)%>%.$GENE_SYMBOL
macrocephaly_gene<-tau_SPM_age_chr_info_hybrid_result%>%filter(gene_symbol%in%macrocephaly_gene)%>%.$gene
macrocephaly_info<-data.frame(table(integrate_stage_result%>%filter(median_exp>=0.2,gene%in%macrocephaly_gene)%>%group_by(gene)%>%slice_max(exp)%>%filter(zscore>=1.2)%>%.$stage))
names(macrocephaly_info)<-c("stage","macrocephaly_num")
macrocephaly_num<-length(unique(integrate_stage_result%>%filter(median_exp>=0.2,gene%in%macrocephaly_gene)%>%.$gene))


DBD_gene<-unique(c(read.table("DBD_LoF_ASD.txt",stringsAsFactors=F)$V1,read.table("/rd/macy/brainspan/zhangyf_data/gencode_v18/exp_results/period/modify_13pcw/disease/DBD_missense_ASD.txt",stringsAsFactors=F)$V1))

DBD_gene<-tau_SPM_age_chr_info_hybrid_result%>%filter(gene_symbol%in%DBD_gene)%>%.$gene

DBD_info<-data.frame(table(integrate_stage_result%>%filter(median_exp>=0.2,gene%in%DBD_gene)%>%group_by(gene)%>%slice_max(exp)%>%filter(zscore>=1.2)%>%.$stage))
names(DBD_info)<-c("stage","DBD_num")
DBD_num<-length(unique(integrate_stage_result%>%filter(median_exp>=0.2,gene%in%DBD_gene)%>%.$gene))


SFARI_gene<-as.character(read.table("SFARI_update.txt",stringsAsFactors=F,sep="\t",h=T)%>%.$ensembl.id)
#SFARI_gene<-tau_SPM_age_chr_info_hybrid_result%>%filter(gene_symbol%in%SFARI_gene)%>%.$gene
SFARI_gene<-tau_SPM_age_chr_info_hybrid_result%>%filter(gene%in%SFARI_gene)%>%.$gene
SFARI_info<-data.frame(table(integrate_stage_result%>%filter(median_exp>=0.2,gene%in%SFARI_gene)%>%group_by(gene)%>%slice_max(exp)%>%filter(zscore>=1.2)%>%.$stage))
names(SFARI_info)<-c("stage","SFARI_num")
SFARI_num<-length(unique(integrate_stage_result%>%filter(median_exp>=0.2,gene%in%SFARI_gene)%>%.$gene))

TADA_gene<-read.table("TADA_gene.txt",stringsAsFactors=F,sep=",")$V1
TADA_gene<-tau_SPM_age_chr_info_hybrid_result%>%filter(gene_symbol%in%TADA_gene)%>%.$gene
TADA_info<-data.frame(table(integrate_stage_result%>%filter(median_exp>=0.2,gene%in%TADA_gene)%>%group_by(gene)%>%slice_max(exp)%>%filter(zscore>=1.2)%>%.$stage))
names(TADA_info)<-c("stage","TADA_num")
TADA_num<-length(unique(integrate_stage_result%>%filter(median_exp>=0.2,gene%in%TADA_gene)%>%.$gene))


BG_info<-data.frame(table(integrate_stage_result%>%filter(median_exp>=0.2)%>%group_by(gene)%>%slice_max(exp)%>%filter(zscore>=1.2)%>%.$stage))
names(BG_info)<-c("stage","BG_num")
BG_num<-length(unique(integrate_stage_result%>%filter(median_exp>=0.2)%>%.$gene))

info<-merge(merge(merge(merge(merge(BG_info,microcephaly_info,by="stage"),macrocephaly_info,by="stage"),DBD_info,by="stage"),SFARI_info,by="stage"),TADA_info,by="stage")

info$microcephaly_count<-microcephaly_num
info$BG_count<-BG_num
info$macrocephaly_count<-macrocephaly_num
info$DBD_count<-DBD_num
info$SFARI_count<-SFARI_num
info$TADA_count<-TADA_num

info<-info%>%mutate(PSG_prop=PSG_num/PSG_count,BG_prop=BG_num/BG_count,positive_prop=positive_num/positive_count,promoter_prop=promoter_num/promoter_count)

info<-info%>%mutate(BG_prop=BG_num/BG_count,microcephaly_prop=microcephaly_num/microcephaly_count,macrocephaly_prop=macrocephaly_num/macrocephaly_count,DBD_prop=DBD_num/DBD_count,SFARI_prop=SFARI_num/SFARI_count,TADA_prop=TADA_num/TADA_count)


info$microcephaly_pvalue<-apply(info,1,function(x){binom.test(as.numeric(x[3]),as.numeric(x[8]),as.numeric(x[2])/as.numeric(x[9]))$p.value})
info$macrocephaly_pvalue<-apply(info,1,function(x){binom.test(as.numeric(x[4]),as.numeric(x[10]),as.numeric(x[2])/as.numeric(x[9]))$p.value})
info$DBD_pvalue<-apply(info,1,function(x){binom.test(as.numeric(x[5]),as.numeric(x[11]),as.numeric(x[2])/as.numeric(x[9]))$p.value})
info$SFARI_pvalue<-apply(info,1,function(x){binom.test(as.numeric(x[6]),as.numeric(x[12]),as.numeric(x[2])/as.numeric(x[9]))$p.value})
info$TADA_pvalue<-apply(info,1,function(x){binom.test(as.numeric(x[7]),as.numeric(x[13]),as.numeric(x[2])/as.numeric(x[9]))$p.value})

info<-info%>%mutate(microcephaly_ratio=microcephaly_prop/BG_prop,macrocephaly_ratio=macrocephaly_prop/BG_prop,DBD_ratio=DBD_prop/BG_prop,SFARI_ratio=SFARI_prop/BG_prop,TADA_ratio=TADA_prop/BG_prop)

time_order<-c("P1","P2","P3","P4","P5","P6","P8","P10","P11","P12","P13","P14")
align_info<-read.table("./period_time_align",stringsAsFactors=F,sep="\t")
align_info<-align_info%>%filter(V1%in%time_order)


info$stage<-factor(info$stage,levels=time_order)
info$stage<-factor(info$stage,levels=time_order,labels=align_info$V2)

pvalue_info<-info[c(1,grep("ratio",names(info)),grep("pvalue",names(info)))]

prop_info<-info[c(1,grep("prop",names(info)))]
BG_info<-prop_info%>%select(stage,BG_prop)
prop_info<-melt(prop_info)

prop_info$variable<-factor(prop_info$variable,levels=c("BG_prop","microcephaly_prop","macrocephaly_prop","DBD_prop","SFARI_prop","TADA_prop"),labels=c("Genomic background","Microcephaly","Macrocephaly","DBD","SFARI","TADA"))



empty_bar<-c("add_1")
to_add<-add_info<-data.frame(stage=rep(empty_bar,each=2),variable=rep(levels(prop_info$variable)),value=rep(NA,2*length(empty_bar)))
colnames(to_add) <- colnames(prop_info)

add_info<-rbind(prop_info,to_add)
add_info$pathway<-factor(add_info$stage,levels=c(levels(prop_info$stage),empty_bar))
add_info$id<-as.numeric(add_info$stage)


pdf("disease_gene_distribution.pdf",width=5,height=5)
ggplot()+theme_minimal()+scale_y_reverse()+theme(legend.title=element_blank(),axis.text.x=element_text(size=15,angle=45,vjust=.5,hjust=.5),axis.text.y = element_blank(),axis.title=element_blank(),legend.position="none",panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank(),panel.grid.major.y=element_line(size=.35,colour="grey75",linetype="dashed"))+geom_bar(data=add_info,aes(x=as.numeric(id), y=value, fill=variable),alpha=.7,position="dodge",stat="identity", width=.85)+geom_point(aes(x=c(as.numeric(add_info$id)-0.22,as.numeric(add_info$id)-0.07,as.numeric(add_info$id)+0.08,as.numeric(add_info$id)+0.23,as.numeric(add_info$id)+0.38),y=-0.03),size=2.5,colour="#DBDBDB")+coord_polar()+scale_fill_manual(values=c("#999999","#87CEFF","#6CA6CD","#FF8247","#CD6839","#8B4726"))
dev.off()
