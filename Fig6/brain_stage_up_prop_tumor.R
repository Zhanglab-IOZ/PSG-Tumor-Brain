library(dplyr)
library(ggplot2)
library(reshape2)
library(egg)

load("integrate_stage_result.rds")



load("tau_SPM_age_chr_info_hybrid_result_18195.rds")

## cycle gene up
load("integrate_DE_info_result_update.rds")
integrate_DE_info_result<-integrate_DE_info_result%>%mutate(up_class=ifelse(cancer_type_positive>=3&cancer_type_positive>3*cancer_type_negative,"up","no_up"))

load("gene_info.rds")
cell_cycle_gene<-as.character(gene_info%>%filter(class!="NC_RNA")%>%.$gene)
cell_cycle_gene<-intersect(integrate_DE_info_result$gene,cell_cycle_gene)



cell_cycle_up_num<-length(intersect(unique(cell_cycle_gene),integrate_DE_info_result%>%filter(up_class=="up")%>%.$gene))
cell_cycle_num<-length(unique(cell_cycle_gene))
cell_cycle_up_prop=cell_cycle_up_num/cell_cycle_num



cell_cycle_info<-data.frame(table(integrate_stage_result%>%filter(median_exp>=0.2,gene%in%cell_cycle_gene)%>%group_by(gene)%>%slice_max(exp)%>%filter(zscore>=1.2)%>%.$stage))
names(cell_cycle_info)<-c("stage","cell_cycle_num")
cell_cycle_up_info<-data.frame(table(integrate_stage_result%>%filter(median_exp>=0.2,gene%in%cell_cycle_gene,gene%in%c(integrate_DE_info_result%>%filter(up_class=="up")%>%.$gene))%>%group_by(gene)%>%slice_max(exp)%>%filter(zscore>=1.2)%>%.$stage))
names(cell_cycle_up_info)<-c("stage","cell_cycle_up_num")

cell_cycle_info<-merge(cell_cycle_info,cell_cycle_up_info,by="stage")
cell_cycle_info<-cell_cycle_info%>%mutate(cell_cycle_prop=cell_cycle_up_num/cell_cycle_num)
cell_cycle_info<-cell_cycle_info%>%mutate(ratio=cell_cycle_prop/cell_cycle_up_prop)
cell_cycle_info$pvalue<-apply(cell_cycle_info,1,function(x){binom.test(as.numeric(x[3]),as.numeric(x[2]),cell_cycle_up_prop)$p.value})



time_order<-c("P1","P2","P3","P4","P5","P6","P8","P10","P11","P12","P13","P14")
align_info<-read.table("./period_time_align",stringsAsFactors=F,sep="\t")
align_info<-align_info%>%filter(V1%in%time_order)
cell_cycle_info$stage<-factor(cell_cycle_info$stage,levels=align_info$V1)



cell_cycle_draw_info<-cell_cycle_info
cell_cycle_draw_info<-cell_cycle_draw_info%>%select(stage,cell_cycle_prop)

empty_bar<-c("add_1","add_2")
#to_add<-add_info<-data.frame(stage=rep(empty_bar,each=2),variable=rep(levels(prop_info$variable)),value=rep(NA,2*length(empty_bar)))
to_add<-add_info<-data.frame(stage=rep(empty_bar,each=1),cell_cycle_prop=rep(NA,1*length(empty_bar)))
colnames(to_add) <- colnames(cell_cycle_draw_info)

add_info<-rbind(cell_cycle_draw_info,to_add)
add_info$stage<-factor(add_info$stage,levels=c(levels(cell_cycle_draw_info$stage),empty_bar))
add_info$id<-as.numeric(add_info$stage)



pdf("cell_cycle_up_brain_stage.pdf",width=5.5,height=5.5)
ggplot()+theme_minimal()+theme(legend.title=element_blank(),axis.text.x=element_text(size=15,angle=45,vjust=.5,hjust=.5),axis.text.y = element_blank(),axis.title=element_blank(),legend.position="none",panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank(),panel.grid.major.y=element_blank())+geom_bar(data=add_info,aes(x=as.numeric(id), y=cell_cycle_prop),alpha=.7,position="dodge",stat="identity", width=.7,fill="#FDE725")+coord_polar()+ylim(c(-0.23,0.65))+geom_segment(aes(x = 0.5, y = cell_cycle_up_prop, xend = nrow(add_info)+0.5, yend = cell_cycle_up_prop), linetype="dashed",colour = "#000000", alpha=.65, size=1.5 , inherit.aes = FALSE )+geom_segment(aes(x = 0.5, y = 0.4, xend = nrow(add_info)+0.5, yend = 0.4), linetype="dashed",colour = "#BFBFBF", alpha=1, size=0.5 , inherit.aes = FALSE )+geom_segment(aes(x = 0.5, y = 0, xend = nrow(add_info)+0.5, yend = 0), linetype="dashed",colour = "#BFBFBF", alpha=1, size=0.5 , inherit.aes = FALSE )+geom_segment(aes(x = 0.5, y = 0.2, xend = nrow(add_info)+0.5, yend = 0.2), linetype="dashed",colour = "#BFBFBF", alpha=1, size=0.5 , inherit.aes = FALSE )+geom_segment(aes(x = 0.5, y = 0.6, xend = nrow(add_info)+0.5, yend = 0.6), linetype="dashed",colour = "#BFBFBF", alpha=1, size=0.5 , inherit.aes = FALSE )
dev.off()


## up prop in each stage
load("integrate_DE_info_result_update.rds")
integrate_DE_info_result<-integrate_DE_info_result%>%mutate(up_class=ifelse(cancer_type_positive>=3&cancer_type_positive>3*cancer_type_negative,"up","no_up"))

BG_up_num<-length(intersect(unique(integrate_stage_result%>%filter(median_exp>=0.2)%>%.$gene),integrate_DE_info_result%>%filter(up_class=="up")%>%.$gene))
BG_num<-length(unique(integrate_stage_result%>%filter(median_exp>=0.2)%>%.$gene))
BG_up_prop=BG_up_num/BG_num

BG_info<-data.frame(table(integrate_stage_result%>%filter(median_exp>=0.2)%>%group_by(gene)%>%slice_max(exp)%>%filter(zscore>=1.2)%>%.$stage))
names(BG_info)<-c("stage","BG_num")
BG_up_info<-data.frame(table(integrate_stage_result%>%filter(median_exp>=0.2,gene%in%c(integrate_DE_info_result%>%filter(up_class=="up")%>%.$gene))%>%group_by(gene)%>%slice_max(exp)%>%filter(zscore>=1.2)%>%.$stage))
names(BG_up_info)<-c("stage","BG_up_num")

BG_info<-merge(BG_info,BG_up_info,by="stage")
BG_info<-rbind(BG_info%>%filter(stage=="P1")%>%mutate(stage="Embryonic"),data.frame(stage="Other_stage",BG_num=sum(BG_info%>%filter(stage!="P1")%>%.$BG_num),BG_up_num=sum(BG_info%>%filter(stage!="P1")%>%.$BG_up_num)))
BG_info<-BG_info%>%mutate(BG_prop=BG_up_num/BG_num)
BG_info<-BG_info%>%mutate(class="All")

BG_info<-rbind(BG_info,data.frame(stage="Genomic background",BG_num=BG_num,BG_up_num=BG_up_num,BG_prop=BG_up_num/BG_num,class="All"))

BG_info$stage<-factor(BG_info$stage,levels=c("Genomic background","Embryonic","Other_stage"),labels=c("Genomic background","Embryonic stage biased","Other stage biased"))

pdf("stage1_other_BG_up.pdf",width=4.5,height=4.5)
ggplot(data=BG_info)+theme_classic()+theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=15,colour="black"),axis.title=element_text(size=15,colour="black"),legend.title=element_blank(),legend.text=element_text(size=15),legend.position="none")+geom_bar(aes(x=stage,y=BG_prop,fill=stage),alpha=.7,stat="identity",width=.7)+scale_fill_manual(values=c("#999999","#87CEEB","#FDE725"))+labs(x="",y="Proportion")+scale_y_continuous(expand=c(0,0),limits=c(0,0.35))
dev.off()