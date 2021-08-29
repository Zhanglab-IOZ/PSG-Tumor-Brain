library(dplyr)
library(ggplot2)
library(reshape2)

## integrate the pan-cancer ssGSEA result
cancer_type<-read.table("./cancer_type",stringsAsFactors=F)$V1

rbind_all<-function(name){
  result<-lapply(name,get)
  names(result)<-name
  result<-do.call(rbind,result)
  return(result)
}
for(i in cancer_type){
  load(paste("./",i,"/ssgsea_result.rds",sep=""))
  assign(paste(i,"_result",sep=""),ssgsea_result)
  rm(i,ssgsea_result)
}
ssgsea_result<-rbind_all(paste(cancer_type,"_result",sep=""))
save(ssgsea_result,file="ssgsea_result.rds")

## rescale the ssGSEA value 
names(ssgsea_result)[5]<-"tissue"
tumor_normal_ssgsea_result<-ssgsea_result
integrate_ssgsea_result<-data.frame(rbind(tumor_normal_ssgsea_result,ESC_ssgsea_result,HPA_ssgsea_result),stringsAsFactors=F,row.names=NULL)
integrate_ssgsea_result<-merge(integrate_ssgsea_result,data.frame(integrate_ssgsea_result%>%group_by(sample,type)%>%summarize(ratio=max(ssgsea_value)-min(ssgsea_value))),by=c("sample","type"))
integrate_ssgsea_result$ssgsea_value=integrate_ssgsea_result$ssgsea_value/integrate_ssgsea_result$ratio
integrate_ssgsea_result<-subset(integrate_ssgsea_result,select=-c(ratio))

save(integrate_ssgsea_result,file="integrate_ssgsea_result.rds")

## analysis the PSGs ssGSEA result
load("ssgsea_result.rds")
names(ssgsea_result)[5]<-"tissue"
PSG_info<-ssgsea_result%>%filter(age_4_group=="Primate")
PSG_info<-PSG_info%>%mutate(scale=max(ssgsea_value)-min(ssgsea_value))
PSG_info$ssgsea_value<-PSG_info$ssgsea_value/PSG_info$scale
PSG_info<-merge(PSG_info,PSG_info%>%filter(type=="normal")%>%group_by(tissue)%>%summarize(scale_value=median(ssgsea_value)),by="tissue")
PSG_info<-data.frame(PSG_info%>%mutate(scale_ssgsea=ssgsea_value-scale_value))
cancer_type_order<-data.frame(PSG_info%>%filter(type=="tumor")%>%group_by(tissue)%>%summarize(value=median(scale_ssgsea)))%>%arrange(desc(value))
PSG_info$tissue<-factor(PSG_info$tissue,levels=cancer_type_order$tissue)
pvalue_info<-data.frame(PSG_info%>%group_by(tissue)%>%summarize(pvalue=wilcox.test(ssgsea_value[which(type=="tumor")],ssgsea_value[which(type=="normal")])$p.value,pvalue_up=wilcox.test(ssgsea_value[which(type=="tumor")],ssgsea_value[which(type=="normal")],a="g")$p.value,pvalue_down=wilcox.test(ssgsea_value[which(type=="tumor")],ssgsea_value[which(type=="normal")],a="l")$p.value))

pdf("PSG_ssgsea_violin.pdf",height=4.5,width=6.3)
ggplot()+theme(legend.position="top",legend.background=element_blank(),panel.grid=element_blank(),panel.border=element_rect(colour="black",fill=NA,size=.5),panel.background=element_blank(),legend.title=element_blank(),legend.text=element_text(size=15),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_text(size=15),axis.title=element_text(size=15),strip.text=element_text(size=13.5))+geom_violin(data=PSG_info,aes(x=type,y=scale_ssgsea,fill=type))+geom_boxplot(data=PSG_info,aes(x=type,y=scale_ssgsea),width=.07,colour="black",fill="black",outlier.colour=NA)+stat_summary(data=PSG_info,aes(x=type,y=scale_ssgsea),fun=median,geom="point",fill="white",shape=21,size=2)+facet_wrap(.~tissue,ncol=7)+scale_y_continuous(limits=c(-0.2,0.4))+scale_fill_manual(labels=c("Normal sample","Tumor sample"),values=c("mediumpurple4","goldenrod3"))+labs(x="",y="Enrichment score")
dev.off()


## analysis the phylo_age ssGSEA result 

load("integrate_ssgsea_result.rds")
tumor_normal_result<-data.frame(integrate_ssgsea_result%>%filter(type%in%c("tumor","normal"))%>%group_by(age_4_group,type,tissue)%>%summarize(ssgsea_value=median(ssgsea_value)))
tumor_normal_result<-merge(tumor_normal_result,data.frame(tumor_normal_result%>%filter(type=="normal")%>%group_by(age_4_group)%>%summarize(adjust_ssgsea_value=median(ssgsea_value))),by="age_4_group")
tumor_normal_result$rescale_ssgsea_value<-tumor_normal_result$ssgsea_value-tumor_normal_result$adjust_ssgsea_value
tumor_normal_result$type<-factor(tumor_normal_result$type,levels=c("normal","tumor"))
pvalue_result<-dcast(tumor_normal_result,age_4_group+tissue~type,value.var="rescale_ssgsea_value")
pvalue_up_result<-data.frame(pvalue_result%>%filter(!age_4_group%in%c("no_info","Vertebrate"))%>%group_by(age_4_group)%>%summarize(pvalue=wilcox.test(tumor,normal,paired=T,alternative="greater")$p.value))
write.table(pvalue_up_result,file="pvalue_up_result.txt",row.names=F,col.names=T,sep="\t",quote=F)
pvalue_up_result<-pvalue_up_result%>%filter(pvalue<=0.05)
pvalue_down_result<-data.frame(pvalue_result%>%filter(!age_4_group%in%c("no_info","Vertebrate"))%>%group_by(age_4_group)%>%summarize(pvalue=wilcox.test(tumor,normal,paired=T,alternative="less")$p.value))
write.table(pvalue_down_result,file="pvalue_down_result.txt",row.names=F,col.names=T,sep="\t",quote=F)
pvalue_down_result<-pvalue_down_result%>%filter(pvalue<=0.05)

appender <- function(string, prefix = "") paste0(prefix, string)
tumor_normal_result<-tumor_normal_result%>%filter(!age_4_group%in%c("no_info","Vertebrate"))
new_levels<-as.character(unique(tumor_normal_result$age_4_group)[as.numeric(na.omit(match(levels(tumor_normal_result$age_4_group),unique(tumor_normal_result$age_4_group))))])

tumor_normal_result$age_4_group<-factor(tumor_normal_result$age_4_group,levels=new_levels)
tumor_normal_result$age_4_group_class<-as.numeric(tumor_normal_result$age_4_group)
tumor_normal_result$age_4_group_class<-factor(tumor_normal_result$age_4_group_class,levels=c(1:length(unique(tumor_normal_result$age_4_group))))

pvalue_up_result<-pvalue_up_result%>%mutate(age_4_group_class=match(pvalue_up_result$age_4_group,levels(tumor_normal_result$age_4_group)))

pvalue_up_result$age_4_group_class<-factor(pvalue_up_result$age_4_group_class,levels=levels(tumor_normal_result$age_4_group_class))
pvalue_down_result<-pvalue_down_result%>%mutate(age_4_group_class=match(pvalue_down_result$age_4_group,levels(tumor_normal_result$age_4_group)))
pvalue_down_result$age_4_group_class<-factor(pvalue_down_result$age_4_group_class,levels=levels(tumor_normal_result$age_4_group_class))


pdf("tumor_vs_normal_phylo_age.pdf",width=11,height=4)
ggplot()+theme(legend.background=element_blank(),legend.text=element_text(size=15),legend.title=element_blank(),axis.title=element_text(size=20),axis.ticks.x=element_blank(),strip.background=element_blank(),strip.text=element_text(size=15,vjust=0.5,hjust=0.5),panel.grid=element_blank(),panel.border=element_rect(colour="grey80",fill=NA,size=1),panel.background=element_blank())+geom_boxplot(data=tumor_normal_result,aes(x=type,y=rescale_ssgsea_value,colour=type),outlier.shape=NA)+geom_text(x=1.5,y=0.095,label="*",size=8,colour="red",data=pvalue_up_result)+geom_text(x=1.5,y=0.095,label="*",size=8,colour="black",data=pvalue_down_result)+facet_grid(.~age_4_group_class,labeller = as_labeller(appender),scales="free_x",switch="x")+theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_text(size=15),axis.text.y=element_text(size=15))+xlab("")+scale_colour_manual(labels=c("Normal","Tumor"),values=c("mediumpurple4","goldenrod3"))+ylab("Scaled single sample geneset\nanalysis (ssgsea) enrichment score")
dev.off()


##PSG broadly transcribed vs tissue-biased 
load("ssgsea_result.rds")
broad_tissue_info<-ssgsea_result;rm(ssgsea_result)
ssgsea_result<-rbind(ssgsea_result%>%filter(!grepl("PSG",age_4_group)),broad_tissue_info%>%filter(grepl("PSG",age_4_group)))
ssgsea_result$age_4_group<-as.character(ssgsea_result$age_4_group)
ssgsea_result$ssgsea_value<-ssgsea_result$ssgsea_value/(max(ssgsea_result$ssgsea_value)-min(ssgsea_result$ssgsea_value))
ssgsea_result<-ssgsea_result%>%filter(grepl("PSG|UC",age_4_group))
parse_result<-data.frame(ssgsea_result%>%group_by(age_4_group,cancer_type)%>%summarize(pvalue_up=wilcox.test(ssgsea_value[which(type=="tumor")],ssgsea_value[which(type=="normal")],a="g")$p.value,pvalue_down=wilcox.test(ssgsea_value[which(type=="tumor")],ssgsea_value[which(type=="normal")],a="l")$p.value))
parse_result<-parse_result%>%mutate(class=ifelse(pvalue_up<=0.05,"up",ifelse(pvalue_down<=0.05,"down","no")))
parse_result$age_4_group<-factor(parse_result$age_4_group,levels=c("UC","PSG/tissue","PSG/broad"))
parse_result$class<-factor(parse_result$class,levels=c("up","down","no"))
pdf("PSG_broad_tissue_UC.pdf",width=7,height=4.5)
ggplot()+theme(panel.grid.major=element_line(colour="palegreen3",linetype="dotted",size=.7),panel.border=element_rect(colour="palegreen3",fill=NA),panel.background=element_blank(),axis.ticks=element_blank(),legend.key=element_blank(),axis.text.y=element_text(size=14),axis.text.x=element_text(size=14,angle=45,hjust=.5,vjust=.5),legend.title=element_blank(),legend.text=element_text(size=14),legend.position="top")+geom_point(data=parse_result,aes(x=cancer_type,y=age_4_group,colour=class),size=7)+labs(x="",y="")+scale_y_discrete(labels=c("UC","Tissue\nbiased PSGs\n(T-PSGs)","Broadly\nexpressed PSGs\n(B-PSGs)"))+scale_colour_manual(values=c("magenta","cadetblue3","gray70"),labels=c("Higher in tumor\n(P<=0.05)","Lower in tumor\n(P<=0.05)","Not significant"))
dev.off()
