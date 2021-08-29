library(dplyr)
library(ggplot2)
library(reshape2)

## parse ddx11 expression across species 

four_species<-read.table("./four_species",stringsAsFactors=F)$V1

get_exp<-function(species){
  setwd(paste("./",species,"_data/",sep=""))
  sample_info<-read.table("sample_relation",stringsAsFactors=F)
  names(sample_info)<-c("sample","tissue")
  sample_info$tissue<-as.character(lapply(strsplit(sample_info$tissue,"_"),`[`,2))
  data<-data.frame(do.call(rbind,lapply(sample_info$sample,get_sample_value,species)))
  data<-merge(data,sample_info,by="sample")
  return(data)
}
get_sample_value<-function(sample,species){
  target_transcript<-as.character(read.table(paste("./",species,"_data/",species,"_reference/ddx11_select_info",sep=""))%>%filter(V3=="DDX11")%>%.$V2)
  setwd(paste("./",species,"_data/",sample,sep=""))
  
  #exp_data<-read.table("rsem.genes.results",stringsAsFactors=F,h=T)%>%select(gene_id,FPKM)
  exp_data<-read.table("rsem.isoforms.results",stringsAsFactors=F,h=T)%>%filter(transcript_id%in%target_transcript)%>%group_by(gene_id)%>%summarize(FPKM=sum(FPKM))
  names(exp_data)[1]<-"gene"
  exp_data$species<-species
  exp_data$sample<-sample
  return(exp_data)
}

integrate_data<-data.frame(do.call(rbind,lapply(four_species,get_exp)))

ddx11_exp<-integrate_data

save(ddx11_exp,file="ddx11_rsem_fpkm_sample.rds")


## compare ddx11 expression before/after duplication
load("ddx11_rsem_fpkm_sample.rds")
load("ortholog_exp.rds")


integrate_info<-integrate_info%>%mutate(log_ortholog_exp=log(FPKM+1,base=2))
ddx11_exp<-ddx11_exp%>%mutate(log_ddx11_exp=log(FPKM+1,base=2))

ortholog_exp<-integrate_info%>%group_by(sample)%>%summarize(log_ortholog_exp=median(log_ortholog_exp))
info<-merge(ddx11_exp,ortholog_exp,by="sample")


info<-info%>%mutate(class=ifelse(species%in%c("human","chimp"),"Post","Pre"))

info<-info%>%mutate(ratio=log_ddx11_exp/log_ortholog_exp)

info<-info%>%group_by(species,tissue,class)%>%summarize(ratio=median(ratio))



info$class<-factor(info$class,levels=c("Pre","Post"))

info$x_value<-jitter(as.numeric(info$class))

pdf("cross_species_ddx11.pdf",height=4,width=5)
ggplot(data=info,aes(x=class,y=ratio))+theme_classic()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.title=element_blank())+geom_violin(fill="pal
egreen3",alpha=.5)+geom_boxplot(width=.1,fill="palegreen3",outlier.colour=NA)+geom_point(aes(x=x_value,y=ratio,shape=tissue),size=3,colour="#326489")+scale_shape_manual(
  values=c(15,19,17,18))+labs(x="",y="Rescaled expression")
dev.off()


wilcox.test(info%>%filter(class=="Post")%>%.$ratio,info%>%filter(class=="Pre")%>%.$ratio,a="g")
