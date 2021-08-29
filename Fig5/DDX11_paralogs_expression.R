library(dplyr)
library(reshape2)
library(ggplot2)

species_name<-c("human","chimp")


get_exp<-function(name,gene_info){
  file_name<-paste("./",name,"/rsem.isoforms.results",sep="")
  exp_data<-read.table(file_name,stringsAsFactors=F,h=T,sep="\t");
  exp_data<-merge(exp_data,gene_info,by="transcript_id")
  exp_result<-data.frame(exp_data%>%group_by(gene)%>%summarize(exp=sum(FPKM)));
  exp_result$sample<-name
  return(exp_result)
}


get_info<-function(species_name){
  dir<-paste("./",species_name,"_data",sep="")
  setwd(dir)
  sample_info<-read.table("sample_relation",stringsAsFactors=F);names(sample_info)<-c("sample","info")
  sample_info$tissue<-as.character(lapply(strsplit(sample_info$info,"_"),`[`,2))
  gene_info<-read.table(paste("./",species_name,"_reference/ddx11_select_info",sep=""))
  gene_info<-gene_info[c(2,3)]
  names(gene_info)<-c("transcript_id","gene")

  integrate_exp<-do.call(rbind,lapply(sample_info$sample,get_exp,gene_info))

  integrate_exp<-merge(integrate_exp,sample_info,by="sample")

  #integrate_exp<-integrate_exp%>%group_by(tissue,gene)%>%summarize(exp=median(exp))
  integrate_exp$species<-species_name
  integrate_exp<-data.frame(integrate_exp)
  return(integrate_exp)
}


ddx11_exp<-data.frame(do.call(rbind,lapply(species_name,get_info)))

ddx11_exp$gene<-as.character(ddx11_exp$gene)

save(ddx11_exp,file="ddx11_compare_rsem_fpkm.rds")


select_ddx11_exp<-ddx11_exp%>%filter(species=="human")

#select_ddx11_exp<-select_ddx11_exp%>%mutate(class=paste(species,gene,sep="_"))

#select_ddx11_exp$class<-factor(select_ddx11_exp$class,levels=c("mouse_DDX11","rhesus_DDX11","chimp_DDX12P","human_DDX12P","chimp_DDX11","human_DDX11"))

#select_ddx11_exp<-select_ddx11_exp%>%filter(grepl("DDX11",class))
select_ddx11_exp<-select_ddx11_exp%>%mutate(log_exp=log((exp+1),base=2))
#select_ddx11_exp$species<-factor(select_ddx11_exp$species,levels=c("human","chimp"),labels=c("Human","Chimpanzee"))
select_ddx11_exp$gene<-factor(select_ddx11_exp$gene,levels=c("AC009533.1","DDX12P","DDX11"),labels=c("LOC642846","DDX12P","DDX11"))
#select_ddx11_exp<-select_ddx11_exp%>%mutate(log_exp=log((exp+1),base=2))

#select_ddx11_exp<-select_ddx11_exp%>%group_by(tissue)%>%mutate(relative_exp=log_exp/max(log_exp))

select_ddx11_exp$x_value<-jitter(as.numeric(select_ddx11_exp$gene),factor=1)

pdf("DDX11_compare_rsem_fpkm.pdf",width=5,height=4)
ggplot(data=select_ddx11_exp,aes(x=gene,y=log_exp))+theme_classic()+theme(legend.title=element_blank())+geom_violin(fill="palegreen3",alpha=.5)+geom_boxplot(width=.1,fill="palegreen3",outlier.colour=NA)+geom_point(aes(x=x_value,y=log_exp,shape=tissue),size=3,colour="#326489")+scale_shape_manual(values=c(15,19,17,18))+labs(x="",y="Log2(FPKM+1)")
dev.off()


temp<-dcast(select_ddx11_exp,sample~gene,value.var="log_exp")

wilcox.test(temp$LOC642846,temp$DDX11,paired=T,a="l")
wilcox.test(temp$DDX12P,temp$DDX11,paired=T,a="l")

