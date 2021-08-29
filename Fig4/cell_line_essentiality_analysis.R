library(dplyr)
library(ggplot2)

load("process_summary.rds")
comm_essential_gene<-read.table("common_essentials_gene_symbol.csv",stringsAsFactors=F)$V1
parse_result<-gene_per_cell_line_result%>%filter(cancer_type%in%cancer_type_select)%>%group_by(gene,gene_symbol,tau_tissue,phylo_age)%>%summarize(num=sum(rank_info<=0.117))

select_gene<-read.table("./select_gene.txt",stringsAsFactors=F)$V1

load("integrate_DE_info_result.rds")
PSG_up<-integrate_DE_info_result%>%filter(cancer_type_positive>=3,cancer_type_positive>3*cancer_type_negative,phylo_age=="Primate")%>%.$gene
load("integrate_DE_info_result.rds")
up_gene<-integrate_DE_info_result%>%filter(cancer_type_positive>=3,cancer_type_positive>3*cancer_type_negative)%>%.$gene


integrate_info<-rbind(parse_result%>%select(gene,num)%>%mutate(class="BG"),parse_result%>%filter(gene%in%up_gene)%>%select(gene,num)%>%mutate(class="UP"),parse_result%>%filter(gene_symbol%in%select_gene)%>%select(gene,num)%>%mutate(class="Oncogenic_PSGs"))
integrate_info$class<-factor(integrate_info$class,levels=c("Oncogenic_PSGs","UP","BG"))

integrate_info<-integrate_info%>%mutate(log_num=log(num+1,base=2))


pdf("dependent_log_cell_line.pdf",height=4,width=4.5)
ggplot(data=integrate_info,aes(x=class,y=log_num))+theme_classic()+theme(axis.title=element_text(size=15),legend.position="none",axis.text=element_text(size=15))+geom_violin(aes(colour=class,fill=class))+geom_boxplot(width=.1,colour="black",fill="black",outlier.colour=NA)+stat_summary(fun=median,geom="point",fill="white",shape=21,size=2)+scale_x_discrete(labels=c("Oncogenic\nPSGs","Upregulated\ngenes","Genomic\nbackground"))+labs(x="",y="log2(cell line number+1)")+scale_fill_manual(values=c("violetred3","deepskyblue3","grey20"))+scale_colour_manual(values=c("violetred3","deepskyblue3","grey20"))
dev.off()

wilcox.test(integrate_info%>%filter(class=="Oncogenic_PSGs")%>%.$num,integrate_info%>%filter(class=="UP")%>%.$num)$p.value
wilcox.test(integrate_info%>%filter(class=="BG")%>%.$num,integrate_info%>%filter(class=="UP")%>%.$num)$p.value

wilcox.test(integrate_info%>%filter(class=="Oncogenic_PSGs")%>%.$log_num,integrate_info%>%filter(class=="UP")%>%.$log_num)$p.value
wilcox.test(integrate_info%>%filter(class=="BG")%>%.$log_num,integrate_info%>%filter(class=="UP")%>%.$log_num)$p.value

cutoff_essentials<-as.numeric(quantile(rank_summary$value,1031/16773))


final_info<-data.frame(cbind(select_gene,as.numeric(rank_summary%>%filter(gene_symbol%in%select_gene)%>%.$value)))
names(final_info)[2]<-"value"
final_info$select_gene<-as.character(final_info$select_gene)
final_info$value<-as.numeric(as.character(final_info$value))
final_info$pvalue<-ecdf(rank_summary$value)(final_info$value)

density_info<-data.frame(cbind(density(rank_summary$value)$x,density(rank_summary$value)$y))
names(density_info)<-c("x_value","y_value")

final_info$yvalue<-as.numeric(lapply(final_info$value,function(x){density_info$y_value[which.min(abs(density_info$x_value-x))]}))

final_info<-final_info%>%mutate(class=ifelse(value<cutoff_essentials,"yes","no"))

pdf("DDX11_common_essentials.pdf",height=4,width=5)
ggplot()+theme_classic()+theme(legend.position="none",panel.grid=element_blank(),axis.text=element_text(size=15),axis.title=element_text(size=15))+geom_line(data=rank_summary,aes(x=value),stat="density")+geom_point(data=final_info,aes(x=value,y=yvalue,shape=class,colour=class,size=class))+scale_size_manual(values=c(2,3))+scale_colour_manual(values=c("#0C0CEA","#FF3E96"))+geom_vline(xintercept=cutoff_essentials,colour="red")+labs(x="Rank of gene in 90th percentile\nleast dependent cell line",y="Density")

dev.off()


#target_pvalue
