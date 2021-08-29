library(dplyr)
library(ggplot2)
library(scatterpie)
library(egg)


load("process_summary.rds")
parse_info<-data.frame(gene_per_cell_line_result%>%filter(cancer_type%in%cancer_type_select)%>%group_by(gene,gene_symbol,tau_tissue,phylo_age,cancer_type)%>%summarize(num=sum(rank_info<=0.117),num_all=n(),prop=num/num_all))

integrate_limma_result<-read.table("final_info.txt",stringsAsFactors=F,h=T)
names(integrate_limma_result)[4:5]<-c("limma_logFC","limma_fdr")
parse_info<-merge(parse_info,integrate_limma_result%>%select(gene,cancer_type,limma_logFC,limma_fdr),by=c("gene","cancer_type"))

load("~/TCGA/kallisto_gencode_18_resource/TCGA/tumor/DE_united/filter_version/all_gene/integrate_DE_info_result_update.rds")
PSG_up<-integrate_DE_info_result%>%filter(cancer_type_positive>=3,cancer_type_positive>3*cancer_type_negative,phylo_age=="Primate")%>%.$gene
parse_info<-parse_info%>%mutate(score=ifelse((limma_logFC>=0.4&limma_fdr<=0.05)&(num>=3),2,ifelse((limma_logFC>=0.4&limma_fdr<=0.05)|(num>=3),1,0)))

select_gene<-unique(parse_info%>%filter(gene%in%PSG_up,score==2)%>%.$gene_symbol)
cat(select_gene,file="select_gene.txt",sep="\n")
select_gene_order<-rev(parse_info%>%filter(gene%in%PSG_up,score==2)%>%group_by(gene_symbol)%>%summarize(num=n())%>%arrange(num)%>%.$gene_symbol)
select_data<-parse_info%>%filter(gene_symbol%in%select_gene)

select_data_summary<-data.frame(select_data%>%group_by(gene_symbol)%>%summarize(num=sum(score==2)))
select_data_summary$gene_symbol<-factor(select_data_summary$gene_symbol,levels=select_gene_order)

bar_plot<-ggplot()+theme_classic()+theme(axis.title=element_text(size=15),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=15),axis.line=element_line(colour="grey70"))+geom_bar(data=select_data_summary,aes(x=gene_symbol,y=num),stat="identity",width=.7,fill="palegreen3")+scale_y_continuous(position="left",expand=c(0,0),limits=c(0,9.5),breaks=c(0,2,4,6,8))+labs(x="",y="Cancer type\nwith 2 evidences")
#  ggplot()+theme_classic()+theme(axis.title=element_text(size=15),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.text.x=element_text(size=15),axis.line=element_line(colour="grey70"))+geom_bar(data=select_data_summary,aes(x=gene_symbol,y=num),stat="identity",width=.7,fill="palegreen3")+coord_flip()+scale_y_continuous(position="right",expand=c(0,0),limits=c(0,9.5),breaks=c(0,4,8))+labs(x="",y="Cancer type\nwith 2 evidences")

select_data$gene_symbol<-factor(select_data$gene_symbol,levels=select_gene_order)

select_data$score<-ifelse(select_data$score==0,1,ifelse(select_data$score==1,2,5))
select_data$cancer_type<-factor(select_data$cancer_type,levels=rev(unique(select_data$cancer_type)))


dot_plot<-ggplot()+theme(axis.ticks=element_blank(),legend.key=element_blank(),legend.text=element_text(size=15),legend.title=element_text(size=15),axis.text.x=element_text(size=11,angle=90,hjust=.5,vjust=0.5),axis.text.y=element_text(size=15),panel.grid.major=element_line(colour="palegreen3",linetype="dotted",size=.7),panel.border=element_rect(colour="palegreen3",fill=NA),panel.background=element_blank(),legend.position="top")+geom_point(data=select_data,aes(y=cancer_type,x=gene_symbol,size=score),colour="violetred1",fill="violetred1")+scale_size_area(breaks=c(1,2,5),label=c(0,1,2),name="Evidence level")+labs(x="",y="")

ggsave(dot_plot,file="TCGA_depmap_dotplot.pdf",height=4.2,width=6.2)








