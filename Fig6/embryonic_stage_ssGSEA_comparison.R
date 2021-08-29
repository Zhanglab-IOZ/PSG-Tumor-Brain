library(dplyr)
library(ggplot2)


load("cell_cycle_all_three.rds")

ssgsea_result<-ssgsea_result%>%mutate(organ_class=ifelse(organ=="forebrain","Cerebrum",ifelse(organ=="hindbrain","Cerebellum","Other organs")))

ssgsea_result$organ_class<-factor(ssgsea_result$organ_class,levels=c("Cerebrum","Cerebellum","Other organs"))


pdf("stage_1_cell_cycle_ssgsea.pdf",height=4,width=4)
ggplot(data=ssgsea_result%>%filter(stage=="P1"),aes(x=organ_class,y=ssgsea_value))+theme_classic()+theme(legend.text=element_text(size=14),legend.title=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=14),axis.title=element_text(size=14),legend.position="bottom")+geom_violin(aes(fill=organ_class),colour=NA,alpha=.7)+geom_boxplot(width=.2,fill="black",outlier.colour=NA,size=.3)+stat_summary(fun.y=median,geom="point",fill="white",shape=21,size=3)+labs(x="",y="Enrichment score")+scale_fill_viridis_d()
dev.off()


wilcox.test(ssgsea_result%>%filter(stage=="P1",organ_class=="Cerebrum")%>%.$ssgsea_value,ssgsea_result%>%filter(stage=="P1",organ_class=="Cerebellum")%>%.$ssgsea_value)$p.value

wilcox.test(ssgsea_result%>%filter(stage=="P1",organ_class=="Cerebellum")%>%.$ssgsea_value,ssgsea_result%>%filter(stage=="P1",organ_class=="Other organs")%>%.$ssgsea_value)$p.value

