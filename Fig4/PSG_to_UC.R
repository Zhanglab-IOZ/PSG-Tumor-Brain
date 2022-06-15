library(dplyr)
library(ggplot2)

data<-data.frame(all_num=c(272,34,11),signal_num=c(101,23,9))
data<-data %>% mutate(value=signal_num/all_num)
ata$class<-c("All PSGs","Cell cycle corrlation","Cell cycle 666")
data$class<-factor(data$class,levels = c("Cell cycle 666","Cell cycle corrlation","All PSGs"))

ggplot()+theme_classic()+
  theme(axis.title=element_text(size=15),axis.text.x=element_text(size=14),legend.position="none",axis.ticks.y=element_blank(),axis.text.y =element_blank())+
  geom_bar(data=data,aes(x=class,y=value,fill=class),stat="identity",width=0.75)+coord_flip()+scale_y_continuous(position = "right",expand = c(0,0),breaks=c(0,0.4,0.8),limits=c(0,0.84),labels = c(0,0.4,0.8))+
  labs(x="Class",y="Proportion of UC derived genes")
