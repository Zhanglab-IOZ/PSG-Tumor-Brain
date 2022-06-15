library(ggplot2)
library(dplyr)

num=c(236,4,11,666)
total=c(1682,2633,138,18195)
Pro<-data.frame(num=num,total=total,type=c("UC","EM","PSG","All"))
Pro[,3]<-factor(Pro[,3],levels=c("All","UC","EM","PSG"))
ggplot(data=Pro,aes(x=type,y=num/total))+geom_bar(stat="identity",width=0.8,alpha=.6)+xlab("type")+ylab("Proportion")+theme_classic()+theme(legend.position="top",legend.title=element_blank(),legend.text=element_text(size=15),axis.ticks.x=element_blank(),axis.text.x=element_text(size=15,colour="black"),axis.text.y=element_text(size=15,colour="black"),axis.title.x=element_blank(),axis.title.y=element_text(size=15,colour="black"))+ylim(c(0,0.17))
ggsave("666propor1.pdf",width=4,height=4)