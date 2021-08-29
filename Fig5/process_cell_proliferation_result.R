library(dplyr)
library(reshape2)
library(ggplot2)

get_info<-function(input_file){
  time<-strsplit(input_file,"\\.|_")[[1]][1]
  cell_line<-strsplit(input_file,"\\.|_")[[1]][2]
  data<-read.table(input_file,stringsAsFactors=F,h=T,sep="\t")
  names(data)<-c("batch","Elapsed","Blank","CTL","NC","DDX11","AS1","AS2","AC","Blank_se","CTL_se","NC_se","DDX11_se","AS1_se","AS2_se","AC_se")
  data$time<-time
  data$cell_line<-cell_line
  data$Elapsed<-round(data$Elapsed)
  data<-data.frame(data%>%filter(Elapsed<=72))
  data<-data%>%filter(Elapsed%%6==0)
  ## CTL_blank_pvalue<-as.numeric(t.test(data$Blank,data$CTL,paired=T)$p.value)
  ## DDX11_blank_pvalue<-as.numeric(t.test(data$Blank,data$DDX11,paired=T)$p.value)
  ## DDX11_CTL_pvalue<-as.numeric(t.test(data$CTL,data$DDX11,paired=T)$p.value)
  data$Blank_scale<-data$Blank/data$Blank[1]
  data$CTL_scale<-data$CTL/data$CTL[1]
  data$NC_scale<-data$NC/data$NC[1]
  data$DDX11_scale<-data$DDX11/data$DDX11[1]
  CTL_blank_normalization_pvalue<-as.numeric(t.test(data$Blank_scale,data$CTL_scale,paired=T)$p.value)
  CTL_blank_pvalue<-as.numeric(t.test(data$Blank,data$CTL,paired=T)$p.value)
  DDX11_blank_normalization_pvalue<-as.numeric(t.test(data$Blank_scale,data$DDX11_scale,paired=T)$p.value)
  DDX11_blank_pvalue<-as.numeric(t.test(data$Blank,data$DDX11,paired=T)$p.value)
  DDX11_CTL_normalization_pvalue<-as.numeric(t.test(data$CTL_scale,data$DDX11_scale,paired=T)$p.value)
  DDX11_CTL_pvalue<-as.numeric(t.test(data$CTL,data$DDX11,paired=T)$p.value)
  #return(data)
  return(c(time,DDX11_blank_pvalue,DDX11_blank_normalization_pvalue,DDX11_CTL_pvalue,DDX11_CTL_normalization_pvalue,CTL_blank_pvalue,CTL_blank_normalization_pvalue))
}

files<-list.files(pattern="txt")
data<-data.frame(do.call(rbind,lapply(files,get_info)))
names(data)<-c("time","DDX11_blank","DDX11_blank_nor","DDX11_CTL","DDX11_CTL_nor","CTL_blank","CTL_blank_nor")

#data$Elapsed<-round(data$Elapsed)
#data<-data%>%filter(Elapsed<=72)
#data<-data%>%filter(Elapsed%%6==0)
#data$Elapsed<-factor(data$Elapsed)





get_plot<-function(input_file){
  time<-strsplit(input_file,"\\.|_")[[1]][1]
  cell_line<-strsplit(input_file,"\\.|_")[[1]][2]
  data<-read.table(input_file,stringsAsFactors=F,h=T,sep="\t")
  names(data)<-c("batch","Elapsed","Blank","CTL","NC","DDX11","AS1","AS2","AC","Blank_se","CTL_se","NC_se","DDX11_se","AS1_se","AS2_se","AC_se")
  data$time<-time
  data$cell_line<-cell_line
  data$Elapsed<-round(data$Elapsed)
  data<-data.frame(data%>%filter(Elapsed<=72))
  data<-data%>%filter(Elapsed%%6==0)
  data$Elapsed<-factor(data$Elapsed)

  point_data<-data%>%select(Elapsed,Blank,CTL,DDX11)
  point_data<-melt(point_data)

  bar_data<-data%>%select(Elapsed,Blank_se,CTL_se,DDX11_se)
  bar_data<-melt(bar_data)
  bar_data$variable<-factor(bar_data$variable,labels=c("Blank","CTL","DDX11"))
  names(bar_data)[3]<-"se_value"

  data<-merge(point_data,bar_data,by=c("Elapsed","variable"))
  data$variable<-factor(data$variable,labels=c("Untransfected","siRNA-Ctrl","siRNA-DDX11"))

  output_file<-paste(time,"_",cell_line,".pdf",sep="")
  #pdf(output_file)
  plot_result<-ggplot()+theme_classic()+theme(panel.background=element_blank(),axis.title=element_text(size=10),axis.text=element_text(size=10),legend.position=c(0.35,0.8),legend.title=element_blank(),legend.text=element_text(size=10))+geom_errorbar(data=data,aes(x=Elapsed,ymin=value-se_value,ymax=value+se_value),colour="black",width=.1)+geom_smooth(method="loess",data=data,aes(x=Elapsed,y=value,group=variable,colour=variable),se=F)+geom_point(data=data,aes(x=Elapsed,y=value,colour=variable,shape=variable),size=2.5)+scale_shape_manual(values=c(19,15,17))+labs(x="Time (hours)",y="Phase Object Confluence (%)")+scale_x_discrete(breaks=c(0,12,24,36,48,60,72))
  ggsave(plot_result,file=output_file,height=2.8,width=2.8)
  #dev.off()
}

lapply(files,get_plot)

