
#tissue considered
tissue_select<-c("adrenal","appendix","bladder","bonemarrow","brain","breast","colon","duodenum","endometrium","esophagus","gallbladder","heart","kidney","liver","lung","lymphnode","ovary","pancreas","placenta","prostate","salivarygland","skin","smallintestine","spleen","stomach","testis","thyroid")
cutoff<-log((1+1),base=2)


get_mean_max<-function(x){
  cancer_type<-x
  tmp_file<-paste(cancer_type,"_exp_result",sep="")	
  tmp_data<-read.table(tmp_file,stringsAsFactors=F,header=T,row.name=1)
  tmp_mean<-apply(tmp_data,1,mean)
  tmp_max<-apply(tmp_data,1,max)
  tmp_median<-apply(tmp_data,1,median)
  tmp_data$mean<-tmp_mean
  tmp_data$median<-tmp_median
  tmp_data$max<-tmp_max
  return(tmp_data)
}

cbind_all<-function(x){
  obj.list <- lapply(x, get)
  names(obj.list) <- x
  integrate<-do.call(cbind, obj.list)
  return(integrate)
}

get_binary<-function(x){
  x<-as.numeric(x)
  result<-c()
  max_value<-max(x)
  min_value<-min(x)
  gap<-(max_value-min_value)/9
  for(i in 1:length(x)){
    if(max_value>0){
      if(x[i]>=min_value & x[i]<min_value+1*gap){
        result<-append(result,0)
      }else if(x[i]>=min_value+1*gap& x[i]<min_value+2*gap){
        result<-append(result,1)
      }else if(x[i]>=min_value+2*gap& x[i]<min_value+3*gap){
        result<-append(result,2)
      }else if (x[i]>=min_value+3*gap& x[i]<min_value+4*gap){
        result<-append(result,3)
      }else if (x[i]>=min_value+4*gap& x[i]<min_value+5*gap){
        result<-append(result,4)
      }else if (x[i]>=min_value+5*gap& x[i]<min_value+6*gap){
        result<-append(result,5)
      }else if (x[i]>=min_value+6*gap& x[i]<min_value+7*gap){
        result<-append(result,6)
      }else if (x[i]>=min_value+7*gap& x[i]<min_value+8*gap){
        result<-append(result,7)
      }else if (x[i]>=min_value+8*gap& x[i]<min_value+9*gap){
        result<-append(result,8)
      }else {
        result<-append(result,9)
      }
    }else{
      result<-rep(0,length(x))
    }
  }
  return(result)
}


get_tau<-function(x){
  max_value=max(x)
  if(max_value > 0){
    tau<-sum(1-x/max_value)/(length(x)-1)
    #index<-which.max(x)
    #final<-list(max_value,tau,index)
    final<-list(tau)
    return(final)
  }else{
    tau<-0
    #index<-sample(1:length(x),1)
    max_value<-0
    #final<-list(max_value,tau,index)
    final<-list(tau)
    return(final)
  }
}


get_expressed_tissue_info<-function(x,y){
  tmp_max_data<-as.numeric(x)
  tmp_max_data<-round(tmp_max_data,4)
  tmp_mean_data<-as.numeric(y)
  tmp_mean_data<-round(tmp_mean_data,4)
  tmp_tissue<-names(max_value)
  
  tmp_mean_index<-which(tmp_mean_data>=cutoff)
  tmp_mean_tissue<-names(max_value)[tmp_mean_index]
  tmp_mean_num<-length(tmp_mean_tissue)
  
  tmp_max_index<-which(tmp_max_data>=cutoff)
  tmp_max_tissue<-names(max_value)[tmp_max_index]
  tmp_max_num<-length(tmp_max_tissue)
  
  
  tmp_max_value<-tmp_max_data
  tmp_mean_value<-tmp_mean_data
  tmp_tissue<-tmp_tissue[order(tmp_mean_value,decreasing=T)]
  tmp_mean_tissue<-tmp_tissue[sort(match(tmp_mean_tissue,tmp_tissue))]
  tmp_max_tissue<-tmp_tissue[sort(match(tmp_max_tissue,tmp_tissue))]
  tmp_max_value<-tmp_max_value[order(tmp_mean_value,decreasing=T)]
  tmp_mean_value<-tmp_mean_value[order(tmp_mean_value,decreasing=T)]
  tmp_info<-paste(paste(tmp_tissue,tmp_max_value,tmp_mean_value,sep=":"),collapse=";")
  tmp_mean_tissue<-paste(tmp_mean_tissue,collapse=";")
  tmp_max_tissue<-paste(tmp_max_tissue,collapse=";")
  return(list(tmp_mean_num,tmp_mean_tissue,tmp_max_num,tmp_max_tissue,tmp_info))
}

#get the mean/median/max expression level for each tissue
for(i in tissue_select){tmp_exp<-get_mean_max(i);tmp_output<-paste("./",i,"_log_integrate_mean_max_pc_result",sep="");write.table(tmp_exp,file=tmp_output,quote=F,sep="\t",row.names=T,col.names=T);}


## use median value to calculate tau value; 
## use max value to determine the tissue with the highest expression level;
median_value<-read.table("integrate_median.txt",stringsAsFactors=F,header=T,row.names=1)
median_value<-median_value[select_tissue]

max_value<-read.table("integrate_max.txt",stringsAsFactors=F,header=T,row.names=1)
max_value<-max_value[select_tissue]

median_value_max<-as.numeric(apply(median_value,1,max))

tau_tissue<-names(median_value)[apply(median_value,1,which.max)]

max_value_max<-c();for(i in 1:length(tau_tissue)){tmp_tissue<-tau_tissue[i];tmp_index<-which(names(max_value)==tmp_tissue);max_value_max[i]<-as.numeric(max_value[i,][tmp_index]);rm(i,tmp_tissue,tmp_index)}

## get the binary
median_binary<-as.data.frame(t(apply(median_value,1,get_binary)))
names(median_binary)<-names(median_value)
gene_name<-c();tau_value<-c();
for(i in 1:nrow(median_binary)){gene_name[i]<-rownames(median_binary)[i];final<-get_tau(as.numeric(median_binary[i,]));tau_value[i]<-final[[1]];rm(final,i)}

gene_name<-as.data.frame(gene_name)
tau_value<-as.data.frame(tau_value)
tau_tissue<-as.data.frame(tau_tissue)
names(tau_tissue)<-"tau_tissue"

tau_result<-cbind_all(c("gene_name","tau_value","tau_tissue"))
tau_result<-as.data.frame(tau_result,stringsAsFactors=F)
tau_result$tau_tissue_mean_exp<-median_value_max
tau_result$tau_tissue_max_exp<-max_value_max
tau_result$tau_tissue_mean_exp<-as.numeric(tau_result$tau_tissue_mean_exp)
tau_result$tau_tissue_max_exp<-as.numeric(tau_result$tau_tissue_max_exp)
tau_result$tau_value<-as.numeric(tau_result$tau_value)
tau_result$tau_value<-round(tau_result$tau_value,4)
tau_result$tau_tissue_mean_exp_result<-ifelse(tau_result$tau_tissue_mean_exp>=cutoff,"yes","no")
tau_result$tau_tissue_max_exp_result<-ifelse(tau_result$tau_tissue_max_exp>=cutoff,"yes","no")
tau_result$tau_tissue_mean_exp<-round(tau_result$tau_tissue_mean_exp,4)
tau_result$tau_tissue_max_exp<-round(tau_result$tau_tissue_max_exp,4)


tissue_gene_name<-c();tissue_mean_num<-c();tissue_mean<-c();tissue_max_num<-c();tissue_max<-c();tissue_info<-c();for(i in 1:nrow(median_value)){tissue_gene_name[i]<-rownames(median_value)[i];tmp_max<-as.numeric(max_value[i,]);tmp_mean<-as.numeric(median_value[i,]);final<-get_expressed_tissue_info(tmp_max,tmp_mean);tissue_mean_num[i]<-final[[1]];tissue_mean[i]<-final[[2]];tissue_max_num[i]<-final[[3]];tissue_max[i]<-final[[4]];tissue_info[i]<-final[[5]];rm(final,i,tmp_max,tmp_mean)}

tissue_gene_name<-as.data.frame(tissue_gene_name)
tissue_mean<-as.data.frame(tissue_mean)
tissue_max<-as.data.frame(tissue_max)
tissue_mean_num<-as.data.frame(tissue_mean_num)
tissue_max_num<-as.data.frame(tissue_max_num)
tissue_info<-as.data.frame(tissue_info)


tissue_exp_info<-cbind_all(c("tissue_gene_name","tissue_mean_num","tissue_max_num"))

names(tissue_exp_info)[1]<-"gene_name"

tau_result<-merge(tau_result,tissue_exp_info,by="gene_name")

tissue_exp_info<-cbind_all(c("tissue_gene_name","tissue_mean","tissue_mean_num","tissue_max","tissue_max_num","tissue_info"))
names(tissue_exp_info)[1]<-"gene_name"

## save the tau result and expression information for each gene; 
write.table(tau_result,file="./tau_result.txt",quote=F,sep="\t",row.names=F,col.names=T)
write.table(tissue_exp_info,file="./tissue_exp_info.txt",quote=F,sep="\t",row.names=F,col.names=T)



