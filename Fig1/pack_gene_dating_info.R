library(dplyr)

setwd("~/TCGA/kallisto_gencode_18_resource/tissue_specificity/before_clean_version/SPM_scale/integrate_tau_SPM")

#load the tau results
load("tau_SPM_age_chr_info_result.rds")

#load the phylostratigraphy data
load("phylostratigraphy.rds")
names(phylo_age_info)<-c("gene","phylo_age")

tau_SPM_age_chr_info_hybrid_result<-merge(tau_SPM_age_chr_info_result,phylo_age_info,by="gene",all.x=T)

rm(tau_SPM_age_chr_info_result)

# split into age groups
tau_SPM_age_chr_info_hybrid_result$phylo_age_class<-NA
tau_SPM_age_chr_info_hybrid_result$phylo_age_class[which(tau_SPM_age_chr_info_hybrid_result$phylo_age==1&tau_SPM_age_chr_info_hybrid_result$age_info<5)]<-"Cellular_organisms"
tau_SPM_age_chr_info_hybrid_result$phylo_age_class[which(tau_SPM_age_chr_info_hybrid_result$phylo_age==2&tau_SPM_age_chr_info_hybrid_result$age_info<5)]<-"Eukaryota"
tau_SPM_age_chr_info_hybrid_result$phylo_age_class[which(tau_SPM_age_chr_info_hybrid_result$phylo_age==3&tau_SPM_age_chr_info_hybrid_result$age_info<5)]<-"Opisthokonta"
tau_SPM_age_chr_info_hybrid_result$phylo_age_class[which(tau_SPM_age_chr_info_hybrid_result$phylo_age==4&tau_SPM_age_chr_info_hybrid_result$age_info<5)]<-"Metazoa"
tau_SPM_age_chr_info_hybrid_result$phylo_age_class[which(tau_SPM_age_chr_info_hybrid_result$phylo_age==5&tau_SPM_age_chr_info_hybrid_result$age_info<5)]<-"Eumetazoa"
tau_SPM_age_chr_info_hybrid_result$phylo_age_class[which(tau_SPM_age_chr_info_hybrid_result$phylo_age==6&tau_SPM_age_chr_info_hybrid_result$age_info<5)]<-"Bilateria"
tau_SPM_age_chr_info_hybrid_result$phylo_age_class[which(tau_SPM_age_chr_info_hybrid_result$phylo_age==7&tau_SPM_age_chr_info_hybrid_result$age_info<5)]<-"Chordata"
tau_SPM_age_chr_info_hybrid_result$phylo_age_class[which(tau_SPM_age_chr_info_hybrid_result$phylo_age==8&tau_SPM_age_chr_info_hybrid_result$age_info<5)]<-"Euteleostomi"
tau_SPM_age_chr_info_hybrid_result$phylo_age_class[which(tau_SPM_age_chr_info_hybrid_result$phylo_age==9&tau_SPM_age_chr_info_hybrid_result$age_info<5)]<-"Ammiota"
tau_SPM_age_chr_info_hybrid_result$phylo_age_class[which(tau_SPM_age_chr_info_hybrid_result$phylo_age==10&tau_SPM_age_chr_info_hybrid_result$age_info<5)]<-"Mammalia"
tau_SPM_age_chr_info_hybrid_result$phylo_age_class[which(tau_SPM_age_chr_info_hybrid_result$phylo_age==11&tau_SPM_age_chr_info_hybrid_result$age_info<5)]<-"Theria"
tau_SPM_age_chr_info_hybrid_result$phylo_age_class[which(tau_SPM_age_chr_info_hybrid_result$age_info>=5&tau_SPM_age_chr_info_hybrid_result$age_info<=6)]<-"Eutheria"
tau_SPM_age_chr_info_hybrid_result$phylo_age_class[which(tau_SPM_age_chr_info_hybrid_result$age_info==7)]<-"Euarchontoglires"

tau_SPM_age_chr_info_hybrid_result$phylo_age_class[which(tau_SPM_age_chr_info_hybrid_result$age_info>=8)]<-"Primate"

tau_SPM_age_chr_info_hybrid_result$phylo_age_class[is.na(tau_SPM_age_chr_info_hybrid_result$phylo_age_class)]<-"no_info"

tau_SPM_age_chr_info_hybrid_result<-subset(tau_SPM_age_chr_info_hybrid_result,select=-c(phylo_age))
names(tau_SPM_age_chr_info_hybrid_result)[which(names(tau_SPM_age_chr_info_hybrid_result)=="phylo_age_class")]<-"phylo_age"

tau_SPM_age_chr_info_hybrid_result$phylo_age<-factor(tau_SPM_age_chr_info_hybrid_result$phylo_age,levels=c("Cellular_organisms","Eukaryota","Opisthokonta","Metazoa","Eumetazoa","Bilateria","Chordata","Euteleostomi","Ammiota","Mammalia","Theria","Eutheria","Euarchontoglires","Primate","no_info"))

save(tau_SPM_age_chr_info_hybrid_result,file="tau_SPM_age_chr_info_hybrid_result_update.rds")

