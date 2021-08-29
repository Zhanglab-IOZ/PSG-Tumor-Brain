Args<-commandArgs()

library(limma)
library(edgeR)
library(dplyr)
library(ggplot2)

cancer_type<-Args[6]


load(paste("./",cancer_type,"/normal_count_result.rds",sep=""))
load(paste("./",cancer_type,"/tumor_count_result.rds",sep=""))
stopifnot(all.equal(tumor_final$gene,normal_expression$gene))

load("tau_SPM_age_chr_info_result.rds")

tumor_result<-tumor_final[-1]
names(tumor_result)<-paste(names(tumor_result),"_tumor",sep="")
normal_result<-normal_expression[-1]
names(normal_result)<-paste(names(normal_result),"_normal",sep="")
integrate_result<-cbind(tumor_result,normal_result)



sample_info<-data.frame(names(integrate_result))
names(sample_info)<-"sample_name"
sample_info$type<-"normal"
sample_info$type[grep("tumor$",sample_info$sample_name)]<-"tumor"
sample_info$type<-factor(sample_info$type,levels=c("normal","tumor"))
design <- model.matrix(~ type, data = sample_info)



dge <- DGEList(counts=integrate_result)
dge<-calcNormFactors(dge,method="upperquartile")
pdf(paste(cancer_type,"_voom.pdf",sep=""))
v <- voom(dge, design, plot=TRUE)
dev.off()
voom_result<-v$E
save(voom_result,file="voom_result.rds")

fit <- lmFit(v, design)
fit <- eBayes(fit)
limma_result<-topTable(fit,n=Inf,sort="none",coef=2)
limma_result$gene<-rownames(limma_result)
limma_age_result<-merge(limma_result,tau_SPM_age_chr_info_result,by="gene")
save(limma_age_result,file="limma_age_result.rds")
