`+` <- function(e1, e2) {
  if (is.character(e1) | is.character(e2)) {
    paste0(e1, e2)
  } else {
    base::`+`(e1, e2)
  }
}

library(survival)
library(survminer)
library(dplyr)

cancer1=scan(file = '/Users/macbook/Desktop/APOP_PAN_CANCER/top_ten/genelists/cancers.csv', what = 'character', sep = ',')

out_file='/Users/macbook/Desktop/APOP_PAN_CANCER/top_ten/topten_cutp.csv'
write.table(cbind("Cancer","HR","p-value","wald-p","logrank-p","C","%95 CI lower","%95 CI upper","min(PI)","max(PI)","cutoff"),
            file=out_file,row.names=F,col.names=F,sep = ',')
for(k in seq(from=1, to=length(cancer1), by=1))
{





cancer=cancer1[k]
beta_file='/Users/macbook/Desktop/APOP_PAN_CANCER/HR_APOP/'+cancer+'_HRmedian_raw.csv'
tcga_datafile="/Users/macbook/Desktop/webserver-pan-can/data_basic/RDS_files/"+cancer+".RDS"

beta<-read.csv(file=beta_file,header =TRUE, sep = ",", dec = ".")
data1<-readRDS(tcga_datafile)
genes <- scan(file = '/Users/macbook/Desktop/APOP_PAN_CANCER/top_ten/genelists/'+cancer+'_genes.csv', what = 'character', sep = ',')

#select(LUSC_C, gep[i])[,1]
PI=0
for(i in seq(from=1, to=length(genes), by=1))#no of cancers
{
  k=beta[beta$Gene==genes[i],][2][,1]
  if (k>0)
  {
  b=1*(select(data1, genes[i])[,1]>median(select(data1, genes[i])[,1]))}
  
  if (k<0)
  {b=1*(select(data1, genes[i])[,1]<median(select(data1, genes[i])[,1]))}
  
  PI=PI+b
}

library("survMisc")
fit1<-Surv(time = data1$OS.time/30, event = data1$vital_status)
fit2 <- coxph(fit1~PI, data = data1)
cut <- cutp(fit2)$PI

g=cut[1,1][,1]

rm(fit1)

surv_object <- Surv(time = data1$OS.time/30, event = data1$vital_status)

#survival analysis: fits cox ph model to find HR for mean cut
fit1 <- survfit(surv_object~PI>=g$PI); 
#summary(fit1);
ggsurvplot(fit1, data=data1)


fit1.coxph <- coxph(surv_object~PI>=g$PI)
# summary(fit1.coxph);
first <- coef(summary(fit1.coxph))
#first
a=summary(fit1.coxph)
#a
 #max U and Q, min p value

fit1.coxph$concordance[6]
  gg<-ggsurvplot(fit1, data=data1,legend=c(0.8,0.8),xlab = "Time (Months)",axes.offset=FALSE,legend.labs = c("Low Risk", "High Risk"))

  ggsave('/Users/macbook/Desktop/APOP_PAN_CANCER/top_ten/KMplots/'+cancer+'cutp.jpg', plot = print(gg),height = 5 ,width = 7, dpi = 1000)


  {write.table(cbind(cancer,first[2],first[5],a$waldtest[3],a$logtest[3],a$concordance[1],a$conf.int[3],a$conf.int[4],min(PI),max(PI),g$PI),
               file=out_file,row.names=F,col.names=F,sep = ',',append = T);#output file
  }

}

