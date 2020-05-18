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

out_file='/Users/macbook/Desktop/APOP_PAN_CANCER/top_ten/topten.csv'
write.table(cbind("Cancer","HR","p-value","wald-p","logrank-p","C","%95 CI lower","%95 CI upper","min(PI)","max(PI)","mean(PI)","median(PI)"),
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
  b=beta[beta$Gene==genes[i],][2][,1]
  PI=PI+b*select(data1, genes[i])[,1]
}

surv_object <- Surv(time = data1$OS.time/30, event = data1$vital_status)

#survival analysis: fits cox ph model to find HR for mean cut
fit1 <- survfit(surv_object~PI>median(PI)); 
#summary(fit1);
#ggsurvplot(fit1, data=test1)


fit1.coxph <- coxph(surv_object~PI>median(PI))
# summary(fit1.coxph);
first <- coef(summary(fit1.coxph))
#first
a=summary(fit1.coxph)
#fit1.coxph$concordance[6]
gg<-ggsurvplot(fit1, data=data1,legend=c(0.8,0.8),xlab = "Time (Months)",axes.offset=FALSE,legend.labs = c("Low Risk", "High Risk"))

ggsave('/Users/macbook/Desktop/APOP_PAN_CANCER/top_ten/KMplots/'+cancer+'.jpg', plot = print(gg),height = 5 ,width = 7, dpi = 1000)


{write.table(cbind(cancer,first[2],first[5],a$waldtest[3],a$logtest[3],a$concordance[1],a$conf.int[3],a$conf.int[4],min(PI),max(PI),mean(PI),median(PI)),
             file=out_file,row.names=F,col.names=F,sep = ',',append = T);#output file
}

}
#write.table(cbind("Cancer","HR","p-value","wald-p","logrank-p","C","%95 CI lower","%95 CI upper"),
#            file=out_file,row.names=F,col.names=F,sep = ',')


## -1	-0.9	-0.8	0.3	-0.9	-0.7	-0.8 HR=3.36,p=10e-12
# 1,1,1,0,0,0,0  HR=4.13,p=10e-11 ;359,90 (cutoff=1)
#0,1,1,1,0,0,1	HR=4.004254941	p=5.78E-11	;362,87 (cutoff=0)
#0,1,1,1,0,0,1	HR=6.400338	p=2.488283e-15	;390,59 (cutoff=1)

ggsurv<-ggsurvplot(
  fit1,                            # survfit object with calculated statistics.
  data = data1,                   # data used to fit survival curves.
  risk.table = TRUE,               # show risk table.
  #title=" \t \t BCL2+BCLXL-BAX-BAK ",
  conf.int = FALSE,                 # show confidence intervals for
  # palette = c("#E7B800", "#2E9FDF"),
  xlim = c(0,405),                  # present narrower X axis, but not affect
  #ylim=c(0.1,1),
  xlab = "Time in months",         # customize X axis label.
  # break.time.by = 20,             # break X axis in time intervals by 500.
  ggtheme = theme_light(),         # customize plot and risk table with a theme.
  risk.table.height = 0.25,         # the height of the risk table
  risk.table.y.text = FALSE,        # show bars instead of names in text annotations
  # # in legend of risk table.
  conf.int.style = "step",         # customize style of confidence intervals
  legend=c(0.8,0.8),
  legend.labs = c("High Risk", "Low Risk"),  # change legend labels.
  # ncensor.plot = TRUE,               # plot the number of censored subjects at time t
  # ncensor.plot.height = 0.25,
  # pval = TRUE,                       # show p-value of log-rank test.
  axes.offset=FALSE,
)
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate("label", x = 300, y = 0.50, # x and y coordinates of the text
                    label = "HR = 6.40, p-val= 2.49e-15 \n 95%CI(4.04 - 10.14)
Wald p= 2e-15 & logrank, p= <2e-16", size = 4)
ggsurv
