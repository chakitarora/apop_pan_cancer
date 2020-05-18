library("glmnet")
library("survival")


tcga_datafile="/Users/macbook/Desktop/webserver-pan-can/data_basic/RDS_files/BRCA.RDS"
patient.data <- readRDS(tcga_datafile)
patient.data <- subset(patient.data, OS.time > 0)


OS<-patient.data$OS.time
cens<-patient.data$vital_status
data<-patient.data[,c(16:ncol(patient.data))]

coxlasso.fit <- cv.glmnet(
  x = data.matrix(patient.data[,c(16:ncol(patient.data))]),
  y = Surv(patient.data$OS.time, patient.data$vital_status),
  family = 'cox', maxit=2000)


coxlasso <- glmnet(
  x = data.matrix(patient.data[,c(16:ncol(patient.data))]),
  y = Surv(patient.data$OS.time, patient.data$vital_status),
  family = 'cox', maxit=2000)

plot(coxlasso.fit)

Coefficients <- coef(coxlasso, s = coxlasso.fit$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
