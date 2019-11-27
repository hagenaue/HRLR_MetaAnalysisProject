#Bmp4 over development qRT-PCR Analysis
##Megan Hagenauer, Data from Kathryn Hilde and Alex Stefanov (using their own tissue and Pam Patterson's tissue)
##2-22-2019


#R-Studio Version 1.0.153 – © 2009-2017 RStudio, Inc.
#Mozilla/5.0 (Macintosh; Intel Mac OS X 10_12_6) AppleWebKit/605.1.15 (KHTML, like Gecko)

#Analysis Notes:

#The calibration curves were analyzed separately in Excel and revealed efficiencies close to 1 (i.e., each unit of Cq represents approximately a doubling of the initial concentration) for both Bmp4 and reference houskeeping gene (GAPDH) within concentrations between 0.1-2 uL. Therefore, a simple Livak method calculation should be appropriate for examining group differences in Bmp4.

#The samples for this analysis (and the calibration analysis) were originally run in triplicate but were already averaged in Excel before loading into R.


setwd("~/Documents/Microarray Gen/HRLR/Bmp4PCR")


P14_Data<-read.csv("P14_ForR.csv", header=T, stringsAsFactors = T)
str(P14_Data)
# 'data.frame':	12 obs. of  4 variables:
#   $ Well            : int  1 2 3 4 5 6 7 8 9 10 ...
# $ TG.Cq.Mean      : num  27.8 26.5 25.8 26.5 26.2 ...
# $ HKG.Cq.Mean     : num  18.1 18.7 18.4 18.4 18.6 ...
# $ Group.Assignment: Factor w/ 2 levels "HR","LR": 1 1 1 1 1 1 2 2 2 2 ...

table(P14_Data[,4])
# HR LR 
# 6  6 

mean(P14_Data[,2])
#[1] 24.35944
mean(P14_Data[,3])
#[1] 18.585


P90_Data<-read.csv("P90_ForR.csv", header=T, stringsAsFactors = T)
str(P90_Data)
# 'data.frame':	12 obs. of  4 variables:
#   $ Well            : Factor w/ 12 levels "Bmp4 HR 1","Bmp4 HR 2",..: 1 2 3 4 5 6 7 8 9 10 ...
# $ TG.Cq.Mean      : num  25.6 26.2 25.9 25.2 25.7 ...
# $ HKG.Cq.Mean     : num  16.4 16.7 16.5 16.5 15.8 ...
# $ Group.Assignment: Factor w/ 2 levels "HR","LR": 1 1 1 1 1 1 2 2 2 2 ...

table(P90_Data[,4])
# HR LR 
# 6  6 
mean(P90_Data[,2])
#[1] 24.29083
mean(P90_Data[,3])
#[1] 16.56583


######################


#Double-checking the housekeeping gene:

t.test(P14_Data[,3]~P14_Data[,4])

# Welch Two Sample t-test
# 
# data:  P14_Data[, 3] by P14_Data[, 4]
# t = -0.23227, df = 8.4661, p-value = 0.8218
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.7463256  0.6085478
# sample estimates:
#   mean in group HR mean in group LR 
# 18.55056         18.61944 


#Calculating Delta Cq (Target-Housekeeping):
P14_DeltaCq<-P14_Data[,2]-P14_Data[,3]

t.test(P14_DeltaCq~P14_Data[,4])

# Welch Two Sample t-test
# 
# data:  P14_DeltaCq by P14_Data[, 4]
# t = 6.0961, df = 5.5923, p-value = 0.00115
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   2.211218 5.266560
# sample estimates:
#   mean in group HR mean in group LR 
# 7.643889         3.905000 

#Plotting:
#Most people would prefer I plot as -DeltaDeltaCT.  
P14_DeltaDeltaCq=(P14_DeltaCq-mean(P14_DeltaCq[which(P14_Data[,4]=="LR")]))*-1

#-DeltaDeltaCq:
mean(P14_DeltaDeltaCq[which(P14_Data[,4]=="HR")])
#[1] -3.738889

pdf("Boxplot_P14_Bmp4PCR.pdf", width=3, height=5)
boxplot(P14_DeltaDeltaCq~P14_Data[,4], xlab="Phenotype", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(3,2), main="P14", ylim=c(-6,2))
stripchart(P14_DeltaDeltaCq~P14_Data[,4], vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1.7, col = 'black')
mtext(expression(paste("Bmp4 Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

####################


######


#Double-checking the housekeeping gene:

t.test(P90_Data[,3]~P90_Data[,4])

# Welch Two Sample t-test
# 
# data:  P90_Data[, 3] by P90_Data[, 4]
# t = -1.8317, df = 8.0979, p-value = 0.1039
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -1.074256  0.122034
# sample estimates:
#   mean in group HR mean in group LR 
# 16.32778         16.80389 


#Calculating Delta Cq (Target-Housekeeping):
P90_DeltaCq<-P90_Data[,2]-P90_Data[,3]

t.test(P90_DeltaCq~P90_Data[,4])

# Welch Two Sample t-test
# 
# data:  P90_DeltaCq by P90_Data[, 4]
# t = 6.8696, df = 8.737, p-value = 8.439e-05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   1.762928 3.505961
# sample estimates:
#   mean in group HR mean in group LR 
# 9.042222         6.407778 


#Plotting:
#Most people would prefer I plot as -DeltaDeltaCT.  
P90_DeltaDeltaCq=(P90_DeltaCq-mean(P90_DeltaCq[which(P90_Data[,4]=="LR")]))*-1

#-DeltaDeltaCq:
mean(P90_DeltaDeltaCq[which(P90_Data[,4]=="HR")])
#[1] -2.634444

pdf("Boxplot_P90_Bmp4PCR.pdf", width=3, height=5)
boxplot(P90_DeltaDeltaCq~P90_Data[,4], xlab="Phenotype", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(3,2), main="P90", ylim=c(-6,2))
stripchart(P90_DeltaDeltaCq~P90_Data[,4], vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1.7, col = 'black')
mtext(expression(paste("Bmp4 Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()
