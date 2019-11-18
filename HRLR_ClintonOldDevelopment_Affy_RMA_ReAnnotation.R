#This is the code that actually worked for reannotating Sarah's developmental data:

#Note: after reading more about this, I decided that the Hippocampal (HC) data should be run separately from the NACC data, so I separated out the HC .cel files and re-ran the code.

#This is where I obtained the updated custom .cdf for defining the probesets:
http://nmg-r.bioinformatics.nl/NuGO_R.html


# "Annotation libraries for the MBNI Custom CDF annotation files (mirrored from MBNI; VERSION 19.0.0 Updated November 2015 / Bioconductor 3.2. Here, also the CDF and PROBE libraries are available for download.To obtain full microarray support, you need to download and install all three libraries for each chip!"

#Here's more instructions about using custom .cdf:
#http://brainarray.mbni.med.umich.edu/brainarray/database/customcdf/cdfreadme.htm

#Ah -  I was originally using the wrong set of annotation. Now fixed!

library(org.Rn.eg.db)
library(plyr)

install.packages(pkgs = c("rgu34arnentrezg.db", "rgu34arnentrezgcdf", "rgu34arnentrezgprobe"), repos = "http://nmg-r.bioinformatics.nl/bioc")

library(rgu34arnentrezg.db)
library(rgu34arnentrezgcdf)
library(rgu34arnentrezgprobe)

# Read in the CEL files in the directory, then normalize the data
#trying the simple version first -setting the working directory to the folder with Sarah's old .cel files (downloaded from NCBI GEO):

data2<-ReadAffy(cdfname ="rgu34arnentrezg")
eset2 <- rma(data2)
write.exprs(eset2,file="data_customCDF.txt")
RMAExpression_customCDF<-read.delim("data_customCDF.txt", sep="\t")
str(RMAExpression_customCDF)
#'data.frame':	4649 obs. of  72 variables:

write.csv(RMAExpression_customCDF, "RMAExpression_customCDF.csv")


#Alright, now I need the additional annotation:

head(RMAExpression_customCDF)
RMAExpression_EntrezID<-sub("_at", "", RMAExpression_customCDF[,1])
head(RMAExpression_EntrezID)
RMAExpression_customCDFAnnotation<-data.frame(RMAExpression_customCDF[,1], RMAExpression_EntrezID, stringsAsFactors = F )
colnames(RMAExpression_customCDFAnnotation)<-c("ProbesetID", "EntrezGeneID")



x <- org.Rn.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

xx[1]
[1] "Slc39a4l"
xx$"100360501"
[1] "Rnh1"

GeneSymbol<-unlist(xx, use.names=FALSE)
EntrezGeneID<-rep(names(xx), lengths(xx))
table(lengths(xx))
# 1 
# 42306
#A 1:1 mapping. Nice...

EntrezVsGeneSymbol<-data.frame(EntrezGeneID, GeneSymbol, stringsAsFactors=F)

RMAExpression_customCDFAnnotation2<-join(RMAExpression_customCDFAnnotation, EntrezVsGeneSymbol, by="EntrezGeneID", type="left")

sum(is.na(RMAExpression_customCDFAnnotation2[,3])==F)
#[1] 4588
dim(RMAExpression_customCDFAnnotation2)
#[1] 4649    3

#So almost all of the results have gene symbols now. Yes...

write.csv(RMAExpression_customCDFAnnotation2, "RMAExpression_customCDFAnnotation2.csv")


#*************************

#Making a couple of extra plots (starting with an empty workspace)

setwd("~/Documents/Microarray Gen/HRLR/HRLR_Studies/Sarah Old Development/RMAJustHC")

RMAExpression<-read.csv("RMAExpression_customCDF.csv", header=T, stringsAsFactors = F)
str(RMAExpression)
# $ ProbesetID                     : chr  "100034253_at" "100036582_at" "100125372_at" "100125373_at" ...
# $ GSM731605_P7_HR_HPC_H01.CEL.gz : num  4.4 4.16 3.23 5.41 7.98 ...
# $ GSM731606_P7_HR_HPC_H02.CEL.gz : num  4.78 4.25 3.18 5.66 7.96 ...
# $ GSM731607_P7_HR_HPC_H03.CEL.gz : num  5.18 4.07 3.2 5.96 8.12 ...

colnames(RMAExpression)
# [1] "ProbesetID"                      "GSM731605_P7_HR_HPC_H01.CEL.gz"  "GSM731606_P7_HR_HPC_H02.CEL.gz"  "GSM731607_P7_HR_HPC_H03.CEL.gz" 
# [5] "GSM731608_P7_HR_HPC_H04.CEL.gz"  "GSM731609_P7_HR_HPC_H05.CEL.gz"  "GSM731610_P7_HR_HPC_H06.CEL.gz"  "GSM731611_P14_HR_HPC_H01.CEL.gz"
# [9] "GSM731612_P14_HR_HPC_H02.CEL.gz" "GSM731613_P14_HR_HPC_H03.CEL.gz" "GSM731614_P14_HR_HPC_H04.CEL.gz" "GSM731615_P14_HR_HPC_H05.CEL.gz"
# [13] "GSM731616_P14_HR_HPC_H06.CEL.gz" "GSM731617_P21_HR_HPC_H01.CEL.gz" "GSM731618_P21_HR_HPC_H02.CEL.gz" "GSM731619_P21_HR_HPC_H03.CEL.gz"
# [17] "GSM731620_P21_HR_HPC_H04.CEL.gz" "GSM731621_P21_HR_HPC_H05.CEL.gz" "GSM731622_P21_HR_HPC_H06.CEL.gz" "GSM731623_P7_LR_HPC_L01.CEL.gz" 
# [21] "GSM731624_P7_LR_HPC_L02.CEL.gz"  "GSM731625_P7_LR_HPC_L03.CEL.gz"  "GSM731626_P7_LR_HPC_L04.CEL.gz"  "GSM731627_P7_LR_HPC_L05.CEL.gz" 
# [25] "GSM731628_P7_LR_HPC_L06.CEL.gz"  "GSM731629_P14_LR_HPC_L01.CEL.gz" "GSM731630_P14_LR_HPC_L02.CEL.gz" "GSM731631_P14_LR_HPC_L03.CEL.gz"
# [29] "GSM731632_P14_LR_HPC_L04.CEL.gz" "GSM731633_P14_LR_HPC_L05.CEL.gz" "GSM731634_P14_LR_HPC_L06.CEL.gz" "GSM731635_P21_LR_HPC_L01.CEL.gz"
# [33] "GSM731636_P21_LR_HPC_L02.CEL.gz" "GSM731637_P21_LR_HPC_L03.CEL.gz" "GSM731638_P21_LR_HPC_L04.CEL.gz" "GSM731639_P21_LR_HPC_L05.CEL.gz"
# [37] "GSM731640_P21_LR_HPC_L06.CEL.gz"

#Notes from QC say no need to remove outliers.

#For plotting:
HRvsLR<-as.factor(c(rep("HR", 18), rep("LR", 18)))
HRvsLR<-relevel(HRvsLR, ref="LR")

Age<-as.factor(c(rep("P7", 6), rep("P14", 6), rep("P21", 6), rep("P7", 6), rep("P14", 6), rep("P21", 6)))
levels(Age)
Age<-relevel(Age, ref="P7")

RMA_Annotation<-read.csv("RMAExpression_customCDFAnnotation2.csv", header=T, stringsAsFactors = F)
str(RMA_Annotation)
# 'data.frame':	4649 obs. of  3 variables:
#   $ ProbesetID  : chr  "100034253_at" "100036582_at" "100125372_at" "100125373_at" ...
# $ EntrezGeneID: chr  "100034253" "100036582" "100125372" "100125373" ...
# $ GeneSymbol  : chr  "Gnl3l" "Olr1867" "Ces1f" "Pnrc2" ...


pdf("Neurocan_VsAgeHRLR.pdf", width=4, height=4)
boxplot(as.numeric(RMAExpression[which(RMA_Annotation[,3]=="Ncan"),-1])~Age*HRvsLR, col=c(2,2,2, 3,3,3), ylab="Ncan: Log(2) Signal", main="MBNI_AffymetrixRgU34A_F6", las=2, cex.lab=1.3)
#stripplot(RMAExpression[which(RMA_Annotation[,3]=="Ncan"),-1])~Age*HRvsLR, jitter.data=T)
dev.off()

pdf("Bmp4_VsAgeHRLR.pdf", width=4, height=4)
boxplot(as.numeric(RMAExpression[which(RMA_Annotation[,3]=="Bmp4"),-1])~Age*HRvsLR, col=c(2,2,2, 3,3,3), ylab="Bmp4: Log(2) Signal", main="MBNI_AffymetrixRgU34A_F6", las=2, cex.lab=1.3)
dev.off()

#Can't plot Sox9, Sox2, or Mki67 (not in dataset)

pdf("Cav1_VsAgeHRLR.pdf", width=4, height=4)
boxplot(as.numeric(RMAExpression[which(RMA_Annotation[,3]=="Cav1"),-1])~Age*HRvsLR, col=c(2,2,2, 3,3,3), ylab="Cav1: Log(2) Signal", main="MBNI_AffymetrixRgU34A_F6", las=2, cex.lab=1.3)
dev.off()

pdf("Hes5_VsAgeHRLR.pdf", width=4, height=4)
boxplot(as.numeric(RMAExpression[which(RMA_Annotation[,3]=="Hes5"),-1])~Age*HRvsLR, col=c(2,2,2, 3,3,3), ylab="Hes5: Log(2) Signal", main="MBNI_AffymetrixRgU34A_F6", las=2, cex.lab=1.3)
dev.off()


#For comparison: Sarah's new colonydevelopment study

setwd("~/Documents/Microarray Gen/HRLR/HRLR_Studies/New Sarah Development Data/Ibirt Server")

NC_P7<-read.csv("NC_annotatedlogdata_P7.csv", header=T, stringsAsFactors = F)
str(NC_P7)
#' #'data.frame':	10 obs. of  10677 variables:
#' $ X              : chr  "HR_HPC_P07_rep1.txt" "HR_HPC_P07_rep2.txt" "HR_HPC_P07_rep3.txt" "HR_HPC_P07_rep4.txt" ...
#' $ phenotype      : chr  "HR" "HR" "HR" "HR" ...
#' $ age            : chr  "P7" "P7" "P7" "P7" ...
#' $ A1bg           : num  6.6 6.13 6.42 5.71 5.96 ...
#' $ A1cf           : num  5.58 6.08 6.15 5.78 5.62 ...

NC_P14<-read.csv("NC_annotatedlogdata_P14.csv", header=T, stringsAsFactors = F)
str(NC_P14)

NC_P21<-read.csv("NC_annotatedlogdata_P21.csv", header=T, stringsAsFactors = F)
str(NC_P21)

NC_Adult<-read.csv("NC_annotatedlogdata_Adult.csv", header=T, stringsAsFactors = F)
str(NC_Adult)

NC_Expression<-cbind(t(NC_P7[,-c(1:3)]), t(NC_P14[,-c(1:3)]), t(NC_P21[,-c(1:3)]), t(NC_Adult[,-c(1:3)]))
str(NC_Expression)
colnames(NC_Expression)<-c(NC_P7$X, NC_P14$X, NC_P21$X, NC_Adult$X)
HRvsLR<-as.factor(c(NC_P7$phenotype, NC_P14$phenotype, NC_P21$phenotype, NC_Adult$phenotype))
levels(HRvsLR)
HRvsLR<-relevel(HRvsLR, ref="LR")
length(HRvsLR)
Age<-as.factor(c(NC_P7$age, NC_P14$age, NC_P21$age, NC_Adult$age))
length(Age)
levels(Age)
Age<-reorder(Age, c(rep(1, 10), rep(2, 10), rep(3, 9), rep(4, 10)))
relevel(Age, ref="P7")

pdf("Neurocan_VsAgeHRLR.pdf", width=5, height=4)
boxplot(as.numeric(NC_Expression[which(row.names(NC_Expression)=="Ncan"),])~Age*HRvsLR, col=c(2,2,2,2, 3,3,3,3), ylab="Ncan: Log(2) Signal", main="Alabama_NimbleGen_F34", las=2, cex.lab=1.3)
#stripplot(RMAExpression[which(RMA_Annotation[,3]=="Ncan"),-1])~Age*HRvsLR, jitter.data=T)
dev.off()

pdf("Bmp4_VsAgeHRLR.pdf", width=5, height=4)
boxplot(as.numeric(NC_Expression[which(row.names(NC_Expression)=="Bmp4"),])~Age*HRvsLR, col=c(2,2,2,2, 3,3,3,3), ylab="Bmp4: Log(2) Signal", main="Alabama_NimbleGen_F34", las=2, cex.lab=1.3)
dev.off()

pdf("Cav1_VsAgeHRLR.pdf", width=5, height=4)
boxplot(as.numeric(NC_Expression[which(row.names(NC_Expression)=="Cav1"),])~Age*HRvsLR, col=c(2,2,2,2, 3,3,3,3), ylab="Cav1: Log(2) Signal", main="Alabama_NimbleGen_F34", las=2, cex.lab=1.3)
dev.off()

pdf("Hes5_VsAgeHRLR.pdf", width=5, height=4)
boxplot(as.numeric(NC_Expression[which(row.names(NC_Expression)=="Hes5"),])~Age*HRvsLR, col=c(2,2,2,2, 3,3,3,3), ylab="Hes5: Log(2) Signal", main="Alabama_NimbleGen_F34", las=2, cex.lab=1.3)
dev.off()

#No Sox9, no Sox2, no Mki67 - too bad!

pdf("Mki67_VsAgeHRLR.pdf", width=5, height=4)
boxplot(as.numeric(NC_Expression[which(row.names(NC_Expression)=="Mki67"),])~Age*HRvsLR, col=c(2,2,2,2, 3,3,3,3), ylab="Mki67: Log(2) Signal", main="Alabama_NimbleGen_F34", las=2, cex.lab=1.3)
dev.off()

#I think the boxplots maybe making the results look more variable than they are (because they keep treating individual points as outliers - silly! It is a sample size of 5!)  I'm going to see if violin plots look slightly more truthful.

library(ggplot2)

DataForPlotting<-cbind.data.frame(as.numeric(NC_Expression[which(row.names(NC_Expression)=="Ncan"),]), Age, HRvsLR)
colnames(DataForPlotting)<-"Ncan"

p<-ggplot(DataForPlotting, aes(x=Age, y=Ncan, fill=HRvsLR)) + geom_violin(trim=FALSE)+
 
  #geom_boxplot(width=0.1)+
  labs(title="Alabama_NimbleGen_F34", x="Age", y ="Ncan: Log(2) Signal")
  #p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
p + geom_dotplot(binaxis='y', stackdir='center',
                 position=position_dodge(1), dotsize=0.5)

#For some reason, adding dots and summary stats seem to not overlay properly on the violins.  I'm not sure that I have the energy to troubleshoot this right now, I'm just going to stick with boxplots for a while.
