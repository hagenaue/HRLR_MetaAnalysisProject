#This is the code that actually worked for reannotating Sarah's P14 Affy 230 data:


#This is where I obtained the updated custom .cdf for defining the probesets:
http://nmg-r.bioinformatics.nl/NuGO_R.html


# "Annotation libraries for the MBNI Custom CDF annotation files (mirrored from MBNI; VERSION 19.0.0 Updated November 2015 / Bioconductor 3.2. Here, also the CDF and PROBE libraries are available for download.To obtain full microarray support, you need to download and install all three libraries for each chip!"

#Here's more instructions about using custom .cdf:
#http://brainarray.mbni.med.umich.edu/brainarray/database/customcdf/cdfreadme.htm

#Ah -  I was originally using the wrong set of annotation. Now fixed!

library(org.Rn.eg.db)
library(plyr)

install.packages(pkgs = c("rae230arnentrezg.db", "rae230arnentrezgcdf", "rae230arnentrezgprobe"), repos = "http://nmg-r.bioinformatics.nl/bioc")

library(rae230arnentrezg.db)
library(rae230arnentrezgcdf)
library(rae230arnentrezgprobe)

# Read in the CEL files in the directory, then normalize the data
#trying the simple version first -setting the working directory to the folder with Sarah's old .cel files (downloaded from NCBI GEO):

data2<-ReadAffy(cdfname ="rae230arnentrezg")
eset2 <- rma(data2)
write.exprs(eset2,file="data_customCDF.txt")
RMAExpression_customCDF<-read.delim("data_customCDF.txt", sep="\t")
str(RMAExpression_customCDF)
#'data.frame':	10019 obs. of  13 variables:

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
#[1] 9960
dim(RMAExpression_customCDFAnnotation2)
#[1] 10019     3

#So almost all of the results have gene symbols now. Yes...

write.csv(RMAExpression_customCDFAnnotation2, "RMAExpression_customCDFAnnotation2.csv")

2^5
[1] 32

#Ok, that worked, but after reading more about it, I believe that the HC data should be run separately from the NACC data, so I separated out the HC .cel files and re-ran the code.


#############I ran this code later to check for outliers##################

tempExpression<-read.csv("RMAExpression_customCDF.csv", header=T)
str(tempExpression)
# $ ProbesetID                  : Factor w/ 10019 levels "100036582_at",..: 1 2 3 4 5 6 7 8 9 10 ...
# $ X1_RN230_HC_P14_F15_H01.CEL : num  3.99 4.79 5 4.7 3.84 ...
# $ X10_RN230_HC_P14_F15_L05.CEL: num  3.99 4.61 4.79 4.78 3.99 ...
# $ X11_RN230_HC_P14_F15_H06.CEL: num  3.82 4.63 4.8 4.54 4.09 ...
# $ X12_RN230_HC_P14_F15_L06.CEL: num  3.59 4.53 4.8 4.62 4.37 ...
# $ X2_RN230_HC_P14_F15_L01.CEL : num  3.78 4.77 4.91 4.72 4.03 ...
# $ X3_RN230_HC_P14_F15_H02.CEL : num  3.88 4.53 5.05 4.79 4.27 ...
# $ X4_RN230_HC_P14_F15_L02.CEL : num  3.87 4.49 4.65 4.8 4.16 ...
# $ X5_RN230_HC_P14_F15_H03.CEL : num  3.98 4.73 4.9 4.67 4.27 ...
# $ X6_RN230_HC_P14_F15_L03.CEL : num  4.96 4.9 3.97 5.75 5.05 ...
# $ X7_RN230_HC_P14_F15_H04.CEL : num  3.74 4.46 5.02 4.84 3.96 ...
# $ X8_RN230_HC_P14_F15_L04.CEL : num  4.05 4.37 4.73 4.87 4.05 ...
# $ X9_RN230_HC_P14_F15_H05.CEL : num  3.9 4.43 4.86 4.85 3.97 ...

tempExpressionMatrix<-as.matrix(tempExpression[,-1])
str(tempExpressionMatrix)
# num [1:10019, 1:12] 3.99 4.79 5 4.7 3.84 ...
# - attr(*, "dimnames")=List of 2
# ..$ : NULL
# ..$ : chr [1:12] "X1_RN230_HC_P14_F15_H01.CEL" "X10_RN230_HC_P14_F15_L05.CEL" "X11_RN230_HC_P14_F15_H06.CEL" "X12_RN230_HC_P14_F15_L06.CEL" ...

write.csv(cor(tempExpressionMatrix), "SamplevsSample_CorMatrix.csv")
