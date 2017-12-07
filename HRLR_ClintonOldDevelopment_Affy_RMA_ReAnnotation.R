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




