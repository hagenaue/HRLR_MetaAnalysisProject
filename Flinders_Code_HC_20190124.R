#Code for re-analyzing the Flinders Sensitive vs. Flinders Resistant Hippocampal microarray data from Blaveri et al. 2010
#Megan Hagenauer, 1/24/2019
#Affymetrix GeneChip Rat Genome 230 2.0


setwd("~/Documents/Microarray Gen/Flinders_Blaveri2010")

library(GEOquery)
gse <- getGEO("GSE20388", GSEMatrix = FALSE)
#lots of warnings. I wonder why. They all look like this:
#Warning messages:
#1: In readLines(con, n = chunksize) :
  #seek on a gzfile connection returned an internal error

head(Meta(gse))
str(Meta(gse))
Meta(GSMList(gse)$GSM510523)$characteristics_ch1
# [1] "strain: FRL"                  "brain region: Frontal Cortex" "cohort: Cohort 1"   

sub("strain: ","", Meta(GSMList(gse)$GSM510523)$characteristics_ch1[1])
#[1] "FRL"
sub("brain region: ","", Meta(GSMList(gse)$GSM510523)$characteristics_ch1[2])
#[1] "Frontal Cortex"
sub("cohort: ","", Meta(GSMList(gse)$GSM510523)$characteristics_ch1[3])
#[1] "Cohort 1"

SampleID<-as.matrix(names(GSMList(gse)))
length(GSMList(gse))
#[1] 78
Strain<-matrix("a", nrow=78, ncol=1)
Tissue<-matrix("a", nrow=78, ncol=1)
Cohort<-matrix("a", nrow=78, ncol=1)


i<-1
for(i in c(1:78)){
  Strain[i]<-sub("strain: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[1])
  Tissue[i]<-sub("brain region: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[2])
  Cohort[i]<-sub("cohort: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[3])
}

Strain<-as.factor(Strain)
#Has FRL as reference

Cohort<-as.factor(Cohort)
#Has Cohort 1 as reference

SampleCharacteristics<-data.frame(SampleID, Strain, Tissue, Cohort, stringsAsFactors=F)
str(SampleCharacteristics)

# 'data.frame':	78 obs. of  4 variables:
#   $ SampleID: chr  "GSM510523" "GSM510524" "GSM510526" "GSM510527" ...
# $ Strain  : Factor w/ 2 levels "FRL","FSL": 1 1 1 1 1 1 1 1 1 1 ...
# $ Tissue  : chr  "Frontal Cortex" "Frontal Cortex" "Frontal Cortex" "Frontal Cortex" ...
# $ Cohort  : Factor w/ 2 levels "Cohort 1","Cohort 2": 1 1 1 1 1 1 1 1 1 1 ...

#Looks good.

write.csv(SampleCharacteristics, "Flinders_SampleCharacteristics.csv")

table(SampleCharacteristics$Tissue)

#for the HC-Only Analyses:
SampleCharacteristics_HC<-SampleCharacteristics[Tissue=="Hippocampus",]
write.csv(SampleCharacteristics_HC, "SampleCharacteristics_HC.csv")

SampleCharacteristics<-SampleCharacteristics_HC
SampleID<-SampleCharacteristics_HC$SampleID
Strain<-SampleCharacteristics_HC$Strain
Tissue<-SampleCharacteristics_HC$Tissue
Cohort<-SampleCharacteristics_HC$Cohort 

Tissue
# [1] "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus"
# [9] "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus"
# [17] "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus"
# [25] "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus"
# [33] "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus" "Hippocampus"

#Perfect

rm(gse)


#Now let's re-run RMA on their data:

library(org.Hs.eg.db)
library(plyr)
library(affy)

#This is where I obtained the updated custom .cdf for defining the probesets:
http://nmg-r.bioinformatics.nl/NuGO_R.html

#cdf and the chip: Rattus norvegicus rat2302

#hmm... those are pretty old. Let's see if we can get something newer:

source("http://www.bioconductor.org/biocLite.R")
biocLite(c("rat2302rnentrezg.db", "rat2302rnentrezgcdf", "rat2302rnentrezgprobe"))
#Error: failed to update BiocInstaller:
#namespace ‘BiocInstaller’ is imported by ‘affy’ so cannot be unloaded

library(org.Rn.eg.db)

install.packages(pkgs = c("rat2302rnentrezg.db", "rat2302rnentrezgcdf", "rat2302rnentrezgprobe"), repos = "http://nmg-r.bioinformatics.nl/bioc")

#Hmm... these are still relatively old annotation:
# rat2302rnentrezg.db_19.0.2.tar.gz
# rat2302rnentrezgcdf_19.0.0.tar.gz
# rat2302rnentrezgprobe_19.0.0.tar.gz
#On the website, they say that this "Current Custom CDF is version 19.0.0 (Nov 33, 2014)."

#This seems to be the location of the newer files:
http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/24.0.0/entrezg.asp
#I downloaded the appropriate files by hand.

#Scrap installation code that didn't work:
# install.packages(pkgs = c("rat2302rnentrezg.db", "rat2302rnentrezgcdf", "rat2302rnentrezgprobe"), repos = "http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/24.0.0/entrezg.asp")
# # Warning in install.packages :
# #   unable to access index for repository http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/24.0.0/entrezg.asp/src/contrib:
# #   Line starting '<html> ...' is malformed!
# #   Warning in install.packages :
# #   packages ‘rat2302rnentrezg.db’, ‘rat2302rnentrezgcdf’, ‘rat2302rnentrezgprobe’ are not available (for R version 3.4.1)
# # Warning in install.packages :
# #   unable to access index for repository http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/24.0.0/entrezg.asp/bin/macosx/el-capitan/contrib/3.4:
# #   Line starting '<html> ...' is malformed!
# 
# 
# #Using the gui:
# install.packages("~/Documents/Microarray Gen/Flinders_Blaveri2010/CEL/rat2302_24/rat2302rnentrezg.db_24.0.0.tar.gz", repos = NULL, type = "source")
# 
# 
# # Warning in install.packages :
# #   installation of package ‘/Users/mhh/Documents/Microarray Gen/Flinders_Blaveri2010/CEL/rat2302_24/rat2302rnentrezg.db_24.0.0.tar.gz’ had non-zero exit status
# 
# install.packages(pkgs = c("rat2302rnentrezg.db", "rat2302rnentrezgcdf", "rat2302rnentrezgprobe"), repos ="~/Documents/Microarray Gen/Flinders_Blaveri2010/CEL/rat2302_24/" )
# 
# # Warning in install.packages :
# #   unable to access index for repository ~/Documents/Microarray Gen/Flinders_Blaveri2010/CEL/rat2302_24/src/contrib:
# #   cannot open URL '~/Documents/Microarray Gen/Flinders_Blaveri2010/CEL/rat2302_24/src/contrib/PACKAGES'
# # Warning in install.packages :
# #   packages ‘rat2302rnentrezg.db’, ‘rat2302rnentrezgcdf’, ‘rat2302rnentrezgprobe’ are not available (for R version 3.4.1)
# # Warning in install.packages :
# #   unable to access index for repository ~/Documents/Microarray Gen/Flinders_Blaveri2010/CEL/rat2302_24/bin/macosx/el-capitan/contrib/3.4:
# #   cannot open URL '~/Documents/Microarray Gen/Flinders_Blaveri2010/CEL/rat2302_24/bin/macosx/el-capitan/contrib/3.4/PACKAGES'
# 
# install.packages("C:Users/mhh/Documents/Microarray Gen/Flinders_Blaveri2010/CEL/rat2302_24/rat2302rnentrezg.db_24.0.0.tar.gz", repos = NULL, type ="source")
# 
# # Warning in install.packages :
# #   installation of package ‘C:Users/mhh/Documents/Microarray Gen/Flinders_Blaveri2010/CEL/rat2302_24/rat2302rnentrezg.db_24.0.0.tar.gz’ had non-zero exit status
# 
# install.packages("C:Users/mhh/Documents/Microarray Gen/Flinders_Blaveri2010/CEL/rat2302_24/Rat2302_Rn_ENTREZG_24.0.0.zip", repos = NULL, type ="win.binary")
# #Error in install.packages : cannot install Windows binary packages on this platform
# 
# 
# install.packages("C:Users/mhh/Documents/Microarray Gen/Flinders_Blaveri2010/CEL/rat2302_24/Rat2302_Rn_ENTREZG_24.0.0/", 
#                  repos = NULL, 
#                  type = "source")
# # Warning in install.packages :
# #   installation of package ‘C:Users/mhh/Documents/Microarray Gen/Flinders_Blaveri2010/CEL/rat2302_24/Rat2302_Rn_ENTREZG_24.0.0/’ had non-zero exit status
# 
# #I also tried installing via the terminal:
# Megans-MBP:~ mhh$ R CMD INSTALL -l  "C:Users/mhh/Documents/Microarray Gen/Flinders_Blaveri2010/CEL/rat2302_24/rat2302rnentrezg.db_24.0.0.tar.gz"
# Error: ERROR: no packages specified
# 
# R CMD INSTALL "C:Users/mhh/Documents/Microarray Gen/Flinders_Blaveri2010/CEL/rat2302_24/rat2302rnentrezg.db_24.0.0.tar.gz"
# Error: ERROR: no packages specified
# 
# install.packages("http://mbni.org/customcdf/24.0.0/entrezg.download/rat2302rnentrezg.db_24.0.0.tar.gz")
# # Warning in install.packages :
# #   package ‘http://mbni.org/customcdf/24.0.0/entrezg.download/rat2302rnentrezg.db_24.0.0.tar.gz’ is not available (for R version 3.4.1)

#After a long series of troubleshooting, it turned out that earlier versions of the package (v22 and earlier) would install fine, but not v23 or v24 (???)
#Here's the code that finally worked:
install.packages("http://mbni.org/customcdf/22.0.0/entrezg.download/rat2302rnentrezg.db_22.0.0.tar.gz", type="source", repos=NULL)
install.packages("http://mbni.org/customcdf/22.0.0/entrezg.download/rat2302rnentrezgcdf_22.0.0.tar.gz", type="source", repos=NULL)
install.packages("http://mbni.org/customcdf/22.0.0/entrezg.download/rat2302rnentrezgprobe_22.0.0.tar.gz", type="source", repos=NULL)

library(rat2302rnentrezg.db)
library(rat2302rnentrezgcdf)
library(rat2302rnentrezgprobe)


#Changed working directory to where the .cel files are located
setwd("~/Documents/Microarray Gen/Flinders_Blaveri2010/CEL/GSE20388_RAW")


data2<-ReadAffy(cdfname ="rat2302rnentrezg")
str(data2)
data2
# AffyBatch object
# size of arrays=834x834 features (41 kb)
# cdf=rat2302rnentrezg (14132 affyids)
# number of samples=39
# number of genes=14132
# annotation=rat2302rnentrezg
# notes=

eset2 <- rma(data2)

setwd("~/Documents/Microarray Gen/Flinders_Blaveri2010")
write.exprs(eset2,file="data_customCDF.txt")
RMAExpression_customCDF<-read.delim("data_customCDF.txt", sep="\t", stringsAsFactors = F)
str(RMAExpression_customCDF)
#'data.frame':	14132 obs. of  40 variables:
write.csv(RMAExpression_customCDF, "RMAExpression_customCDF.csv")


ScanDate<-protocolData(data2)$ScanDate
#Yep, there are definitely different scan dates here. 

library(reshape2)
ScanDate_Split<-colsplit(ScanDate, pattern=" ", c("ScanDate", "ScanTime"))
table(ScanDate_Split$ScanDate)
# 02/23/06 02/24/06 02/28/06 06/30/06 07/07/06 
# 1       10       10        9        9 

#So several potential batches. 2/23-2/28 are close enough together that they may be the same batch.

table(ScanDate_Split$ScanDate, Strain, Cohort)

#Looks like the batches are partially redundant with cohort, and almost evenly balanced by line (hurray!)

# , , Cohort = Cohort 1
# 
# Strain
#           FRL FSL
# 02/23/06   0   1
# 02/24/06   5   5
# 02/28/06   7   3
# 06/30/06   0   0
# 07/07/06   0   0
# 
# , , Cohort = Cohort 2
# 
# Strain
#           FRL FSL
# 02/23/06   0   0
# 02/24/06   0   0
# 02/28/06   0   0
# 06/30/06   5   4
# 07/07/06   5   4



#####################################################

#Adding annotation:

head(RMAExpression_customCDF)
RMAExpression_EntrezID<-sub("_at", "", RMAExpression_customCDF[,1])
head(RMAExpression_EntrezID)
RMAExpression_customCDFAnnotation<-data.frame(RMAExpression_customCDF[,1], RMAExpression_EntrezID, stringsAsFactors = F )
colnames(RMAExpression_customCDFAnnotation)<-c("ProbesetID", "EntrezGeneID")

library(org.Rn.eg.db)
x <- org.Rn.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

GeneSymbol<-unlist(xx, use.names=FALSE)
EntrezGeneID<-rep(names(xx), lengths(xx))
table(lengths(xx))
# 1 
# 47850 

EntrezVsGeneSymbol<-data.frame(EntrezGeneID, GeneSymbol, stringsAsFactors=F)

RMAExpression_customCDFAnnotation2<-join(RMAExpression_customCDFAnnotation, EntrezVsGeneSymbol, by="EntrezGeneID", type="left")
head(RMAExpression_customCDFAnnotation2)

sum(is.na(RMAExpression_customCDFAnnotation2[,3])==F)
#[1] 14074
dim(RMAExpression_customCDFAnnotation2)
#[1] 14132     3
#So almost all of the results have gene symbols.

write.csv(RMAExpression_customCDFAnnotation2, "RMAExpression_customCDFAnnotation2.csv")


SignalSortedNoNA3<-as.matrix(RMAExpression_customCDF[,-1])
#I know this is a weird variable name - I just set things up to recycle some code

#Double-checking that the Sample info and Expression data are in the same order:
cbind(SampleID, colnames(SignalSortedNoNA3))
# Yep, they are both in ascending numeric order.


SampleCharacteristics<-data.frame(SampleCharacteristics, ScanDate_Split)


####Quality Control################

setwd("~/Documents/Microarray Gen/Flinders_Blaveri2010/QC")

png("Boxplot_RMAExpression_customCDF.png", width=2000, height=400)
boxplot(SignalSortedNoNA3)
dev.off()
#Beautiful - almost looks quantile normalized it is so regular. 


#Let's check PCA and see if there are any other obvious problem samples:

################################
# #Run principal components analysis (PCA) to determine which major gradients of sample-sample correlations exist in the data (i.e. who is similar to whom):

pcaNormFilterednoOutliers<-prcomp(t(SignalSortedNoNA3))
tmp<-pcaNormFilterednoOutliers$x[,1:4]
write.table(tmp, "PCA_1_4.txt", sep="\t")


PCeigenvectors<-pcaNormFilterednoOutliers$rotation[ ,c(1:4)]
PCeigenvectors2<-cbind(PCeigenvectors, RMAExpression_customCDFAnnotation2)
write.csv(PCeigenvectors2, "PCeigenvectors.csv")

PC1noOutliers<-pcaNormFilterednoOutliers$x[,1]
PC2noOutliers<-pcaNormFilterednoOutliers$x[,2]

PC3noOutliers<-pcaNormFilterednoOutliers$x[,3]
PC4noOutliers<-pcaNormFilterednoOutliers$x[,4]

# #Output a scree plot for the PCA (no outliers):
png("10 PCA Scree Plot1.png")
plot(summary(pcaNormFilterednoOutliers)$importance[2,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off()

png("10 PCA Scree Plot2.png")
plot(summary(pcaNormFilterednoOutliers)$importance[3,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Cumulative Proportion of Variance Explained", col=3)
dev.off()


# #Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
png("10 PC1 vs PC2_byStrain.png")
plot(PC1noOutliers~PC2noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers", col=as.factor(SampleCharacteristics$Strain))
dev.off()
#Strain isn't strongly related to PC1 & PC2
#There don't appear to be any  major outliers.

png("10 PC1 vs PC2_byCohort.png")
plot(PC1noOutliers~PC2noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers", col=as.factor(SampleCharacteristics$Cohort))
dev.off()
#Cohort absolutely drives PC1 and PC2

png("10 PC1 vs PC2_byScanDate.png")
plot(PC1noOutliers~PC2noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers", col=as.factor(SampleCharacteristics$ScanDate))
dev.off()
#The scandates within the cohorts don't seem to matter much.


# #Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
png("10 PC3 vs PC4_byStrain.png")
plot(PC3noOutliers~PC4noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers", col=as.factor(SampleCharacteristics$Strain))
dev.off()
#Strain is related to PC3 and PC4

png("10 PC3 vs PC4_byCohort.png")
plot(PC3noOutliers~PC4noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers", col=as.factor(SampleCharacteristics$Cohort))
dev.off()
#Cohort is related to PC3... kind-of

png("10 PC3 vs PC4_byScanDate.png")
plot(PC3noOutliers~PC4noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers", col=as.factor(SampleCharacteristics$ScanDate))
dev.off()
#The scandates within the cohorts don't seem to matter much.
#Phew, that means we can probably throw that variable out.


SubjectPCA<-cbind(PC1noOutliers, PC2noOutliers, PC3noOutliers, PC4noOutliers)

PCAoutput<-data.frame(SampleCharacteristics, SubjectPCA)
write.csv(PCAoutput, "PCAoutput.csv")


#######################################

#Visualize the sample-sample correlations using a heatmap:
png("09 Sample Sample Correlations Heatmap.png")
image(cor(SignalSortedNoNA3), main="Visualizing the correlations between entire samples (by index#)", xlab="Red=Less correlated, Light yellow=Highly correlated")
dev.off()
#Note that the heatmap can be tailored to focus on a certain level of correlation by using the command zlim=c(lower limit, upper limit)

#Visualize the sample-sample correlations using a boxplot:
png("09 Boxplot Sample Sample Correlations.png", width=1000, height=600)
boxplot(data.frame(cor(SignalSortedNoNA3)), cex=0.25, las=3, par(cex.axis=0.75, mar=c(15,4,4,4)), main="Boxplot of sample-sample correlations", xlab="Subject", ylab="Sample-Sample Correlations")
Median10thQuantile<-median(apply((cor(SignalSortedNoNA3)), 1, quantile, 0.1))
MedianQuantile<-median(apply((cor(SignalSortedNoNA3)), 1, quantile, 0.5))
abline(a=Median10thQuantile, b=0, col=2)
abline(a=MedianQuantile, b=0, col=3)
mtext(paste("Median Sample-Sample Correlation=", round(MedianQuantile, digits=3), sep=" ")) 
dev.off()
#Beautiful data - no outliers.  


########################################################################

#Model Output:
setwd("~/Documents/Microarray Gen/Flinders_Blaveri2010/ModelOutput")

GeneBySubjVar2_Pvalues<-matrix(0, length(SignalSortedNoNA3[,1]), 3)
GeneBySubjVar2_Betas<-matrix(0, length(SignalSortedNoNA3[,1]), 3)
GeneBySubjVar2_Tstat<-matrix(0, length(SignalSortedNoNA3[,1]), 3)
colnames(GeneBySubjVar2_Pvalues)<-c("Intercept",  "Strain", "Cohort")
colnames(GeneBySubjVar2_Betas)<-c("Intercept",  "Strain", "Cohort")
colnames(GeneBySubjVar2_Tstat)<-c("Intercept",  "Strain", "Cohort")
row.names(GeneBySubjVar2_Pvalues)<-row.names(SignalSortedNoNA3)
row.names(GeneBySubjVar2_Betas)<-row.names(SignalSortedNoNA3)
row.names(GeneBySubjVar2_Tstat)<-row.names(SignalSortedNoNA3)

head(GeneBySubjVar2_Pvalues)


for(i in c(1:length(SignalSortedNoNA3[,1]))){
  
  temp<-summary.lm(lm(SignalSortedNoNA3[i,]~Strain+Cohort))
  
  GeneBySubjVar2_Betas[i,]<-temp$coefficients[,1]
  GeneBySubjVar2_Tstat[i,]<-temp$coefficients[,3]
  GeneBySubjVar2_Pvalues[i,]<-temp$coefficients[,4]
  
}

GeneBySubjVar2_Pvalues2<-cbind(RMAExpression_customCDFAnnotation2, GeneBySubjVar2_Pvalues)
GeneBySubjVar2_Betas2<-cbind(RMAExpression_customCDFAnnotation2, GeneBySubjVar2_Betas)

write.csv(GeneBySubjVar2_Pvalues2, "GeneBySubjVar2_Pvalues.csv")
write.csv(GeneBySubjVar2_Betas2, "GeneBySubjVar2_Betas.csv")

GeneBySubjVar2_Tstat2<-cbind(RMAExpression_customCDFAnnotation2, GeneBySubjVar2_Tstat)

write.csv(GeneBySubjVar2_Tstat2, "GeneBySubjVar2_Tstat.csv")


for (i in c(1:length(GeneBySubjVar2_Pvalues[1,]))){
  png(paste(paste("17 Histogram of Raw Pvalues for", colnames(GeneBySubjVar2_Pvalues)[i], sep="  "), "png", sep="."))	
  hist(GeneBySubjVar2_Pvalues[,i], breaks=100, col=i, main=paste("Raw P-values for", colnames(GeneBySubjVar2_Pvalues)[i], sep="  "), xlab="Unadjusted p-value", ylab="Count")
  abline(a=(length(GeneBySubjVar2_Pvalues[,1])/100), b=0)
  dev.off()		
}	


GeneBySubjVar2_PvaluesAdj<-matrix(0, length(GeneBySubjVar2_Pvalues[,1]), length(GeneBySubjVar2_Pvalues[1,]))
colnames(GeneBySubjVar2_PvaluesAdj)<-colnames(GeneBySubjVar2_Pvalues)
row.names(GeneBySubjVar2_PvaluesAdj)<-row.names(GeneBySubjVar2_Pvalues)


library(multtest)

for (i in c(1:length(GeneBySubjVar2_Pvalues[1,]))){
  
  #Applying multiple-comparison corrections to the raw p-values (Benjamini-Hochberg):
  TempPvalAdj<-mt.rawp2adjp(GeneBySubjVar2_Pvalues[,i], proc=c("BH"))
  GeneBySubjVar2_PvaluesAdj[,i]<-TempPvalAdj$adjp[order(TempPvalAdj$index),2]
  
}

GeneBySubjVar2_PvaluesAdj2<-cbind(RMAExpression_customCDFAnnotation2, GeneBySubjVar2_PvaluesAdj)
write.csv(GeneBySubjVar2_PvaluesAdj2, "GeneBySubjVar2_PvaluesAdj.csv")

GeneBySubjVar2DF<-as.data.frame(cbind(GeneBySubjVar2_Betas, GeneBySubjVar2_Tstat, GeneBySubjVar2_Pvalues, GeneBySubjVar2_PvaluesAdj))

temp<-cbind(RMAExpression_customCDFAnnotation2, GeneBySubjVar2DF)
write.csv(temp, "GeneBySubjVar2DF.csv" )

#Looking at the results, they resemble the original analysis (e.g., FSL have low Col6a1 and Tmem176a)

###########################


#####################
#Out of utter curiousity...

#Running Cell Type analyses:
setwd("~/Documents/Microarray Gen/Flinders_Blaveri2010/CellType")


install.packages("devtools")

library("devtools")

install_github("hagenaue/BrainInABlender")


library("BrainInABlender")

temp<-data.frame(as.character(RMAExpression_customCDFAnnotation2[,3]), SignalSortedNoNA3, stringsAsFactors=F)
write.csv(temp, "UserInput.csv")
str(temp)
table(table(RMAExpression_customCDFAnnotation2[,3]))
SirUnMixALotOutput<-Sir_UnMixALot(userInput=temp, dataColumns=c(2:40), geneColumn=1, species="mouse")

PublicationSpecific_CellTypeIndex<-SirUnMixALotOutput$PublicationSpecific_CellTypeIndex
AveragePrimary_CellTypeIndex<-SirUnMixALotOutput$AveragePrimary_CellTypeIndex

str(PublicationSpecific_CellTypeIndex)

#Wow- the output is really beautiful. Beautiful clusters.




################################################################
#Recycling some code to examine cell type vs. line:

CellTypeBySubjVar2_Pvalues<-matrix(0, length(AveragePrimary_CellTypeIndex[,1]), 3)
CellTypeBySubjVar2_Betas<-matrix(0, length(AveragePrimary_CellTypeIndex[,1]), 3)
CellTypeBySubjVar2_Tstat<-matrix(0, length(AveragePrimary_CellTypeIndex[,1]), 3)
colnames(CellTypeBySubjVar2_Pvalues)<-c("Intercept",  "Strain", "Cohort")
colnames(CellTypeBySubjVar2_Betas)<-c("Intercept",  "Strain", "Cohort")
colnames(CellTypeBySubjVar2_Tstat)<-c("Intercept",  "Strain", "Cohort")
row.names(CellTypeBySubjVar2_Pvalues)<-row.names(AveragePrimary_CellTypeIndex)
row.names(CellTypeBySubjVar2_Betas)<-row.names(AveragePrimary_CellTypeIndex)
row.names(CellTypeBySubjVar2_Tstat)<-row.names(AveragePrimary_CellTypeIndex)

head(CellTypeBySubjVar2_Pvalues)


for(i in c(1:length(AveragePrimary_CellTypeIndex[,1]))){
  
  temp<-summary.lm(lm(AveragePrimary_CellTypeIndex[i,]~Strain+Cohort))
  
  CellTypeBySubjVar2_Betas[i,]<-temp$coefficients[,1]
  CellTypeBySubjVar2_Tstat[i,]<-temp$coefficients[,3]
  CellTypeBySubjVar2_Pvalues[i,]<-temp$coefficients[,4]
  
}

write.csv(data.frame(CellTypeBySubjVar2_Betas, CellTypeBySubjVar2_Tstat, CellTypeBySubjVar2_Pvalues),"CellTypeBySubjVar2_Results.csv")

##############################

#Code for a general comparison of hippocampal differential expression across rats bred for internalizing or externalizing behavior.
#This code uses the differentially expressed gene list (Strain FDR<0.05) in the code above for the results for Blaveri et al. 2010 (PLOS One) 
#It also uses lists of differentially expressed genes that come directly from publications. To make these gene lists useful:
#I removed irregular formatting (e.g. footnotes on gene symbols)
#I removed gene symbols that were NA or multi-annotated. 
# I updated the gene symbol annotation for three of the studies that appeared to have older annotation (Zhang et al. 2005, Sabriego et al. 2013, Garafola and Henn et al. 2014). 
#Then I compiled all of the differentially expressed genes from each publication into a master excel sheet


#Reading in the master sheet:

setwd("~/Documents/Microarray Gen/HRLR/SecondaryAnalyses")

SelectivelyBredRats_HippocampalDE<-read.csv("SelectivelyBredRats_HippocampalDE.csv", header=T, stringsAsFactors = F)
str(SelectivelyBredRats_HippocampalDE)

str(SelectivelyBredRats_HippocampalDE)

table(SelectivelyBredRats_HippocampalDE$Citation)

# Andrus et al. 2012 (Molecular Psychiatry)                                       Birt et al.  submitted 
# 471                                                          192 
# Blaveri et al. 2010 (PLOS One)                                        Diaz-Moran et al. 2013 
# 908                                                          209 
# Garafola and Henn 2014 (Brain Research)                                Meckes et al. 2018 (PLOS ONE) 
# 13                                                          973 
# Raghavan et al. 2018 (Frontiers in Endocrinology)          Sabariego et al. 2013 (Behavioural Brain Research)  
# 46                                                           10 
# Wilhelm et al. 2013 (Pharmacology Biochemistry and Behavior)                 Zhang et al. 2005 (Genes Brain and Behavior) 
# 224                                                           38 


table(SelectivelyBredRats_HippocampalDE$Model.Comparison)
# bLR vs. bHR Congenitally Helpless (cLH) and Helpless Resistant (cNLH) 
# 192                                                        13 
# Flinder Sensitive vs. Flinders Resistant                     Flinders Sensitive vs. Sprague-Dawley 
# 908                                                       224 
# NIH-HS High Anxiety vs. NIH-HS Low Anxiety                                               RLA vs. RHA 
# 209                                                        10 
# SLA vs. SHA                                              WKY vs. F344 
# 38                                                       973 
# WMI vs. WLI (females)                                       WMI vs. WLI (males) 
# 46                                                       471 

#I still need to combine the results reflecting the same symbol for each publication.

unique(SelectivelyBredRats_HippocampalDE$GeneSymbol_Updated[SelectivelyBredRats_HippocampalDE$Citation=="Andrus et al. 2012 (Molecular Psychiatry)"])

Number_UniqueGeneSymbols_HippocampalDE<-matrix(0, 10, 1)
row.names(Number_UniqueGeneSymbols_HippocampalDE)<-names(table(SelectivelyBredRats_HippocampalDE$Model.Comparison))
colnames(Number_UniqueGeneSymbols_HippocampalDE)<-"Number_UniqueGeneSymbols_HippocampalDE"

for(i in c(1:10)){
  Number_UniqueGeneSymbols_HippocampalDE[i,1]<-length(unique(SelectivelyBredRats_HippocampalDE$GeneSymbol_Updated[SelectivelyBredRats_HippocampalDE$Model.Comparison==row.names(Number_UniqueGeneSymbols_HippocampalDE)[i]]))
}

Number_UniqueGeneSymbols_HippocampalDE

# Number_UniqueGeneSymbols_HippocampalDE
# bLR vs. bHR                                                                                  192
# Congenitally Helpless (cLH) and Helpless Resistant (cNLH)                                     13
# Flinder Sensitive vs. Flinders Resistant                                                     908
# Flinders Sensitive vs. Sprague-Dawley                                                        210
# NIH-HS High Anxiety vs. NIH-HS Low Anxiety                                                   209
# RLA vs. RHA                                                                                   10
# SLA vs. SHA                                                                                   36
# WKY vs. F344                                                                                 864
# WMI vs. WLI (females)                                                                         46
# WMI vs. WLI (males)                                                                          460

OverlapInHippocampalDE<-matrix(0, 10, 10)
colnames(OverlapInHippocampalDE)<-names(table(SelectivelyBredRats_HippocampalDE$Model.Comparison))
row.names(OverlapInHippocampalDE)<-names(table(SelectivelyBredRats_HippocampalDE$Model.Comparison))

for(i in c(1:10)){
  for(j in c(1:10)){
    OverlapInHippocampalDE[i,j]<-sum(unique(SelectivelyBredRats_HippocampalDE$GeneSymbol_Updated[SelectivelyBredRats_HippocampalDE$Model.Comparison==colnames(OverlapInHippocampalDE)[i]])%in%unique(SelectivelyBredRats_HippocampalDE$GeneSymbol_Updated[SelectivelyBredRats_HippocampalDE$Model.Comparison==colnames(OverlapInHippocampalDE)[j]]))
  }}
    
write.csv(OverlapInHippocampalDE, "OverlapInHippocampalDE.csv")

#The results are interesting - there is actually quite a bit of overlap. However, it is hard to intepret because each study uses different statistical thresholds, different sample sizes, and different transcriptional profiling platforms. So what we are seeing likely represents an enrichment over what would be expected randomly, but the magnitude of that enrichment is hard to quantify.


#So let's just determine which genes are popping up over and over again.

GeneSymbols_HippocampalDE_AcrossDatasets<-matrix(0, length(unique(SelectivelyBredRats_HippocampalDE$GeneSymbol_Updated)), 10)
row.names(GeneSymbols_HippocampalDE_AcrossDatasets)<-unique(SelectivelyBredRats_HippocampalDE$GeneSymbol_Updated)
colnames(GeneSymbols_HippocampalDE_AcrossDatasets)<-names(table(SelectivelyBredRats_HippocampalDE$Model.Comparison))

for(j in c(1:10)){
  for(i in c(1:nrow(GeneSymbols_HippocampalDE_AcrossDatasets))){
GeneSymbols_HippocampalDE_AcrossDatasets[i,j]<-row.names(GeneSymbols_HippocampalDE_AcrossDatasets)[i]%in%SelectivelyBredRats_HippocampalDE$GeneSymbol_Updated[SelectivelyBredRats_HippocampalDE$Model.Comparison==colnames(GeneSymbols_HippocampalDE_AcrossDatasets)[j]]
  }}

head(GeneSymbols_HippocampalDE_AcrossDatasets)

GeneSymbols_HippocampalDE_AcrossDatasets_Sum<-apply(GeneSymbols_HippocampalDE_AcrossDatasets, 1, sum)
max(GeneSymbols_HippocampalDE_AcrossDatasets_Sum)
#[1] 4

write.csv(cbind(GeneSymbols_HippocampalDE_AcrossDatasets, GeneSymbols_HippocampalDE_AcrossDatasets_Sum), "GeneSymbols_HippocampalDE_AcrossDatasets.csv")

#Ooh... I have a thought - let's add in direction of effect!

#the fact that there are multiple results for the same gene symbol in some datasets is causing problems.
# for(i in number of models...)
# tapply(SelectivelyBredRats_HippocampalDE$GeneSymbol_Updated average)
# need someplace to save it

str(SelectivelyBredRats_HippocampalDE)

# data.frame':	3084 obs. of  6 variables:
#  $ GeneSymbol                                                          : chr  "Tmem144" "Asb15" "Kif15" "Pkhd1l1" ...
# $ GeneSymbol_Updated                                                  : chr  "Tmem144" "Asb15" "Kif15" "Pkhd1l1" ...
# $ DirectionOfEffect..Up.Upregulated.in.Internalizing.Phenotype.       : chr  "Up" "Down" "Up" "Up" ...
# $ DirectionOfEffect_Numeric..1.Upregulated.in.Internalizing.Phenotype.: int  1 -1 1 1 -1 1 1 -1 1 1 ...
# $ Model.Comparison                                                    : chr  "bLR vs. bHR" "bLR vs. bHR" "bLR vs. bHR" "bLR vs. bHR" ...
# $ Citation                                                            : chr  "Birt et al.  submitted" "Birt et al.  submitted" "Birt et al.  submitted" "Birt et al.  submitted" ...

length(unique(SelectivelyBredRats_HippocampalDE$GeneSymbol_Updated))
#[1] 2598

SelectivelyBredRats_HippocampalDE_AveragedByGeneSymbol<-data.frame(GeneSymbol="NA", AverageDirectionOfEffect=0, Model.Comparison="NA", stringsAsFactors = F)
str(SelectivelyBredRats_HippocampalDE_AveragedByGeneSymbol)

for(j in c(1:10)){
  
  temp<-SelectivelyBredRats_HippocampalDE[SelectivelyBredRats_HippocampalDE$Model.Comparison==colnames(GeneSymbols_HippocampalDE_AcrossDatasets)[j],]  
  
  temp2<-tapply(temp$DirectionOfEffect_Numeric..1.Upregulated.in.Internalizing.Phenotype. , temp$GeneSymbol_Updated, mean)
  
  temp3<-data.frame(names(temp2), unname(temp2), rep(colnames(GeneSymbols_HippocampalDE_AcrossDatasets)[j], length(temp2)))
  
  colnames(temp3)<-c("GeneSymbol", "AverageDirectionOfEffect", "Model.Comparison")
  
  SelectivelyBredRats_HippocampalDE_AveragedByGeneSymbol<-rbind.data.frame(SelectivelyBredRats_HippocampalDE_AveragedByGeneSymbol, temp3)
  
  rm(temp, temp2, temp3)
}

str(SelectivelyBredRats_HippocampalDE_AveragedByGeneSymbol)

# 'data.frame':	2949 obs. of  3 variables:
#   $ GeneSymbol              : chr  "NA" "Aar2" "Abhd10" "Acad11" ...
# $ AverageDirectionOfEffect: num  0 -1 1 -1 1 -1 1 1 1 -1 ...
# $ Model.Comparison        : chr  "NA" "bLR vs. bHR" "bLR vs. bHR" "bLR vs. bHR" ...

head(SelectivelyBredRats_HippocampalDE_AveragedByGeneSymbol)
SelectivelyBredRats_HippocampalDE_AveragedByGeneSymbol<-SelectivelyBredRats_HippocampalDE_AveragedByGeneSymbol[-1,]

table(SelectivelyBredRats_HippocampalDE_AveragedByGeneSymbol$AverageDirectionOfEffect)
# -1 -0.666666666666667 -0.333333333333333                  0                  1 
# 1628                  1                  3                 20               1296 

#...so there are only 24 genes that had contradictory effects for different probes within the same dataset.

#Here's an example of one:
SelectivelyBredRats_HippocampalDE_AveragedByGeneSymbol[row.names(SelectivelyBredRats_HippocampalDE_AveragedByGeneSymbol)=="Uba5",]
# AverageDirectionOfEffect                      Model.Comparison
# Uba5                        0 Flinders Sensitive vs. Sprague-Dawley


GeneSymbols_HippocampalDE_AcrossDatasets2<-matrix(0, length(unique(SelectivelyBredRats_HippocampalDE_AveragedByGeneSymbol$GeneSymbol)), 10)
row.names(GeneSymbols_HippocampalDE_AcrossDatasets2)<-unique(SelectivelyBredRats_HippocampalDE_AveragedByGeneSymbol$GeneSymbol)
colnames(GeneSymbols_HippocampalDE_AcrossDatasets2)<-names(table(SelectivelyBredRats_HippocampalDE$Model.Comparison))

str(GeneSymbols_HippocampalDE_AcrossDatasets2)
# num [1:2598, 1:10] 0 0 0 0 0 0 0 0 0 0 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:2598] "Aar2" "Abhd10" "Acad11" "Acbd4" ...
# ..$ : chr [1:10] "bLR vs. bHR" "Congenitally Helpless (cLH) and Helpless Resistant (cNLH)" "Flinder Sensitive vs. Flinders Resistant" "Flinders Sensitive vs. Sprague-Dawley" ...

for(j in c(1:10)){
  for(i in c(1:nrow(GeneSymbols_HippocampalDE_AcrossDatasets2))){
    
    #Grab the data for a single model comparison j
  temp<-SelectivelyBredRats_HippocampalDE_AveragedByGeneSymbol[SelectivelyBredRats_HippocampalDE_AveragedByGeneSymbol$Model.Comparison==colnames(GeneSymbols_HippocampalDE_AcrossDatasets2)[j],]
  
  if(row.names(GeneSymbols_HippocampalDE_AcrossDatasets2)[i]%in%temp$GeneSymbol){
    
    GeneSymbols_HippocampalDE_AcrossDatasets2[i,j]<-temp$AverageDirectionOfEffect[temp$GeneSymbol==row.names(GeneSymbols_HippocampalDE_AcrossDatasets2)[i]]
    rm(temp)
  }else{
  GeneSymbols_HippocampalDE_AcrossDatasets2[i,j]<-NA
  rm(temp)}
  }
}

str(GeneSymbols_HippocampalDE_AcrossDatasets2)
head(GeneSymbols_HippocampalDE_AcrossDatasets2)


write.csv(GeneSymbols_HippocampalDE_AcrossDatasets2, "GeneSymbols_HippocampalDE_AcrossDatasets_wAVEDirection.csv")


#Hmm... after making that data.frame, we can also play with correlations between the directions of effects in the top genes from the different datasets:

write.csv(cor(GeneSymbols_HippocampalDE_AcrossDatasets2, use="pairwise.complete.obs"), "CorrelationMatrix_HippocampalDE_AcrossDatasets_wAVEDirection.csv")

dev.off()

colnames(GeneSymbols_HippocampalDE_AcrossDatasets2)
# [1] "bLR vs. bHR"                                               "Congenitally Helpless (cLH) and Helpless Resistant (cNLH)"
# [3] "Flinder Sensitive vs. Flinders Resistant"                  "Flinders Sensitive vs. Sprague-Dawley"                    
# [5] "NIH-HS High Anxiety vs. NIH-HS Low Anxiety"                "RLA vs. RHA"                                              
# [7] "SLA vs. SHA"                                               "WKY vs. F344"                                             
# [9] "WMI vs. WLI (females)"                                     "WMI vs. WLI (males)" 


#HR vs. LR compared to Flinders Sensitive vs. Resistant:
table(GeneSymbols_HippocampalDE_AcrossDatasets2[,1],GeneSymbols_HippocampalDE_AcrossDatasets2[,3])
#   -1   1
# -1 14  2
# 1  10  4
# The down-regulated genes can be consistent, the upregulated ones less so.

#HR vs. LR compared to Flinders Sensitive vs. Sprague:
table(GeneSymbols_HippocampalDE_AcrossDatasets2[,1],GeneSymbols_HippocampalDE_AcrossDatasets2[,4])

#   -1 0 1
# -1  3 0 0
# 1   5 0 2
# The down-regulated genes can be consistent, the upregulated ones less so.

#HR vs. LR compared to WKY-F344:
table(GeneSymbols_HippocampalDE_AcrossDatasets2[,1],GeneSymbols_HippocampalDE_AcrossDatasets2[,8])

#     -1 -0.666666666666667 -0.333333333333333  0  1
# -1  7                  0                  0  0  3
# 1  11                  0                  0  0 10
# The down-regulated genes can be consistent, the upregulated ones less so.

#Flinders Sensitive vs. Resistant vs. Fliders Sensitive vs. Sprague
table(GeneSymbols_HippocampalDE_AcrossDatasets2[,3],GeneSymbols_HippocampalDE_AcrossDatasets2[,4])
#     -1  0  1
# -1 20  0  6
# 1   2  0 10
#So in the Flinders datasets, both the upregulated and down-regulated genes tend to be relatively consistent.

#WKY vs. F344 vs. WMI vs. WLI (males):
table(GeneSymbols_HippocampalDE_AcrossDatasets2[,8],GeneSymbols_HippocampalDE_AcrossDatasets2[,10])
#                   -1  0  1
# -1                 20  0  8
# -0.666666666666667  0  0  0
# -0.333333333333333  0  0  0
# 0                   0  0  0
# 1                   2  0  5
#So in the WKY based models, the genes that are down-regulated in association with internalizing behavior are relatively consistent, but the upregulated genes are not.


##############################

#Let's make some Venn Diagrams illustrating the intersection between different selectively-bred rat HC DE datasets:

library(VennDiagram)

#The intersection between our HRLR genes and the top Flinders Sensitive vs. Flinders Resistant genes was examined first in Excel (spreadsheet GeneBySubjVar2DF.xlsx), spreadsheet: vsHRLR.
#Then I went back and used more consistent data cleaning/filtering across datasets (see code above) and redid the Venn Diagrams. 

#There are three decision points in the illustration of the results:

#1) Use just the 76 FDR<0.05 genes for HRLR (from the figure - including Bmp4, since it was FDR<0.05 in PCR) or use a slightly wider set (192 w/ FDR<0.05).

#2) Include genes that are not represented in both datasets in the counts?  I'm inclined to say yes, simply for practical reasons - we do not have the full data for other comparisons (e.g., w/ the Redei lab data), so it would be nice to use a consistent methodology.

#3) Intersect just two datasets at a time or more?  I'm inclined to do 2 first - I think three+ may get a little overwhelming, although it might be helpful to include the intersection w/ the Flinders Sensitive vs. Sprague results too...

#4) Illustrate intersection just based on FDR cut-offs or also expected direction of effect?  I'm inclined to use the Venn Diagram to illustrate intersection based on FDR and use colored text to illustrate shared direction of effect. Let's see how it looks...

#Old version (before compiling Master Database):
# pdf("VennDiagram_HRvsLR_vs_FlinderSvsR_Reanalyzed_FDR05.pdf", height=4.5, width=4.5)
# grid.newpage()
# draw.pairwise.venn(76, 913, 14, category = c("bLR vs. bHR", "Flinders Sensitive vs. Resistant"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2),cat.fontfamily = rep("sans", 2))
# dev.off()

#For calculating enrichment, it would be good to determine which genes are actually shared between the datasets. This gets complicated when including the developmental data, and also may be biased by negative results in both datasets related to genes that are not actually expressed in the hippocampus being represented in the microarray results.

#To rid of both issues:

#1) Trimmed out genes not found in our HRLR_HC_RNASeq metaanalyses (i.e., genes that are likely to be unexpressed or lowly expressed). There are 10,804 genes with results in both the Flinders and HRLR_HC_RNASeq, 787 of which have FDR<0.05 in Flinders, 72 have FDR<0.05 in the HRLR adult meta-analysis, and 13 are shared between the two. 

#Here's a Venn Diagram illustrating that: 

pdf("VennDiagram_HRvsLR_vs_FlinderSvsR_Reanalyzed_FDR05_OnlyAdultHCExpressed.pdf", height=4.5, width=4.5)
grid.newpage()
draw.pairwise.venn(72, 787, 13, category = c("bLR vs. bHR", "Flinders Sensitive vs. Resistant"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2),cat.fontfamily = rep("sans", 2))
dev.off()

#Formatted for figure-making
pdf("VennDiagram_HRvsLR_vs_FlinderSvsR_Reanalyzed_FDR05_OnlyAdultHCExpressed_NoNames.pdf", height=4.5, width=4.5)
grid.newpage()
draw.pairwise.venn(72, 787, 13, category = c("", ""), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2),cat.fontfamily = rep("sans", 2), cex=1.5)
dev.off()


ContingencyTable<-cbind(c((10804-(787-13)-(72-13)), (787-13)), c((72-13), 13))
colnames(ContingencyTable)<-c("NotSig_HRLR", "Sig_HRLR")
rownames(ContingencyTable)<-c("NotSig_Flinders", "Sig_Flinders")
ContingencyTable

                  #NotSig_HRLR Sig_HRLR
#NotSig_Flinders        9971       59
#Sig_Flinders            774       13

# %NotSigHRLR but Sig in Flinders
774/(9971+774)
#[1] 0.0720335

# %SigHRLR and Sig in Flinders
13/(59+13)
#[1] 0.1805556

fisher.test(ContingencyTable)

# Fisher's Exact Test for Count Data
# 
# data:  ContingencyTable
# p-value = 0.00186
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 1.421384 5.262875
# sample estimates:
# odds ratio 
# 2.838305 


############################


#Alright, before compiling the master database, I did an initial comparison with Wilhelm's Flinders Sensitive vs. Sprague results. To do this, I trimmed out the results that completely lacked any annotation or were annotated to multiple genes. And then removed duplicates, double-checking whether they included effects in opposite directions. That left 211 DE genes.
# pdf("VennDiagram_FlindersSvsSprague_vs_FlinderSvsR_Reanalyzed_FDR05.pdf", height=4.5, width=4.5)
# grid.newpage()
# draw.pairwise.venn(210, 913, 38, category = c("Flinders Sensitive vs. Sprague", "Flinders Sensitive vs. Resistant"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2),cat.fontfamily = rep("sans", 2))
# dev.off()

#Trimming it down using the same method as earlier (to try to only catch genes that we are confident are expressed in the adult hippocampus)
#Only 159 of the FlindersSvsSprague DE genes were found in both the FlindersReanalysis and HRLR_HC_RNASeq. 
#Interestingly, only two of the genes that overlapped in these datasets were not detected as expressed in our HRLR_HC_RNASeq analysis.
#Hmm... but one of the genes that overlapped with the top HRLR genes is no longer in the dataset (Fmo5) - it must not have included on the platform used for the FlindersSvsFlindersR study. 
#And we have no idea what other genes were present in the rest of the analysis because the full results weren't released, so we can't do an enrichment calculation. :(


#Before compiling the master database, I also did an initial comparison of bHR vs. bLR with Wilhelm's Flinders Sensitive vs. Sprague results

# pdf("VennDiagram_FlindersSvsSprague_vs_HRvsLR_FDR05.pdf", height=4.5, width=4.5)
# grid.newpage()
# draw.pairwise.venn(76, 210, 5, category = c("bLR vs. bHR", "Flinders Sensitive vs. Sprague"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), cat.fontfamily = rep("sans", 2))
# dev.off()


#... and then combined them all into a three-way Venn Diagram:
#Note: if we do a three-way, do we want to use FDR<0.05 or FDR<0.10?

# pdf("VennDiagram_bHRBLR_FlindersSvsSprague_vs_FlinderSvsR_Reanalyzed_FDR05.pdf", height=6, width=6)
# grid.newpage()
# overrideTriple=T
# #This override allows it to be scaled by set size
# draw.triple.venn(76, 210, 913, 5, 38, 14, 1, category = c("bLR vs. bHR ", "Flinders Sensitive vs. Sprague-Dawley", "Flinders Sensitive vs. Flinders Resistant"),  lty = rep("blank", 3), fill = c("light blue", "pink", "yellow"), alpha = rep(0.5, 3), cat.pos = c(0, 0, 180), cat.dist = rep(0.025, 3), cat.fontfamily = rep("sans", 3))
# dev.off()
# 
# #Version without names so that I can improve the formatting outside of R...
# 
# pdf("VennDiagram_bHRBLR_FlindersSvsSprague_vs_FlinderSvsR_Reanalyzed_FDR05_NoNames.pdf", height=6, width=6)
# grid.newpage()
# overrideTriple=T
# #This override allows it to be scaled by set size
# draw.triple.venn(76, 210, 913, 5, 38, 14, 1, category = c(" ", " ", " "),  lty = rep("blank", 3), fill = c("light blue", "pink", "yellow"), alpha = rep(0.5, 3), cat.pos = c(0, 0, 180), cat.dist = rep(0.025, 3), cat.fontfamily = rep("sans", 3), cex=1.5)
# dev.off()

#After compiling the master database (with consistent annotation, filtering of results, etc), these numbers changed a small amount:

#New version

pdf("VennDiagram_bHRBLR_FlindersSvsSprague_vs_FlinderSvsR_Reanalyzed_FDR05_NoNames_v2.pdf", height=6, width=6)
grid.newpage()
overrideTriple=T
#This override allows it to be scaled by set size
draw.triple.venn(76, 210, 908, 5, 38, 14, 1, category = c(" ", " ", " "),  lty = rep("blank", 3), fill = c("light blue", "pink", "yellow"), alpha = rep(0.5, 3), cat.pos = c(0, 0, 180), cat.dist = rep(0.025, 3), cat.fontfamily = rep("sans", 3), cex=1.5)
dev.off()


#cut parameters
#reverse=FALSE, euler.d = TRUE, scaled = TRUE, 
#, cat.default.pos="text" 

#Comparison of bHR vs. bLR with Mecke's WKY vs. F34 results

#916 unique gene symbols in Mecke's top results after removing 115 duplicates.

# pdf("VennDiagram_WKYvsF34_vs_HRvsLR_FDR05.pdf", height=4.5, width=4.5)
# grid.newpage()
# draw.pairwise.venn(76, 916, 14, category = c("bHR vs. bLR", "WKY vs. F34"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2),cat.fontfamily = rep("sans", 2))
# dev.off()
# 
# #Version without names so that I can change the formatting outside of R...
# #Also made to match formatting on 3-way:
# pdf("VennDiagram_WKYvsF34_vs_HRvsLR_FDR05_NoNames.pdf", height=5.5, width=5.5)
# grid.newpage()
# draw.pairwise.venn(76, 916, 14, category = c("", ""), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2),cat.fontfamily = rep("sans", 2), cex=1.5)
# dev.off()

#Hmmm... I'm realizing that we should probably throw out the 51 results in Mecke that are multi-annotated (since there is no hope of those matching with our dataset)

pdf("VennDiagram_WKYvsF34_vs_HRvsLR_FDR05_NoNames_NoDoubleAnnoation.pdf", height=5.5, width=5.5)
grid.newpage()
draw.pairwise.venn(76, 864, 14, category = c("", ""), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2),cat.fontfamily = rep("sans", 2), cex=1.5)
dev.off()


#************************

#Comparison of bHR vs. bLR with Andrus' WMI vs. WLI results:

#This initial analysis was run before I compiled the Master Database:
#For comparison, footnotes were removed from the gene symbols.
#They had only 11 duplicated gene names in their DE results, no multi-annotated or unannotated results.
#So that leaves 460 DE genes.

# pdf("VennDiagram_WMIvsWLI_vs_HRvsLR_FDR05.pdf", height=4.5, width=4.5)
# grid.newpage()
# draw.pairwise.venn(76, 460, 1, category = c("bHR vs. bLR", "WMI vs. WLI"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2),cat.fontfamily = rep("sans", 2))
# dev.off()
# 
 
pdf("VennDiagram_WMIvsWLI_vs_HRvsLR_FDR05_NoNames.pdf", height=4.5, width=4.5)
grid.newpage()
draw.pairwise.venn(76, 460, 1, category = c("", ""), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2),cat.fontfamily = rep("sans", 2), cex=1.5)
 dev.off()

#These numbers were unchanged after compiling the master database:


#It turns out there is another WMI vs. WLI study (in females). Let's do a three way overlap with that.
#This version was made before compiling the master database.
 
# pdf("VennDiagram_bHRBLR_WMIvsWLI_MaleFemale_FDR05.pdf", height=5, width=5)
# grid.newpage()
# overrideTriple=TRUE
# #This override should allow it to be scaled by set size... but it doesn't seem to be working
# draw.triple.venn(76, 460, 40, 1, 1, 0, 0, category = c("bLR vs. bHR ", "WMI vs. WLI (Male)", "WMI vs. WLI (Female)"),  lty = rep("blank", 3), euler.d =TRUE, scaled = TRUE, fill = c("light blue", "pink", "yellow"), alpha = rep(0.5, 3), cat.pos = c(0, 0, 180), cat.dist = rep(0.025, 3), cat.fontfamily = rep("sans", 3))
# dev.off()
# 
# 
# pdf("VennDiagram_bHRBLR_WMIvsWLI_MaleFemale_FDR05_NoNames.pdf", height=5, width=5)
# grid.newpage()
# overrideTriple=TRUE
# #This override should allow it to be scaled by set size... but it doesn't seem to be working
# draw.triple.venn(76, 460, 40, 1, 1, 0, 0, category = c("", "", ""),  lty = rep("blank", 3), euler.d =TRUE, scaled = TRUE, fill = c("light blue", "pink", "yellow"), alpha = rep(0.5, 3), cat.pos = c(0, 0, 180), cat.dist = rep(0.025, 3), cat.fontfamily = rep("sans", 3))
# dev.off()

#Fiddling around with the code, I can't get the circles to scale properly. So I just ended up adding in the WMI vs. WLI female circle to the previously-outputted male overlap diagram by hand graphically. Sigh.
#After making the master database, these numbers changed a little (fixed in figure)


#hmmm... some general thoughts: 

#1) It might be helpful to output versions that are black & white and missing labels so that we have more flexibility with downstream formatting.
#Nevermind - looked too awkward.


#2) Overlap with HRvsLR FDR<0.10 might be interesting. It might also be easier to plot for the three way Venn.

#Messed around with it. Looks  prettier... but the number of genes found in the overlap may be overwhelming for the reader. And I don't want to scoop Hozi's meta-analysis story.  So maybe we'll stick with the simpler version.



#O.k. - I'm getting ambitious now: Can we do a 4-way Venn Diagram?

pdf("VennDiagram_4way_WKYvsF34_vs_HRvsLR_vs_Flinders_FDR05_NoNames.pdf", height=8.5, width=8.5)
grid.newpage()

draw.quad.venn(76, 210, 913, 864,  5, 14, 14, 38, n24,
               n34, n123, n124, n134, n234, n1234, 
               
               
               1,
               
               category = rep("",
                              4), lwd = rep(2, 4), lty = rep("solid", 4), col =
                 rep("black", 4), fill = NULL, alpha = rep(0.5, 4),
               label.col = rep("black", 15), cex = rep(1, 15),
               fontface = rep("plain", 15), fontfamily = rep("serif",
                                                             15), cat.pos = c(-15, 15, 0, 0), cat.dist = c(0.22,
                                                                                                           0.22, 0.11, 0.11), cat.col = rep("black", 4), cat.cex
               = rep(1, 4), cat.fontface = rep("plain", 4),
               cat.fontfamily = rep("serif", 4), cat.just =
                 rep(list(c(0.5, 0.5)), 4), rotation.degree = 0,
               rotation.centre = c(0.5, 0.5), ind = TRUE, cex.prop =
                 NULL, print.mode = "raw", sigdigs = 3, direct.area =
                 FALSE, area.vector = 0, ...)


#... and I'm running out of time to finish my revisions, so I guess this more ambitious coding will need to wait until some other time...


