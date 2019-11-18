## Meta-analysis on Adult HR/LR Data

#This script contains code running the meta-analysis of 5 adult HR/LR datasets
#To run the meta-analysis, t-test of phenotype effect on gene expression is read-in for each dataset. Cohen's d is then extracted from t-test values and joined together for the meta-analysis.
#Note: one dataset, the RNA-seq F29 dataset, is the only dataset not re-analyzed prior to the meta-analysis as original data files were unavailable. Instead, Cohen's d is extracted from the p-values obtained through a prior analysis.

### Loading necessary libraries

library(plyr)
library(metafor)
library(compute.es)
library(bmp)

#This package is necessary to run the tes function

source("https://bioconductor.org/biocLite.R")
biocLite("multtest")
library(multtest)


## Reading in data

# Reading in test statistic data

RNAseq_F29 <- read.csv("data/RNAseq F29/LRHR_rgsc34_th2out_g1tx_diff_AdLvsAdH_out_gene_exp.csv")
RNAseq_F37 <- read.csv("output/Blandino RNAseq F37/BlandinoTT.csv")
RNAseq_F43 <- read.csv("output/Aydin RNAseq F43/F43_AdultTT.csv")
Affy_F4    <- read.csv("output/Stead Affy F4/Stead_AdultF4TT.csv")
NimbleGen_F34 <- read.csv("output/Alabama New Colony Development/newdevelopmentdataAdultTT.csv")

# Reading in cell-type data

cellRNAseq_F37 <- read.csv("output/Blandino RNAseq F37/BlandinoCellTypeTT_forMeta.csv")
cellRNAseq_F43 <- read.csv("output/Aydin RNAseq F43/VehCellTypeTT_forMeta.csv")
cellAffy_F4 <- read.csv("output/Stead Affy F4/Stead_AdultF4CellTypeTT_forMeta.csv")
cellNimbleGen_F34 <- read.csv("output/Alabama New Colony Development/NewColonyDevelopmentCellTypeTT_forMeta.csv")


#setting directories

outDir <- "output/MetaAnalysis Adult/"
outForestPlots <- "output/MetaAnalysis Adult/ForestPlots_topGenes_Adult/"
outGOIPlots <- "output/MetaAnalysis Adult/Forest Plots for GOI Adult/"
outCellPlots <- "output/MetaAnalysis Adult/ForestPlots_CellType_Adult/"


## Extracting Effect Size
### The following code converts t-test values into effect sizes d (Cohen's D) and g (Hedge's g) and extracts the corresponding variances for each of the five adult datasets. 


# RNAseq_F29

colnames(RNAseq_F29)

#need to omit data with no value (p-value=1)
rnaseq<-subset(RNAseq_F29, RNAseq_F29$p_value!=1)
dim(rnaseq)

sum(duplicated(rnaseq$gene))
#lots of entries have associated gene IDs but no gene symbol

#removing genes with no gene symbol

rnaseq <- subset(rnaseq, !(rnaseq$gene=="-"))

###Removing extreme values (genes where transcript info was only present for one group HR or LR)

rnaseq_nobad<-subset(rnaseq, rnaseq$test_stat > -1.79769e+308 & rnaseq$test_stat < 1.79769e+308)


sum(duplicated(rnaseq_nobad$gene))
#still lots of duplicated genes. Many of these seem to have very different expression values.

#Averaging across duplicate genes
colnames(rnaseq_nobad)

rnaseq_avg<-tapply(rnaseq_nobad[,11], as.character(rnaseq_nobad[,3]), function(y) mean(y))


#extracting the test statistic averaged by gene
rnaseqTT<-cbind.data.frame(rnaseq_avg)

CohDRNAseq<-matrix(0, length(rnaseqTT[,1]), 4)

for(i in 1:length(rnaseqTT[,1])){
  CohDRNAseq[i, 1]<-tes(rnaseqTT[i,], 2,2)$d
  CohDRNAseq[i, 2]<-tes(rnaseqTT[i,], 2,2)$'var.d'
  CohDRNAseq[i, 3]<-tes(rnaseqTT[i,], 2,2)$g
  CohDRNAseq[i, 4]<-tes(rnaseqTT[i,], 2,2)$'var.g'
}

colnames(CohDRNAseq)[1:4]<-c("d RNAseq F29", "var d RNAseq F29", "g RNAseq F29", "var g RNAseq F29")


#adding in gene symbol information
CohDRNAseq<-cbind.data.frame(CohDRNAseq, Symbol=row.names(rnaseqTT))


#Flipping direction of effect so it matches other datasets (HR positive LR negative)
CohDRNAseq[,1] <- (-1*CohDRNAseq[[1]])



# RNAseq_F37

F37TT <- cbind.data.frame(RNAseq_F37$t.test)

CohDF37<-matrix(0, length(RNAseq_F37[,1]), 4)

for(i in 1:length(F37TT[,1])){
  CohDF37[i, 1]<-tes(F37TT[i,], 6,6)$d
  CohDF37[i, 2]<-tes(F37TT[i,], 6,6)$'var.d'
  CohDF37[i, 3]<-tes(F37TT[i,], 6,6)$g
  CohDF37[i, 4]<-tes(F37TT[i,], 6,6)$'var.g'
}


colnames(CohDF37)<-c("d F37", "var d F37", "g F37", "var g F37")
CohDF37 <- cbind.data.frame(Symbol=RNAseq_F37$X, CohDF37)



#RNAseq_F43

#Getting Cohen's d from Cigdem Aydin's dataset
F43TT<-cbind.data.frame(RNAseq_F43$t.test)

CohDF43<-matrix(0, length(F43TT[,1]), 4)

for(i in 1:length(F43TT[,1])){
  CohDF43[i,1]<-tes(F43TT[i,], 5,5)$d
  CohDF43[i,2]<-tes(F43TT[i,], 5,5)$'var.d'
  CohDF43[i, 3]<-tes(F43TT[i,], 5,5)$g
  CohDF43[i, 4]<-tes(F43TT[i,], 5,5)$'var.g'
}


colnames(CohDF43)<-c("d F43", "var d F43", "g F43", "var g F43")
CohDF43<-cbind.data.frame(Symbol=RNAseq_F43$X, CohDF43)


#Affy_F4

#extracting TTstat for tes function
affyF4TT<-cbind.data.frame(Affy_F4$t.test)


CohDAffyF4<-matrix(0, length(Affy_F4[,1]), 4)

for(i in 1:length(affyF4TT[,1])){
  CohDAffyF4[i, 1]<-tes(affyF4TT[i,], 5,5)$d
  CohDAffyF4[i, 2]<-tes(affyF4TT[i,], 5,5)$'var.d'
  CohDAffyF4[i, 3]<-tes(affyF4TT[i,], 5,5)$g
  CohDAffyF4[i, 4]<-tes(affyF4TT[i,], 5,5)$'var.g'
}

colnames(CohDAffyF4)<-c("d Affy F4", "var d Affy F4", "g Affy F4", "var g Affy F4")
CohDAffyF4<-cbind.data.frame(Symbol=(Affy_F4$'X'), CohDAffyF4)


#NimbleGen_F34

ncF34TT<-cbind.data.frame(NimbleGen_F34[,2])

CohDncF34<-matrix(0, length(ncF34TT[,1]), 4)

for(i in 1:length(ncF34TT[,1])){
  CohDncF34[i, 1]<-tes(ncF34TT[i,], 5,5)$d
  CohDncF34[i, 2]<-tes(ncF34TT[i,], 5,5)$'var.d'
  CohDncF34[i, 3]<-tes(ncF34TT[i,], 5,5)$g
  CohDncF34[i, 4]<-tes(ncF34TT[i,], 5,5)$'var.g'
}

colnames(CohDncF34)<-c("d NC F34", "var d NC F34", "g NC F34", "var g NC F34")
CohDncF34<-cbind.data.frame(Symbol=NimbleGen_F34$X, CohDncF34)


## Meta-analysis
#The following code joins together the effect size dataframes for each of the five adult datasets and runs the meta-analysis

dfs<-list(
  CohDF43,
  CohDF37,
  CohDncF34,
  CohDAffyF4,
  CohDRNAseq
)

meta<-join_all(dfs, by="Symbol", type="full", match="all")
dim(meta)

sum(duplicated(meta$Symbol))


#Checking number of overlapping genes
allFive<-join_all(dfs, by="Symbol", type="inner")
length(!duplicated(allFive$Symbol))
#Number of genes found in ALL FIVE datasets.

colnames(meta)

#####Looking at correlation between cohen's d of different datasets

colnames(allFive) 
#Have to use one with only genes found in all datasets because cor can't handle NAs

#Grabbing only cohen's d columns
temp<-cbind(allFive$`d F43`, allFive$`d F37`, allFive$`d NC F34`, 
            allFive$`d RNAseq`, allFive$`d Affy F4`)

temp1<-as.matrix(temp)
cormatrix<-cor(temp1)

#output cor matrix
write.csv(cormatrix, paste0(outDir, "CohensdCorMatrix_AdultMeta.csv"))


#### Running Meta-Analysis

metaOutput<-matrix(0, length(meta$Symbol), 6)

for(i in 1:length(meta$Symbol)){
  effect<-c(meta$`d F43`[i], meta$`d F37`[i], meta$`d NC F34`[i], meta$`d RNAseq F29`[i], meta$`d Affy F4`[i])
  var<-c(meta$`var d F43`[i], meta$`var d F37`[i], meta$`var d NC F34`[i], meta$`var d RNAseq F29`[i], meta$`var d Affy F4`[i])
  if(sum(is.na(effect))>3){}
  else{
    metaOutput[i, 1]<-rma.mv(effect, var)$b #gives estimate
    metaOutput[i, 2]<-rma.mv(effect, var)$se #gives standard error
    metaOutput[i, 3]<-rma.mv(effect, var)$pval #gives pval
    metaOutput[i, 4]<-rma.mv(effect, var)$ci.lb #gives confidence interval lower bound
    metaOutput[i, 5]<-rma.mv(effect, var)$ci.ub #gives confidence interval upper bound
  }
  metaOutput[i, 6] <- sum(5 - sum(is.na(effect)))
}

colnames(metaOutput)<-c("estimate", "SE", "pval", "CI_lb", "CI_ub", "datasets")

metaOutputSymbol<-cbind.data.frame(GeneSymbol=meta$Symbol, metaOutput)


sum(metaOutputSymbol$pval==0) 
#Some genes are found in only one dataset and 
#therefor have no values associated with them from the meta-analysis

metaOutputSymbol <- subset(metaOutputSymbol, !(metaOutputSymbol$pval==0))
dim(metaOutputSymbol)
#total genes present in at least two datasets

write.csv(metaOutputSymbol, paste0(outDir, "AdultMetaAnalysisOutput.csv"))



## Multiple Comparison corrections
#The following code applies two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli) 
#Meta-analysis output with adjusted p-values is then output along with effect size information.


colnames(metaOutputSymbol)

tempPvalAdjMeta<-mt.rawp2adjp(metaOutputSymbol$pval, proc=c("BH", "BY"))
metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]

#adjusted pvalue object is in same orientation as metaoutput so can simply be binded together


metaAnalysisOutputFDR<-cbind.data.frame(metaPvalAdj, metaOutputSymbol)
dim(metaAnalysisOutputFDR)


write.csv(metaAnalysisOutputFDR, paste0(outDir, "AdultMetaAnalysisOutputFDR.csv"))


#Extracting sig genes (FDR < 0.05)

adultMetaSigGeneList<-subset(metaAnalysisOutputFDR, metaAnalysisOutputFDR$'BH' < 0.05)


write.csv(adultMetaSigGeneList, paste0(outDir, "AdultMetaSigGeneList.csv"))


#putting coh d and adjusted pval info together

#Renaming Symbol column from meta dataframe to join with meta output dataframe
temp<-meta
colnames(temp)[1] <- "GeneSymbol"

cohDandPval<-join(metaAnalysisOutputFDR, temp, by="GeneSymbol", type="inner")

write.csv(cohDandPval, paste0(outDir, "AdultMeta_CohDandPval.csv"))


## Cell Type MetaAnalysis
#The following code performs the meta-analysis on the cell type test statistic values, very similar to previous code.

# Extracting Cohen's D from all

#RNAseq_F43          
#sample size is 5
cellF43TT <- cbind.data.frame(cellRNAseq_F43$t.test)

CohD_CellType_F43<-matrix(0, length(cellF43TT[,1]), 4)

for(i in 1:length(cellF43TT[,1])){
  CohD_CellType_F43[i, 1]<-tes(cellF43TT[i,], 5,5)$d
  CohD_CellType_F43[i, 2]<-tes(cellF43TT[i,], 5,5)$'var.d'
  CohD_CellType_F43[i, 3]<-tes(cellF43TT[i,], 5,5)$g
  CohD_CellType_F43[i, 4]<-tes(cellF43TT[i,], 5,5)$'var.g'
}

row.names(CohD_CellType_F43)<-(cellRNAseq_F43$X)
colnames(CohD_CellType_F43)<-c("d F43", "var d F43", "g F43", "var g F43")



#RNAseq_F37         
#sample size is 6
cellF37TT <- cbind.data.frame(cellRNAseq_F37$t.test)

CohD_CellType_F37<-matrix(0, length(cellF37TT[,1]), 4)

for(i in 1:length(cellF37TT[,1])){
  CohD_CellType_F37[i, 1]<-tes(cellF37TT[i,], 6,6)$d
  CohD_CellType_F37[i, 2]<-tes(cellF37TT[i,], 6,6)$'var.d'
  CohD_CellType_F37[i, 3]<-tes(cellF37TT[i,], 6,6)$g
  CohD_CellType_F37[i, 4]<-tes(cellF37TT[i,], 6,6)$'var.g'
}

row.names(CohD_CellType_F37)<-cellRNAseq_F37$X
colnames(CohD_CellType_F37)<-c("d F37", "var d F37", "g F37", "var g F37")

#New Colony         
#sample size is 5
cellF34TT <- cbind.data.frame(cellNimbleGen_F34$Adult.TT)

CohD_CellType_NCF34<-matrix(0, length(cellF34TT[,1]), 4)

for(i in 1:length(cellF34TT[,1])){
  CohD_CellType_NCF34[i, 1]<-tes(cellF34TT[i,], 5,5)$d
  CohD_CellType_NCF34[i, 2]<-tes(cellF34TT[i,], 5,5)$'var.d'
  CohD_CellType_NCF34[i, 3]<-tes(cellF34TT[i,], 5,5)$g
  CohD_CellType_NCF34[i, 4]<-tes(cellF34TT[i,], 5,5)$'var.g'
}

row.names(CohD_CellType_NCF34)<-cellNimbleGen_F34$X
colnames(CohD_CellType_NCF34)<-c("d NC F34", "var d NC F34", "g NC F34", "var g NC F34")


#John Stead          
#sample size is 5
cellF4TT <- cbind.data.frame(cellAffy_F4$Stead_TT)

CohD_CellType_F4<-matrix(0, length(SteadTT[,1]), 4)

for(i in 1:length(cellF4TT[,1])){
  CohD_CellType_F4[i, 1]<-tes(cellF4TT[i,], 5,5)$d
  CohD_CellType_F4[i, 2]<-tes(cellF4TT[i,], 5,5)$'var.d'
  CohD_CellType_F4[i, 3]<-tes(cellF4TT[i,], 5,5)$g
  CohD_CellType_F4[i, 4]<-tes(cellF4TT[i,], 5,5)$'var.g'
}


row.names(CohD_CellType_F4)<-cellAffy_F4$X
colnames(CohD_CellType_F4)<-c("d Affy F4", "var d Affy F4", "g Affy F4", "var g Affy F4")



#### Cbinding Cohen's D for all datasets
#since each cohen's d dataset is in the same order and for the same 10 cell types, they can simply be binded together

CohD_Adult_CellType<-cbind.data.frame(CohD_CellType_F43, CohD_CellType_F37, CohD_CellType_NCF34, CohD_CellType_F4)

CellTypeMetaAnalysisOutput<-matrix(0, 10, 5)

for(i in 1:10){
  effect<-c(CohD_Adult_CellType$`d F43`[i], CohD_Adult_CellType$`d F37`[i], 
            CohD_Adult_CellType$`d NC F34`[i], CohD_Adult_CellType$`d Affy F4`[i])
  var<-c(CohD_Adult_CellType$`var d F43`[i], CohD_Adult_CellType$`var d F37`[i],
         CohD_Adult_CellType$`var d NC F34`[i], CohD_Adult_CellType$`var d Affy F4`[i])
  CellTypeMetaAnalysisOutput[i, 1]<-rma.mv(effect, var)$b #gives estimate
  CellTypeMetaAnalysisOutput[i, 2]<-rma.mv(effect, var)$se #gives standard error
  CellTypeMetaAnalysisOutput[i, 3]<-rma.mv(effect, var)$pval #gives pval
  CellTypeMetaAnalysisOutput[i, 4]<-rma.mv(effect, var)$ci.lb #gives confidence interval lower bound
  CellTypeMetaAnalysisOutput[i, 5]<-rma.mv(effect, var)$ci.ub #gives confidence interval upper bound
}

colnames(CellTypeMetaAnalysisOutput)<-c("estimate", "SE", "pval", "CI_lb", "CI_ub")
row.names(CellTypeMetaAnalysisOutput)<-row.names(CohD_Adult_CellType)

#joining with Cohen's D info

colnames(CohD_Adult_CellType)
temp<-cbind.data.frame(CellTypeMetaAnalysisOutput, CohD_Adult_CellType)

write.csv(temp, paste0(outDir, "Adult_CellTypeMeta_CohDandPval.csv"))


###############

## Forest Plots

#First exporting all significant genes as .png files and then exporting genes of interests as TIFF files

#Subsetting significant genes
sigMeta<-subset(cohDandPval, cohDandPval$BH < 0.05)
length(sigMeta$GeneSymbol)
#[1] 74


##Forest Plots with official thesis names for datasets


#Forloop over all sig genes and create forest plots for them

for(i in 1:length(sigMeta$GeneSymbol)){
  
  png(paste0(outForestPlots, "Forest Plot ", paste((sigMeta$GeneSymbol)[i]), ".png"))
  
  effect<-c(sigMeta$`d Affy F4`[i], sigMeta$`d RNAseq F29`[i], sigMeta$`d NC F34`[i], 
            sigMeta$`d F37`[i], sigMeta$`d F43`[i])
  
  var<-c(sigMeta$`var d Affy F4`[i], sigMeta$`var d RNAseq F29`[i], sigMeta$`var d NC F34`[i], 
         sigMeta$`var d F37`[i], sigMeta$`var d F43`[i])
  
  forest.rma(rma.mv(effect, var), slab=c("MBNI_AffymetrixRae230_F4", "MBNI_RNASeq_F29",  
                                         "Alabama_NimbleGen_F34", "MBNI_RNASeq_F37 ", 
                                         "MBNI_RNASeq_F43"), xlim=c(-18, 17))
  
  text(-12.7, 5.7, "Investigator & Study", cex=1)
  text(-13, 6.5, "LR", pos=4, cex=1.5)
  mtext(paste(sigMeta$GeneSymbol[i]), line=-1.5, cex=2)
  text(11.7, 5.7, "Cohen's D [95% CI] ", cex=1)
  text(12, 6.5, "HR ", pos=2, cex=1.5)
  
  dev.off()
  
}


### Export high res forest plots for genes of interest

#subsetting genes of interest

temp<-subset(meta, meta$Symbol=="Bmp4" | meta$Symbol=="Ncan" | meta$Symbol=="Trhr" |
               meta$Symbol=="C1qa" | meta$Symbol=="Fos" | meta$Symbol=="C1qb" | 
               meta$Symbol=="C1qc" | meta$Symbol=="Ucp2" | meta$Symbol == "Tubg1" |
               meta$Symbol=="Mfge8" | meta$Symbol=="Fxyd7" | meta$Symbol=="Slc9a3r1")


#Export forest plots as TIFF files
for(i in 1:length(temp$Symbol)){
  
  tiff(paste0(outGOIPlots, "Forest Plot ", 
              paste((temp$Symbol)[i]), ".tiff"), res=300, compression="lzw",
       width = 8.5, height = 8, units = 'in')
  
  effect<-c(temp$`d Affy F4`[i], temp$`d RNAseq F29`[i], temp$`d NC F34`[i], 
            temp$`d F37`[i], temp$`d F43`[i])
  var<-c(temp$`var d Affy F4`[i], temp$`var d RNAseq F29`[i], temp$`var d NC F34`[i], 
         temp$`var d F37`[i], temp$`var d F43`[i])
  
  forest.rma(rma.mv(effect, var), slab=c("MBNI_AffymetrixRae230_F4", "MBNI_RNASeq_F29", 
                                         "Alabama_NimbleGen_F34", "MBNI_RNASeq_F37",
                                         "MBNI_RNASeq_F43"), 
             xlim=c(-18, 15), cex=1.6)
  text(-11, 5.6, "Investigator & Study", cex=1.6)
  text(-14, 6.6, "LR", cex=1.8)
  mtext(paste("Adult\n", temp$Symbol[i]), line=-1.5, cex=2.5)
  text(8, 5.6, "Cohen's D [95% CI]", cex=1.6)
  text(11, 6.6, "HR", cex=1.8)
  
  dev.off()
}

#Megan's version (for quickly extracting a few more forest plots)
#Output folder:
setwd("~/Documents/Microarray Gen/HRLR/ThesisMetaAnalysisOutput/ForestPlots/Prettier")

temp<-subset(AdultMeta, AdultMeta$GeneSymbol=="Cav1" | AdultMeta$GeneSymbol=="Sdc4"| AdultMeta$GeneSymbol=="Mfge8" | AdultMeta$GeneSymbol=="C1qa"| AdultMeta$GeneSymbol=="Trhr"| AdultMeta$GeneSymbol=="Ucp2")

temp<-subset(AdultMeta, AdultMeta$GeneSymbol=="Gpc1" |AdultMeta$GeneSymbol=="Sox2" | AdultMeta$GeneSymbol=="Hes5" | AdultMeta$GeneSymbol=="Cd24")

temp<-subset(AdultMeta, AdultMeta$GeneSymbol=="Mki67" | AdultMeta$GeneSymbol=="Sox9"| AdultMeta$GeneSymbol=="Tek"| AdultMeta$GeneSymbol=="Sp3"|AdultMeta$GeneSymbol=="Uhrf1" | AdultMeta$GeneSymbol=="C2cd3")

for(i in 1:length(temp$GeneSymbol)){
  
  tiff(paste("Forest Plot ", (temp$GeneSymbol)[i], ".tiff"), res=300, compression="lzw",
       width = 8.5, height = 8, units = 'in')
  
  effect<-c(temp$`d.Affy.F4`[i], temp$`d.RNAseq.F29`[i], temp$`d.NC.F34`[i], 
            temp$`d.F37`[i], temp$`d.F43`[i])
  var<-c(temp$`var.d.Affy.F4`[i], temp$`var.d.RNAseq.F29`[i], temp$`var.d.NC.F34`[i], 
         temp$`var.d.F37`[i], temp$`var.d.F43`[i])
  
  forest.rma(rma.mv(effect, var), slab=c("MBNI_AffymetrixRae230_F4", "MBNI_RNASeq_F29", 
                                         "Alabama_NimbleGen_F34", "MBNI_RNASeq_F37",
                                         "MBNI_RNASeq_F43"), 
             xlim=c(-18, 15), cex=1.6)
  text(-11, 6.6, "Investigator & Study", cex=1.6)
  text(-14, 7.6, "bLR", cex=1.8)
  text(8, 6.6, "Cohen's D [95% CI]", cex=1.6)
  text(11, 7.6, "bHR", cex=1.8)
  mtext(paste("Adult\n", temp$GeneSymbol[i]), line=-1.5, cex=2.5)
  dev.off()
}

#Interesting - the formatting isn't quite right if there are different numbers of studies included in the output.
i<-2

#This works for Cav1, Sdc4, Mfge8, C1qa, Trhr, Ucp2 (5 datasets)
text(-11, 6.6, "Investigator & Study", cex=1.6)
text(-14, 7.6, "bLR", cex=1.8)
text(8, 6.6, "Cohen's D [95% CI]", cex=1.6)
text(11, 7.6, "bHR", cex=1.8)

#This works for Sox2, Hes5, Cd24, Gpc1 (which have data from 4 studies) - height = 8,:
text(-11, 5.6, "Investigator & Study", cex=1.6)
text(-14, 6.6, "bLR", cex=1.8)
text(8, 5.6, "Cohen's D [95% CI]", cex=1.6)
text(11, 6.6, "bHR", cex=1.8)

#This works for Sox9, Mki67 (which have data from 3 studies) - height = 6,:
text(-11, 4.6, "Investigator & Study", cex=1.6)
text(-14, 5.6, "bLR", cex=1.8)
text(8, 4.6, "Cohen's D [95% CI]", cex=1.6)
text(11, 5.6, "bHR", cex=1.8)


#Exporting forest plots of GOI found in only 3 datasets

temp<-subset(meta, meta$Symbol=="Etv4" | meta$Symbol=="Sox9")

for(i in 1:length(temp$Symbol)){
  
  tiff(paste0(outGOIPlots, "Forest Plot ", 
              paste((temp$Symbol)[i]), ".tiff"), res=300, compression="lzw",
       width = 8.5, height = 8, units = 'in')
  
  effect<-c(temp$`d Affy F4`[i], temp$`d RNAseq F29`[i], temp$`d NC F34`[i], 
            temp$`d F37`[i], temp$`d F43`[i])
  var<-c(temp$`var d Affy F4`[i], temp$`var d RNAseq F29`[i], temp$`var d NC F34`[i], 
         temp$`var d F37`[i], temp$`var d F43`[i])
  
  forest.rma(rma.mv(effect, var), slab=c("MBNI_AffymetrixRae230_F4", "MBNI_RNASeq_F29", 
                                         "Alabama_NimbleGen_F34", "MBNI_RNASeq_F37",
                                         "MBNI_RNASeq_F43"), 
             xlim=c(-18, 15), cex=1.6)
  text(-11, 3.6, "Investigator & Study", cex=1.6)
  text(-14, 4.6, "bLR", cex=1.8)
  mtext(paste("Adult\n", temp$Symbol[i]), line=-1.5, cex=2.5)
  text(8, 3.6, "Cohen's D [95% CI]", cex=1.6)
  text(11, 4.6, "bHR", cex=1.8)
  
  dev.off()
}



### Cell Type Forest Plots


#Exporting cell type forest plots as TIFF files

for(i in 1:10){
  
  tiff(paste0(outCellPlots, "Cell Type Forest Plot ", 
              paste(row.names(CohD_Adult_CellType)[i]), ".tiff"), res=300, 
       compression="lzw", width = 7.5, height = 7.5, units = 'in')
  
  effect<-c(CohD_Adult_CellType$`d Affy F4`[i], CohD_Adult_CellType$`d NC F34`[i], 
            CohD_Adult_CellType$`d F37`[i], CohD_Adult_CellType$`d F43`[i])
  var<-c(CohD_Adult_CellType$`var d Affy F4`[i], CohD_Adult_CellType$`var d NC F34`[i],
         CohD_Adult_CellType$`var d F37`[i], CohD_Adult_CellType$`var d F43`[i])
  
  forest.rma(rma.mv(effect, var), slab=c("MBNI_AffymetrixRae230_F4", "Alabama_NimbleGen_F34",
                                         "MBNI_RNASeq_F37", "MBNI_RNASeq_F43"), 
             xlim=c(-19, 15), cex=1.6)
  text(-12, 4.7, "Investigator & Study", cex=1.6)
  text(-15, 5.5, "bLR", cex=1.8)
  mtext(paste("Adult\n", row.names(CohD_Adult_CellType)[i]), line=-1.5, cex=2.4)
  text(8, 4.7, "Cohen's D [95% CI]", cex=1.6)
  text(11, 5.5, "bHR", cex=1.8)
  
  dev.off()
  
}




######### Volcano Plots


### Generating Volcano Plots to Observe Overall Spread of Expression Values (using estimate)

colnames(metaAnalysisOutputFDR)


############Adult expression colored by FDR and estimate

tiff(paste0(outDir, "VolcanoPlotAdult.tiff"), width = 5, height = 5, 
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))

par(mai=c(1.02, 1,0.9,0.40))
with(metaAnalysisOutputFDR, plot(estimate, -log10(pval), pch=19, main="Overall Expression", 
                                 xlim=c(-3,3), cex.lab=1.8, cex.main=2, cex=0.6))
# Add colored points: red if padj<0.05, blue of estimate>1)
with(subset(metaAnalysisOutputFDR, abs(estimate)>1), points(estimate, -log10(pval), 
                                                            pch=19, col="red", cex=0.6))
with(subset(metaAnalysisOutputFDR, BH< .05 ), points(estimate, -log10(pval), 
                                                     pch=19, col="blue", cex=0.6))
legend(-1.5, 7.6, legend=c("estimate > 1", "FDR < 0.05"), col=c("red", "blue"), pch=19, cex=1.2)

dev.off()



#############Overall Adult expression colored by dataset

tiff(paste0(outDir, "VolcanoPlotAdult_byDatasets.tiff"), width = 5, height = 5, 
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))

with(metaAnalysisOutputFDR, plot(estimate, -log10(pval), pch=19, main="Overall Expression by\n Dataset Number", 
                                 xlim=c(-3,3), cex.lab=1.8, cex.main=2, cex=0.6))

#color by dataset
with(subset(metaAnalysisOutputFDR, datasets > 4), points(estimate, -log10(pval), pch=19, col="blue", cex=0.6))
with(subset(metaAnalysisOutputFDR, datasets < 5 & datasets > 3), points(estimate, -log10(pval), pch=19, col="red", cex=0.6))
with(subset(metaAnalysisOutputFDR, datasets < 4 & datasets > 2), points(estimate, -log10(pval), pch=19, col="green3", cex=0.6))
with(subset(metaAnalysisOutputFDR, datasets < 3), points(estimate, -log10(pval), pch=19, col="purple", cex=0.6))

legend(-.5, 7.5, legend=c("2", "3", "4", "5"), 
       col=c("purple", "green3", "red", "blue"), pch=19, cex=1.2)


dev.off()




############Overall Adult Expression of cell-type specific genes
## Adding meta genes to cell type specific gene matrix in order to view overall expression of only the genes used in the cell type meta-analysis


#Read in list of cell type specific genes
cellInfo <- read.csv("data/CellTypeSpecificGenes_Master3.csv")
colnames(cellInfo)[5]<-"GeneSymbol" #renamae gene column to join

#Join data with cell type info to get only cell type specific genes
temp<-join(metaAnalysisOutputFDR, cellInfo, by="GeneSymbol", type="inner")

#Remove duplicates
temp<-subset(temp, !duplicated(temp$GeneSymbol))



#Output volcano plot
tiff(paste0(outDir, "VolcanoPlot CellTypeAdult.tiff"), width = 5, height = 5, units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.50))
with(temp, plot(estimate, -log10(pval), pch=19, main="Adult Overall Expression\nof Cell Type Specific Genes", 
                xlim=c(-3,3), cex.lab=1.8, cex.main=2, cex=0.6))

# Add colored points: red if padj<0.05, blue of estimate>1)
with(subset(temp, abs(estimate)>1), points(estimate, -log10(pval), 
                                           pch=19, col="red", cex=0.6))
with(subset(temp, BH< .05 ), points(estimate, -log10(pval), 
                                    pch=19, col="blue", cex=0.6))
legend(-1.5, 7.6, legend=c("estimate > 1", "FDR < 0.05"), col=c("red", "blue"), pch=19, cex=1.2)

dev.off()


#############################################

#Megan's notes (from April 20, 2018)

#I wanted to double check the effect of two things on the secondary results (fGSEA, PGEA):
#1) The dependency of the estimate on the number of datasets
#2) The bias in Cohen's d towards larger estimates in datasets with smaller sample sizes.

setwd("~/Documents/Microarray Gen/HRLR/ThesisMetaAnalysisOutput")

#Megan's notes:  Instead of re-loading everything, I just read this back in:
AdultMeta_CohDandPval<-read.csv("AdultMeta_CohDandPval.csv", header=T, stringsAsFactors = F)

str(AdultMeta_CohDandPval)

#First, let's output some basic information about the similarity of the different datasets (note: Isabelle already started this process, I just wanted some pretty graphs that I could put in my presentation or in the Suppl section of the paper):

CohDMatrix_OrderedByGeneration<-cbind(AdultMeta_CohDandPval$d.Affy.F4, AdultMeta_CohDandPval$d.RNAseq.F29, AdultMeta_CohDandPval$d.NC.F34, AdultMeta_CohDandPval$d.F37, AdultMeta_CohDandPval$d.F43)
str(CohDMatrix_OrderedByGeneration)
#Numeric - good.
colnames(CohDMatrix_OrderedByGeneration)<-c("Affy_F4", "RNASeq_F29", "NimbleGen_F34", "RNASeq_F37", "RNA_SeqF43")

write.csv(cor(CohDMatrix_OrderedByGeneration, use="complete.obs"), "AdultCohensD_CorMatrix_CompleteObs.csv")

write.csv(cor(CohDMatrix_OrderedByGeneration, use="pairwise.complete.obs"), "AdultCohensD_CorMatrix_PairwiseCompleteObs.csv")

HedgesGMatrix_OrderedByGeneration<-cbind(AdultMeta_CohDandPval$g.Affy.F4, AdultMeta_CohDandPval$g.RNAseq.F29, AdultMeta_CohDandPval$g.NC.F34, AdultMeta_CohDandPval$g.F37, AdultMeta_CohDandPval$g.F43)

colnames(HedgesGMatrix_OrderedByGeneration)<-c("Affy_F4", "RNASeq_F29", "NimbleGen_F34", "RNASeq_F37", "RNA_SeqF43")

write.csv(cor(HedgesGMatrix_OrderedByGeneration, use="complete.obs"), "AdultHedgesG_CorMatrix_CompleteObs.csv")
write.csv(cor(HedgesGMatrix_OrderedByGeneration, use="pairwise.complete.obs"), "AdultHedgesG_CorMatrix_PairwiseCompleteObs.csv")
#Interesting: For the smallest RNA-Seq dataset (the n=2 dataset), the correlation with the RNASeq_F43 dataset drops from 0.11 (using Cohen's D) to -0.11 using HedgesG.  So the estimates across datasets don't agree as well.  That doesn't necessarily mean that they are worse estimates - it may be that the bias introduced by cohen's d is similar across all of these dinky little datasets, increasing the positive correlation amongst all of them. Still, interesting.


library(lattice)

pdf("ScatterplotMatrix_CohensD_Adult.pdf", height=8, width=8)
splom(data.frame(CohDMatrix_OrderedByGeneration), col=1)
dev.off()
#Slight plus-sign appearance - not nearly as bad as for logFC

pdf("ScatterplotMatrix_HedgesG_Adult.pdf", height=8, width=8)
splom(data.frame(HedgesGMatrix_OrderedByGeneration), col=1)
dev.off()
#Hedges G looks basically the same. Only the small F29 RNA-Seq dataset looks any different, and barely.


#Let's do the same thing for the P14 data:

P14Meta_CohDandPval<-read.csv("p14Meta_cohDandpval.csv", header=T, stringsAsFactors = F)
str(P14Meta_CohDandPval)

P14CohDMatrix_OrderedByGeneration<-cbind(P14Meta_CohDandPval$d.P14, P14Meta_CohDandPval$d.affy, P14Meta_CohDandPval$d.illumina, P14Meta_CohDandPval$d.RNAseq, P14Meta_CohDandPval$d.new.colony.P14)
colnames(P14CohDMatrix_OrderedByGeneration)<-c("Affy_F6", "Affy_F15", "Illumina_F15", "RNASeq_F29", "NimbleGen_F34")

write.csv(cor(P14CohDMatrix_OrderedByGeneration, use="complete.obs"), "P14CohensD_CorMatrix_CompleteObs.csv")
write.csv(cor(P14CohDMatrix_OrderedByGeneration, use="pairwise.complete.obs"), "P14CohensD_CorMatrix_PairwiseCompleteObs.csv")

pdf("ScatterplotMatrix_CohensD_P14.pdf", height=8, width=8)
splom(data.frame(P14CohDMatrix_OrderedByGeneration), col=1)
dev.off()

P7Meta_CohDandPval<-read.csv("P7cohDandPval.csv", header=T, stringsAsFactors = F)
str(P7Meta_CohDandPval)

pdf("Scatterplot_CohensD_P7.pdf", height=4, width=4)
plot(P7Meta_CohDandPval$d.new.colony.P7~P7Meta_CohDandPval$d.P7, ylab="NimbleGen_F34", xlab="Affy_F6")
dev.off()

P21Meta_CohDandPval<-read.csv("P21cohDandPval.csv", header=T, stringsAsFactors = F)
str(P21Meta_CohDandPval)

pdf("Scatterplot_CohensD_P21.pdf", height=4, width=4)
plot(P21Meta_CohDandPval$d.new.colony.P21~P21Meta_CohDandPval$d.P21, ylab="NimbleGen_F34", xlab="Affy_F6")
dev.off()


#Ok, final snoop: I want to output a meta-analysis just based on the findings from the two latest generation RNA-Seq datasets and see how their results compare to the full adult meta-analysis.


metaOutput_JustLateGen<-matrix(0, length(AdultMeta_CohDandPval$GeneSymbol), 6)

for(i in 1:length(AdultMeta_CohDandPval$GeneSymbol)){
  effect<-c(AdultMeta_CohDandPval$d.F37[i], AdultMeta_CohDandPval$d.F43[i])
  var<-c(AdultMeta_CohDandPval$var.d.F37[i], AdultMeta_CohDandPval$var.d.F43[i])
  if(sum(is.na(effect))>0){}else{
    metaOutput_JustLateGen[i, 1]<-rma.mv(effect, var)$b #gives estimate
    metaOutput_JustLateGen[i, 2]<-rma.mv(effect, var)$se #gives standard error
    metaOutput_JustLateGen[i, 3]<-rma.mv(effect, var)$pval #gives pval
    metaOutput_JustLateGen[i, 4]<-rma.mv(effect, var)$ci.lb #gives confidence interval lower bound
    metaOutput_JustLateGen[i, 5]<-rma.mv(effect, var)$ci.ub #gives confidence interval upper bound
  }
  metaOutput_JustLateGen[i, 6] <- sum(2 - sum(is.na(effect)))
}

colnames(metaOutput_JustLateGen)<-c("estimate", "SE", "pval", "CI_lb", "CI_ub", "datasets")

metaOutput_JustLateGenSymbol<-cbind.data.frame(GeneSymbol=AdultMeta_CohDandPval$GeneSymbol, metaOutput_JustLateGen)


sum(metaOutput_JustLateGenSymbol$pval==0) 
#[1] 2950
#Some genes are found in only one dataset and 
#therefor have no values associated with them from the meta-analysis

metaOutput_JustLateGenSymbol<-subset(metaOutput_JustLateGenSymbol, !(metaOutput_JustLateGenSymbol$pval==0))
dim(metaOutput_JustLateGenSymbol)
#[1] 13319     7
#total genes present in both datasets


#write.csv(metaOutputSymbol, paste0(outDir, "AdultMetaAnalysisOutput.csv"))

## Multiple Comparison corrections
#The following code applies two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli) 
#Meta-analysis output with adjusted p-values is then output along with effect size information.

colnames(metaOutput_JustLateGenSymbol)

library(multtest)
tempPvalAdjMeta<-mt.rawp2adjp(metaOutput_JustLateGenSymbol$pval, proc=c("BH", "BY"))
metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]

#adjusted pvalue object is in same orientation as metaoutput so can simply be binded together

dim(metaOutput_JustLateGenSymbol)
dim(metaPvalAdj)

metaOutput_JustLateGen_FDR<-cbind.data.frame(metaOutput_JustLateGenSymbol, metaPvalAdj)
dim(metaOutput_JustLateGen_FDR)
str(metaOutput_JustLateGen_FDR)

write.csv(metaOutput_JustLateGen_FDR, "metaOutput_JustLateGen_FDR.csv")
#More significant genes, order seems relatively similar, at least at first glance.

#Let's compare these results with those from the full adult meta-analysis:
library(plyr)

AdultMetaAnalysisOutputFDR<-read.csv("AdultMetaAnalysisOutputFDR_forR.csv", header=T, stringsAsFactors = F)
str(AdultMetaAnalysisOutputFDR)

colnames(AdultMetaAnalysisOutputFDR)[-1]<-paste("Full", colnames(AdultMetaAnalysisOutputFDR)[-1], sep="_")

colnames(metaOutput_JustLateGen_FDR)[-1]<-paste("LateGen", colnames(metaOutput_JustLateGen_FDR)[-1], sep="_")

AdultMetaAnalysisOutputFDR_vs_JustLateGen<-join(AdultMetaAnalysisOutputFDR, metaOutput_JustLateGen_FDR, by="GeneSymbol", type="inner")

pdf("Estimate_FullvsOnlyLateGenMetaAnalysis.pdf", width=8, height=8)

par(mfrow=c(2,2))

plot(AdultMetaAnalysisOutputFDR_vs_JustLateGen$Full_estimate~AdultMetaAnalysisOutputFDR_vs_JustLateGen$LateGen_estimate, xlab="Estimate: Only RNA-Seq Data from F37 and F43", ylab="Estimate: Full Meta-Analysis")

plot(AdultMetaAnalysisOutputFDR_vs_JustLateGen$Full_estimate~AdultMetaAnalysisOutputFDR_vs_JustLateGen$LateGen_estimate, col=AdultMetaAnalysisOutputFDR_vs_JustLateGen$Full_datasets, xlab="Estimate: Only RNA-Seq Data from F37 and F43", ylab="Estimate: Full Meta-Analysis")

plot(y=(rank(AdultMetaAnalysisOutputFDR_vs_JustLateGen$Full_estimate)/13319), x=(rank(AdultMetaAnalysisOutputFDR_vs_JustLateGen$LateGen_estimate)/13319), xlab="Percentile Rank: Only RNA-Seq Data from F37 and F43", ylab="Percentile Rank: Full Meta-Analysis")

plot(y=(rank(AdultMetaAnalysisOutputFDR_vs_JustLateGen$Full_estimate)/13319), x=(rank(AdultMetaAnalysisOutputFDR_vs_JustLateGen$LateGen_estimate)/13319), xlab="Percentile Rank: Only RNA-Seq Data from F37 and F43", ylab="Percentile Rank: Full Meta-Analysis", col=AdultMetaAnalysisOutputFDR_vs_JustLateGen$Full_datasets)

dev.off()


#How about a Venn Diagram of top genes?

library(VennDiagram)

AdultMetaAnalysisOutputFDR_vs_JustLateGen$GeneSymbol[AdultMetaAnalysisOutputFDR_vs_JustLateGen$Full_BH<0.05]
length(AdultMetaAnalysisOutputFDR_vs_JustLateGen$GeneSymbol[AdultMetaAnalysisOutputFDR_vs_JustLateGen$Full_BH<0.05])
#[1] 72

AdultMetaAnalysisOutputFDR_vs_JustLateGen$GeneSymbol[AdultMetaAnalysisOutputFDR_vs_JustLateGen$LateGen_BH<0.05]
length(AdultMetaAnalysisOutputFDR_vs_JustLateGen$GeneSymbol[AdultMetaAnalysisOutputFDR_vs_JustLateGen$LateGen_BH<0.05])
#[1] 494

sum(AdultMetaAnalysisOutputFDR_vs_JustLateGen$GeneSymbol[AdultMetaAnalysisOutputFDR_vs_JustLateGen$Full_BH<0.05]%in%AdultMetaAnalysisOutputFDR_vs_JustLateGen$GeneSymbol[AdultMetaAnalysisOutputFDR_vs_JustLateGen$LateGen_BH<0.05])
#[1] 68

72-68
#[1] 4
494-68
#[1] 426


grid.newpage()
draw.pairwise.venn(72, 494, 68, category = c("Full Meta-Analysis", "Only RNA-Seq Data from F37 and F43"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))

AdultMetaAnalysisOutputFDR_vs_JustLateGen$GeneSymbol[AdultMetaAnalysisOutputFDR_vs_JustLateGen$Full_BH<0.10]
length(AdultMetaAnalysisOutputFDR_vs_JustLateGen$GeneSymbol[AdultMetaAnalysisOutputFDR_vs_JustLateGen$Full_BH<0.10])  
#[1] 185

AdultMetaAnalysisOutputFDR_vs_JustLateGen$GeneSymbol[AdultMetaAnalysisOutputFDR_vs_JustLateGen$LateGen_BH<0.10]
length(AdultMetaAnalysisOutputFDR_vs_JustLateGen$GeneSymbol[AdultMetaAnalysisOutputFDR_vs_JustLateGen$LateGen_BH<0.10])
#[1] 892

sum(AdultMetaAnalysisOutputFDR_vs_JustLateGen$GeneSymbol[AdultMetaAnalysisOutputFDR_vs_JustLateGen$Full_BH<0.10]%in%AdultMetaAnalysisOutputFDR_vs_JustLateGen$GeneSymbol[AdultMetaAnalysisOutputFDR_vs_JustLateGen$LateGen_BH<0.10])
#[1] 171

grid.newpage()
draw.pairwise.venn(185, 892, 171, category = c("Full Meta-Analysis", "Only RNA-Seq Data from F37 and F43"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))

par(mfrow=c(1,1))

pdf("VennDiagram_FullvsOnlyLateGenMetaAnalysis_FDR05.pdf", width=4, height=4)
grid.newpage()
draw.pairwise.venn(72, 494, 68, category = c("Full Meta-Analysis", "Only RNA-Seq Data from F37 and F43"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))
dev.off()

pdf("VennDiagram_FullvsOnlyLateGenMetaAnalysis_FDR10.pdf", width=4, height=4)
grid.newpage()
draw.pairwise.venn(185, 892, 171, category = c("Full Meta-Analysis", "Only RNA-Seq Data from F37 and F43"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))
dev.off()


#####################################

#Ok, snooping once again: It turns out the rma.mv function was defaulting to a fixed-effects model (the inverse variance method - even though all of the output says random effects!).  Let's see what happens if we run either an actual random effects model.


metaOutput_RE<-matrix(0, length(AdultMeta_CohDandPval$GeneSymbol), 6)

for(i in 1:length(AdultMeta_CohDandPval$GeneSymbol)){
  effect<-c(AdultMeta_CohDandPval$d.Affy.F4[i], AdultMeta_CohDandPval$d.RNAseq.F29[i], AdultMeta_CohDandPval$d.NC.F34[i], AdultMeta_CohDandPval$d.F37[i], AdultMeta_CohDandPval$d.F43[i])
  var<-c(AdultMeta_CohDandPval$var.d.Affy.F4[i], AdultMeta_CohDandPval$var.d.RNAseq.F29[i], AdultMeta_CohDandPval$var.d.NC.F34[i],  AdultMeta_CohDandPval$var.d.F37[i], AdultMeta_CohDandPval$var.d.F43[i])
  if(sum(is.na(effect))>3){}else{
    temp<-rma(effect, var)
    metaOutput_RE[i, 1]<-temp$b #gives estimate
    metaOutput_RE[i, 2]<-temp$se #gives standard error
    metaOutput_RE[i, 3]<-temp$pval #gives pval
    metaOutput_RE[i, 4]<-temp$ci.lb #gives confidence interval lower bound
    metaOutput_RE[i, 5]<-temp$ci.ub #gives confidence interval upper bound
  }
  metaOutput_RE[i, 6] <- sum(5 - sum(is.na(effect)))
}

colnames(metaOutput_RE)<-c("estimate", "SE", "pval", "CI_lb", "CI_ub", "datasets")

metaOutput_REGenSymbol<-cbind.data.frame(GeneSymbol=AdultMeta_CohDandPval$GeneSymbol, metaOutput_RE)


sum(metaOutput_REGenSymbol$pval==0) 
#[1] 0
#No genes are found in only one dataset and 
#therefor have no values associated with them from the meta-analysis
#Why is that?  There should be 2950. Ah - maybe that is just for the late gene only output.
hist(metaOutput_RE[, 6])
hist(AdultMetaAnalysisOutputFDR$Full_datasets)
length(AdultMetaAnalysisOutputFDR[,1])


#write.csv(metaOutputSymbol, paste0(outDir, "AdultMetaAnalysisOutput.csv"))

## Multiple Comparison corrections
#The following code applies two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli) 
#Meta-analysis output with adjusted p-values is then output along with effect size information.

colnames(metaOutput_REGenSymbol)

library(multtest)
tempPvalAdjMeta<-mt.rawp2adjp(metaOutput_REGenSymbol$pval, proc=c("BH", "BY"))
metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]

#adjusted pvalue object is in same orientation as metaoutput so can simply be binded together

dim(metaOutput_REGenSymbol)
dim(metaPvalAdj)

metaOutput_RE_FDR<-cbind.data.frame(metaOutput_REGenSymbol, metaPvalAdj)
dim(metaOutput_RE_FDR)
str(metaOutput_RE_FDR)

write.csv(metaOutput_RE_FDR, "metaOutput_RE_FDR.csv")
#Looks like the same order as before, but fewer of the genes survive FDR.
#Ok, it's actually not quite the same order - and it looks like the genes that have much greater change with generation are getting bumped down.  Maybe not so ideal.

pdf("PlotRandomVsFixedEffectsEstimate.pdf")
plot(metaOutput_RE_FDR$estimate~AdultMeta_CohDandPval$estimate, xlab="Fixed Effects Estimate", ylab="Random Effects Estimate")
dev.off()
#Basically the same estimate from both models

pdf("PlotRandomVsFixedEffects_LogPval.pdf")
plot(log(metaOutput_RE_FDR$pval, 10)~log(AdultMeta_CohDandPval$pval,10), col=(metaOutput_RE_FDR$datasets<3)+1, xlab="Fixed Effects Log(10)Pval", ylab="Random Effects Log(10)Pval")
dev.off()
#Many genes have similar p-values, some have less significant p-values in the RE model. Depends on the amount of heterogeneity present in the studies.

pdf("PlotRandomVsFixedEffects_Rank.pdf")
plot(rank(metaOutput_RE_FDR$pval)~rank(AdultMeta_CohDandPval$pval), col=(metaOutput_RE_FDR$datasets<3)+1, xlab="Fixed Effects Pval Rank", ylab="Random Effects Pval Rank")
dev.off()


#So let's try including generation as a variable (although note all of the disclaimers associated with that- e.g., platform and generation go hand in hand, and that sort of model will use 2 degrees of freedom - maybe only apply to genes with at least 5 studies?)


Generation<-c(4, 29, 34, 37, 43)
GenerationCentered<-Generation-43

metaOutput_wGen<-matrix(0, length(AdultMeta_CohDandPval$GeneSymbol), 11)

for(i in 1:length(AdultMeta_CohDandPval$GeneSymbol)){
  effect<-c(AdultMeta_CohDandPval$d.Affy.F4[i], AdultMeta_CohDandPval$d.RNAseq.F29[i], AdultMeta_CohDandPval$d.NC.F34[i], AdultMeta_CohDandPval$d.F37[i], AdultMeta_CohDandPval$d.F43[i])
  var<-c(AdultMeta_CohDandPval$var.d.Affy.F4[i], AdultMeta_CohDandPval$var.d.RNAseq.F29[i], AdultMeta_CohDandPval$var.d.NC.F34[i],  AdultMeta_CohDandPval$var.d.F37[i], AdultMeta_CohDandPval$var.d.F43[i])
 # if(sum(is.na(effect))>2){}else{
    if(sum(is.na(effect))>1){}else{
    temp<-rma(effect, var, mods=~GenerationCentered)
    metaOutput_wGen[i, 1]<-temp$b[1] #gives estimate
    metaOutput_wGen[i, 2]<-temp$se[1] #gives standard error
    metaOutput_wGen[i, 3]<-temp$pval[1] #gives pval
    metaOutput_wGen[i, 4]<-temp$ci.lb[1] #gives confidence interval lower bound
    metaOutput_wGen[i, 5]<-temp$ci.ub[1] #gives confidence interval upper bound
    metaOutput_wGen[i, 6]<-temp$b[2] #gives estimate
    metaOutput_wGen[i, 7]<-temp$se[2] #gives standard error
    metaOutput_wGen[i, 8]<-temp$pval[2] #gives pval
    metaOutput_wGen[i, 9]<-temp$ci.lb[2] #gives confidence interval lower bound
    metaOutput_wGen[i, 10]<-temp$ci.ub[2] #gives confidence interval upper bound
    
  }
  metaOutput_wGen[i, 11] <- sum(5 - sum(is.na(effect)))
}

colnames(metaOutput_wGen)<-c("HRLR_estimate", "HRLR_SE", "HRLR_pval", "HRLR_CI_lb", "HRLR_CI_ub", "Gen_estimate", "Gen_SE", "Gen_pval", "Gen_CI_lb", "Gen_CI_ub", "datasets")

metaOutput_wGen_GenSymbol<-cbind.data.frame(GeneSymbol=AdultMeta_CohDandPval$GeneSymbol, metaOutput_wGen)


sum(metaOutput_wGen_GenSymbol$HRLR_pval==0) 
#[1] 2891 - when using min 3 datasets
#[1] 7106 - when using min 4 datasets

hist(metaOutput_wGen[, 11])

metaOutput_wGen_GenSymbol<-metaOutput_wGen_GenSymbol[metaOutput_wGen_GenSymbol$HRLR_pval!=0,]
dim(metaOutput_wGen_GenSymbol)
#[1] 13378    12 - when using min 3 datasets
#[1] 9163   12 - when using min 4 datasets

#write.csv(metaOutputSymbol, paste0(outDir, "AdultMetaAnalysisOutput.csv"))

## Multiple Comparison corrections
#The following code applies two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli) 
#Meta-analysis output with adjusted p-values is then output along with effect size information.

colnames(metaOutput_wGen_GenSymbol)

library(multtest)
tempPvalAdjMeta<-mt.rawp2adjp(metaOutput_wGen_GenSymbol$HRLR_pval, proc=c("BH", "BY"))
metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]

#adjusted pvalue object is in same orientation as metaoutput so can simply be binded together

dim(metaOutput_wGen_GenSymbol)
dim(metaPvalAdj)

metaOutput_wGen_FDR<-cbind.data.frame(metaOutput_wGen_GenSymbol, metaPvalAdj)
dim(metaOutput_wGen_FDR)
str(metaOutput_wGen_FDR)
colnames(metaOutput_wGen_FDR)[c(13:15)]<-c("HRLR_rawp", "HRLR_BH", "HRLR_BY")

tempPvalAdjMeta<-mt.rawp2adjp(metaOutput_wGen_GenSymbol$Gen_pval, proc=c("BH", "BY"))
metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]

metaOutput_wGen_FDR<-cbind.data.frame(metaOutput_wGen_FDR, metaPvalAdj)
dim(metaOutput_wGen_FDR)
str(metaOutput_wGen_FDR)
colnames(metaOutput_wGen_FDR)[c(16:18)]<-c("Gen_rawp", "Gen_BH", "Gen_BY")

write.csv(metaOutput_wGen_FDR, "metaOutput_wGen_FDR.csv")

boxplot(abs(metaOutput_wGen_FDR$HRLR_estimate)~metaOutput_wGen_FDR$datasets)
summary.lm(lm(abs(metaOutput_wGen_FDR$HRLR_estimate)~metaOutput_wGen_FDR$datasets))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                   1.001644   0.023747   42.18   <2e-16 ***
#   metaOutput_wGen_FDR$datasets -0.100002   0.005804  -17.23   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.538 on 13376 degrees of freedom
# Multiple R-squared:  0.02171,	Adjusted R-squared:  0.02164 
# F-statistic: 296.9 on 1 and 13376 DF,  p-value: < 2.2e-16
#looks like the relationship between estimate and #datasets is still there, despite trying to control for generation

#For comparison:
boxplot(abs(AdultMeta_CohDandPval$estimate)~AdultMeta_CohDandPval$datasets)

#I would like to see what happens when I only use the genes in 4 to 5 datasets (so better estimate of slope for Gen)
#Still no genes with effects of generation that surpass FDR<0.30.


####################################################

#And yet one more analysis:  What if I run a model with all ages included, but include age as a term in the model? 
#The coding for this may be funky, because P7 and P21 only have data from a few genes.

Age<-c("P7", "P7", "P14", "P14", "P14", "P14", "P14", "P21", "P21", "Adult", "Adult", "Adult", "Adult", "Adult")
table(Age)

#Oh wait - that's right. Our Cohen's d output is really not set up for this right now (because the output is by age and doesn't include any genes that weren't found in at least 2 studies for each age.)  I'm going to skip this and come back to it if a reviewer (or anyone else) is interested. 
