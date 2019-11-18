library(car)
library(plyr)


setwd("~/Phenotype Project/ibirt")

######################

#Reading in Peter's data:

PeterHR_LRdata <- read.csv("file:///C:/Users/Izzy/Documents/Phenotype Project/PeterBlandino_BasalData/Original Files/LowerCase_Akil5_geneexp_filteredandnormalized_log2fpm.csv")
colnames(PeterHR_LRdata)

#Megan's code for reading in Peter's data:
setwd("~/Documents/Microarray Gen/HRLR/HRLR_Studies/NIH_basalHRLR_RNAseq")

PeterHR_LRdata <- read.csv("LowerCase_Akil5_geneexp_filteredandnormalized_log2fpm.csv", header=T, stringsAsFactors = F)
colnames(PeterHR_LRdata)
# [1] "Gene"   "X10HR"  "X22HR"  "X26HR"  "X2HR"   "X34HR"  "X38HR"  "X149IR" "X158IR" "X213IR" "X248IR" "X261IR" "X293IR" "X16LR"  "X20LR" [16] "X32LR"  "X44LR"  "X48LR"  "X8LR"  
head(PeterHR_LRdata)
# Gene     X10HR      X22HR     X26HR      X2HR      X34HR     X38HR    X149IR    X158IR     X213IR    X248IR    X261IR    X293IR
# 1     A1bg -3.316754 -1.8308858 -6.796308 -1.877379 -2.3757866 -2.811070 -4.192394 -2.010690 -0.9914813 -3.955577 -3.383915 -1.454892
# 2      A2m  4.438876  4.6586755  4.076367  4.447728  4.6967480  4.708129  4.064050  4.814896  4.9379770  4.054531  4.887548  4.317394
# 3  A3galt2  1.314221  0.9163481  1.092436  1.502763  1.0660509  1.453633  1.193037  1.627911  1.6632514  1.716848  1.347027  1.294457
# 4 Aa926063 -2.579788 -2.5227635 -1.510905 -3.157487 -0.7908241 -3.659067 -2.607432 -2.531522 -0.9914813 -1.885188 -2.466377 -3.142948
# 5     Aaas  4.233566  4.1288402  4.180257  4.147597  4.1404387  4.234538  4.303461  4.102895  4.2221812  4.193154  4.030183  4.118625
# 6     Aacs  4.551555  4.7901195  5.218759  4.437104  4.9717609  4.885511  5.229250  4.815666  4.7927113  4.476547  4.881310  5.100226
## No need to check data quality as this data has already been processed
#Moving onto normalization - keeping IRs for behavior analyses

######
## Reading in Peter's sample descriptions/behavior scores

sampleInfo <- read.csv("file:///C:/Users/Izzy/Documents/Phenotype Project/PeterBlandino_BasalData/Original Files/PeterHRLR_Akil5_sample_description.csv")
colnames(sampleInfo)

#Megan's code for reading in behavioral data - note, I made sure that the behavioral data was in the same order as the gene expression data.
setwd("~/Documents/Microarray Gen/HRLR/HRLR_Studies/NIH_basalHRLR_RNAseq/akil-5")

sampleInfo <- read.csv("Akil5_sample_description.csv", header=T, stringsAsFactors = F)
colnames(sampleInfo)
# [1] "Rat_ID"             "FamilyID"           "Generation"         "Lineage"            "RatID"              "TotalLocoScore"    
# [7] "Boli"               "DistanceTraveled"   "TimeImmobile"       "PercentTimeOpenArm" "Sex"    

#################

#Reading in Cigdem's data:

setwd("~/Documents/Microarray Gen/HRLR/HRLR_Studies/Cigdem_HRLR_RNAseq")
Cigdem_HRLRdata<-read.csv("CigdemNew_Vehicle.csv", header=T, stringsAsFactors = F)
colnames(Cigdem_HRLRdata)
#[1] "X"         "HR_CONT_1" "HR_CONT_2" "HR_CONT_3" "HR_CONT_4" "HR_CONT_5" "LR_CONT_1" "LR_CONT_2" "LR_CONT_3" "LR_CONT_4" "LR_CONT_5"
head(Cigdem_HRLRdata)
# X HR_CONT_1 HR_CONT_2 HR_CONT_3 HR_CONT_4 HR_CONT_5 LR_CONT_1 LR_CONT_2 LR_CONT_3 LR_CONT_4 LR_CONT_5
# 1 ENSRNOG00000029897  1.984111  1.611828  1.911635  1.879735  1.978242  1.592064  1.872764  2.052356  2.135435  2.247886
# 2 ENSRNOG00000014303  8.080741  8.114191  8.115799  8.109010  8.074836  8.069525  8.042232  8.146873  8.127431  8.047361
# 3 ENSRNOG00000014330  6.626785  6.711768  6.653282  6.678040  6.657021  6.686488  6.656599  6.742801  6.714538  6.740110
# 4 ENSRNOG00000049505  2.637977  2.645080  2.541030  2.665524  2.627412  2.746949  2.674064  2.697691  2.726000  2.758554
# 5 ENSRNOG00000014916  5.624220  5.667279  5.690249  5.689228  5.776825  5.719697  5.760098  5.726541  5.717119  5.823067
# 6 ENSRNOG00000014996  3.044383  3.196790  3.199910  3.232999  3.202155  3.065995  3.249105  3.086303  2.972160  3.228257

Cigdem_Annotation<-read.csv("GeneAnnotation_forNewData.csv", header=T, stringsAsFactors = F)
colnames(Cigdem_Annotation)
#[1] "ENSEMBLID"  "GeneSymbol" "GeneType" 
head(Cigdem_Annotation)
# ENSEMBLID     GeneSymbol       GeneType
# 1 ENSRNOG00000000007           Gad1 protein_coding
# 2 ENSRNOG00000000008           Alx4 protein_coding
# 3 ENSRNOG00000000010          Cbln1 protein_coding
# 4 ENSRNOG00000000012          Tcf15 protein_coding

Cigdem_sampleInfo <- read.csv("sampledescription_justVeh.csv", header=T, stringsAsFactors = F)
colnames(Cigdem_sampleInfo)
#[1] "Sample_ID"   "Core_ID"     "Sample_Desc" "Rat_Type"    "Treatment"   "SI_Score"    "Group"  
head(Cigdem_sampleInfo)
# Sample_ID Core_ID Sample_Desc Rat_Type Treatment SI_Score   Group
# 1        53   53816   HR_CONT_1       HR      CONT    14.00 HR_CONT
# 2        54   53817   HR_CONT_2       HR      CONT    16.67 HR_CONT
# 3       114   53825   HR_CONT_3       HR      CONT    17.67 HR_CONT
# 4       115   53826   HR_CONT_4       HR      CONT    24.33 HR_CONT
# 5       132   53829   HR_CONT_5       HR      CONT    24.00 HR_CONT
# 6        11   53807   LR_CONT_1       LR      CONT    12.33 LR_CONT
#Note: I already made sure these were in the same order as the gene expression data.


#################  Normalizing data

##Checking for duplicate gene symbols
sum(duplicated(PeterHR_LRdata$Gene))
#[1] 2
#only two. Will still average across gene symbols so there are no duplicates.

#Averaging by gene symbol
AveragebygenesymbolPeterHRLR<-matrix(0, length(names(table(PeterHR_LRdata[,1]))), 18)
for(i in 2:19){#the forloop runs over every column of data and calculates the mean by gene symbol:
  AveragebygenesymbolPeterHRLR[,i-1]<-tapply(PeterHR_LRdata[,i], PeterHR_LRdata[,1], function(y) mean(y))
}

colnames(AveragebygenesymbolPeterHRLR)<-colnames(PeterHR_LRdata[2:19]) #Assigning subject names to columns
row.names(AveragebygenesymbolPeterHRLR)<-(names(table(PeterHR_LRdata[,1]))) #Assigning gene names as row names, however they must be in alphabetical order now because averaging caused them the values to be rearranged
dim(AveragebygenesymbolPeterHRLR) #checking dimensions
#[1] 20907    18

#Megan's notes: I cut out the scaling and centering step, because it is unnecessary and makes the figures harder to understand:
### Scaling and Centering to convert to zscores
# AveragebygenesymbolPeterHRLR_zscore<-t(scale(t(AveragebygenesymbolPeterHRLR), center=TRUE, scale=TRUE))
# str(AveragebygenesymbolPeterHRLR_zscore)
# dim(AveragebygenesymbolPeterHRLR_zscore)
#[1] 20907    18

sum(duplicated(Cigdem_Annotation$GeneSymbol))
#[1] 249
sum(is.na(Cigdem_Annotation$GeneSymbol))
#[1] 0

AveragebygenesymbolCigdemHRLR<-matrix(0, length(names(table(Cigdem_Annotation$GeneSymbol))), 10)
for(i in 2:11){#the forloop runs over every column of data and calculates the mean by gene symbol:
  AveragebygenesymbolCigdemHRLR[,i-1]<-tapply(Cigdem_HRLRdata[,i], Cigdem_Annotation$GeneSymbol, function(y) mean(y))
}

colnames(AveragebygenesymbolCigdemHRLR)<-colnames(Cigdem_HRLRdata[2:11]) #Assigning subject names to columns
row.names(AveragebygenesymbolCigdemHRLR)<-(names(table(Cigdem_Annotation$GeneSymbol))) #Assigning gene names as row names, however they must be in alphabetical order now because averaging caused them the values to be rearranged
dim(AveragebygenesymbolCigdemHRLR) #checking dimensions
#[1] 16394    10

######################

#Creating vector of genes of interest
genesOfInterest <- c("C1qa", "C1qb", "C1qc", "Mfge8", "Ucp2", "Cd4", "Bmp4", "Ncan", "Fxyd7", "Trhr", "Sox9", "Cav1", "Trh", "Mki67", "Uhrf1", "Cav1", "Tek", "Sp3", "C2cd3", "Apba2", "Etv4", "Slc27a1", "Rltpr", "Adamts2", "Fbxo31", "Bdnf", "Fgf2", "Fgf9", "Fgfr3", "Fgfr1", "Itgb4")

geneT <- t(AveragebygenesymbolPeterHRLR)

#pulling out column numbers associated with GOIs
goiCols <- which(colnames(geneT)%in%genesOfInterest==TRUE)

#Subsetting GOIs data
goiData <- cbind.data.frame(geneT[,goiCols])

#Adding column of rat IDs
goiData <- cbind.data.frame(Rat_ID=row.names(goiData), goiData, stringsAsFactors=F)
str(goiData)

library(plyr)
#Joining info and gene expression together
#locoData <- join(sampleInfo, goiData, by="Rat_ID", type="inner")
#Megan's notes: I already made sure these are in the same order so instead:
locoData<-cbind.data.frame(sampleInfo, goiData[,-1], stringsAsFactors=F )

dim(locoData)
#[1] 18 41



## Observing variable effects on gene expression

colnames(locoData)
#[1] "Rat_ID"             "FamilyID"           "Generation"         "Lineage"           
#[5] "RatID"              "TotalLocoScore"     "Boli"               "DistanceTraveled"  
#[9] "TimeImmobile"       "PercentTimeOpenArm" "Sex"                "Bmp4"              
#[13] "C1qa"               "C1qb"               "C1qc"               "Cd4"               
#[17] "Fxyd7"              "Mfge8"              "Ncan"               "Ucp2" 


variableEffects <- matrix(0, ncol=4, nrow=9)


for(i in 12:20){
    
  variableEffects[i-11,1]<-colnames(locoData[i])
  
  variableEffects[i-11,2]<-Anova(lm(locoData[,i]~TotalLocoScore,
                                                  data=locoData),
                                               type=3)$"Pr(>F)"[2]
  
  variableEffects[i-11,3]<-Anova(lm(locoData[,i]~Boli,
                                    data=locoData),
                                 type=3)$"Pr(>F)"[2]
  
  variableEffects[i-11,4]<-Anova(lm(locoData[,i]~PercentTimeOpenArm,
                                    data=locoData),
                                 type=3)$"Pr(>F)"[2]
               }
            
colnames(variableEffects) <- c("Gene", "TotalLocoScore", "Boli", "PercentTimeOpenArm")    


write.csv(variableEffects, "output/Blandino RNAseq F37/BehaviorData_ANOVA3_pvals.csv")



Anova(lm(locoData$PercentTimeOpenArm~locoData$Lineage, data=locoData), type=3)
#Response: locoData$PercentTimeOpenArm
#Sum Sq Df F value    Pr(>F)    
#(Intercept)      7549.4  1 43.2558 8.779e-06 ***
#  locoData$Lineage 2345.7  2  6.7202  0.008245 ** 
#  Residuals        2618.0 15  


Anova(lm(locoData$TotalLocoScore~locoData$Lineage, data=locoData), type=3)
#Response: locoData$TotalLocoScore
#Sum Sq Df F value    Pr(>F)    
#(Intercept)      33257313  1 1627.68 < 2.2e-16 ***
#  locoData$Lineage 18162978  2  444.47 4.463e-14 ***
#  Residuals          306485 15 


Anova(lm(locoData$Lineage~locoData$Boli, data=locoData), type=3)
#Response: locoData$Boli
#Sum Sq Df F value   Pr(>F)   
#(Intercept)       2.667  1  1.2435 0.282341   
#locoData$Lineage 27.444  2  6.3990 0.009786 **
#  Residuals        32.167 15   



#Percent time open arm

out<-c(
  capture.output(
    for (i in 13:25){
      print(paste(colnames(locoData)[i], "vs. PercentTimeOpenArm", sep="  "))
      print(summary.lm(lm(locoData[,i]~locoData$PercentTimeOpenArm, data=locoData)))
      print(" ")
    }
  )
)
cat(out, file="lm of PercentTimeOpenArm_F37Data.txt", sep="\n", append=FALSE)
rm(out)


#total loco scores


out<-c(
  capture.output(
    for (i in 13:25){
      print(paste(colnames(locoData)[i], "vs. TotalLocoScore", sep="  "))
      print(summary.lm(lm(locoData[,i]~locoData$TotalLocoScore, data=locoData)))
      print(" ")
    }
  )
)
cat(out, file="lm of TotalLocoScore_F37Data.txt", sep="\n", append=FALSE)
rm(out)


#LM of boli effects


out<-c(
  capture.output(
    for (i in 13:25){
      print(paste(colnames(locoData)[i], "vs. Boli", sep="  "))
      print(summary.lm(lm(locoData[,i]~locoData$Boli, data=locoData)))
      print(" ")
    }
  )
)
cat(out, file="lm of Boli_F37Data.txt", sep="\n", append=FALSE)
rm(out)





#Examining coefficients for Mfge8
summary(lm(locoData[,18]~TotalLocoScore+Boli+PercentTimeOpenArm,
   data=locoData))


###########################

##  Creating Plots for Effects of Coefficients on Gene Expression


# Peter's data:

#Creating color vector for plots
locoData$Lineage
#[1] HR LR HR LR LR HR HR LR HR HR LR LR IR IR IR IR IR IR

colorPhen<-locoData$Lineage
colorPhen <-gsub("HR", "green3", colorPhen) #Making HRs green
colorPhen <-gsub("LR", "red", colorPhen) #Making LRs red
colorPhen <-gsub("IR", "goldenrod1", colorPhen) #Making IRs gold

#Create folder for boxplots

setwd("~/Documents/Microarray Gen/HRLR/HRLR_Studies/NIH_basalHRLR_RNAseq/GeneBehav_Plots")

#This is me cheating a little for the automatic xlab names:
BehaviorVariableNames<-vector(length=10, mode="character")
BehaviorVariableNames[6]<-"Locomotor Score"
BehaviorVariableNames[7]<-"EPM: Fecal Boli"
BehaviorVariableNames[10]<-"% Time in the Open Arm"

for(i in c(12:41)){
  for(j in c(6,7,10)){
  tiff(paste0(colnames(locoData[i]), 
              "vs", colnames(locoData[j]), "Plot.tiff", sep=" "), res=300, compression="lzw",
              width = 6, height = 6, units = 'in')
    par(mai=c(.9,1.05,0.95,0.4))
    
    plot(locoData[,i]~locoData[,j], data=locoData, col="black", bg=colorPhen, pch=21, cex=1.5, cex.axis=1.5,
       xlab=BehaviorVariableNames[j], cex.lab=2, cex.main=2.2, ylab=paste(colnames(locoData[i]),": Log(2) FPM", sep=""))
  
  temp <- lm(locoData[,i]~locoData[,j])
  abline(temp, lwd=2)
  mtext(paste("R-squared=", signif(summary.lm(temp)$r.squared, 2), ", P-value=", signif(summary.lm(temp)$coefficients[2,4], 2)), cex=1.5, line=1)
  if(temp$coefficients[2] > 0){
    legend("bottomright", legend=c("HR", "IR", "LR"), col=c("green3", "goldenrod1", "red"), 
           pch=19, cex=1.3)
  }else{legend("topright", legend=c("HR", "IR", "LR"), col=c("green3", "goldenrod1", "red"), 
              pch=19, cex=1.3)}
  dev.off()
  }
}


boxplot(locoData$PercentTimeOpenArm~locoData$Lineage)

#Creating boxplot of EPM by phenotype


tiff(paste0("output/Peter Blandino Basal data/Behavior Plots/EPM vs Phenotype.tiff", sep=" "), 
     res=300, compression="lzw",
     width = 4, height = 4.5, units = 'in')
par(mai=c(.8,.6,0.95,0.4))
    
boxplot(locoData$PercentTimeOpenArm~locoData$Lineage, data=locoData, 
        col=c("green3", "goldenrod1", "red"), cex=1, pch=19, 
        main="Percent Time in Open Arm\n by Phenotype", cex.main=1.6)
    
dev.off()


#Updated by MH to match other plots:

tiff(paste0("EPM vs Phenotype.tiff", sep=" "), 
     res=300, compression="lzw",
     width = 4, height = 4.5, units = 'in')
par(mai=c(0.8,1,0.95,0.4))

boxplot(locoData$PercentTimeOpenArm~locoData$Lineage, data=locoData, 
        col=c("green3", "goldenrod1", "red"), cex=1, pch=19, ylab="% Time in the Open Arm", cex.lab=1.3)

dev.off()

tiff(paste0("Locomotor vs Phenotype.tiff", sep=" "), 
     res=300, compression="lzw",
     width = 4, height = 4.5, units = 'in')
par(mai=c(0.8,1,0.95,0.4))

boxplot(locoData$TotalLocoScore~locoData$Lineage, data=locoData, 
        col=c("green3", "goldenrod1", "red"), cex=1, pch=19, ylab="Locomotor Score", cex.lab=1.3)

dev.off()


tiff(paste0("Boli vs Phenotype.tiff", sep=" "), 
     res=300, compression="lzw",
     width = 4, height = 4.5, units = 'in')
par(mai=c(0.8,1,0.95,0.4))

boxplot(locoData$Boli~locoData$Lineage, data=locoData, 
        col=c("green3", "goldenrod1", "red"), cex=1, pch=19, ylab="EPM: Fecal Boli", cex.lab=1.3)

dev.off()
#That's lame


#Outputting all behavioral x gene correlations:


Peter_GeneBehaviorCorrelations<-matrix(0, nrow(AveragebygenesymbolPeterHRLR), ncol=9)
row.names(Peter_GeneBehaviorCorrelations)<-row.names(AveragebygenesymbolPeterHRLR)
colnames(Peter_GeneBehaviorCorrelations)<-c("Locomotor_Beta", "Locomotor_Rsquared", "Locomotor_Pval", "EPMBoli_Beta", "EPMBoli_Rsquared", "EPMBoli_Pval", "PercentOpenArm_Beta", "PercentOpenArm_Rsquared", "PercentOpenArm_Pval")

for(i in c(1:nrow(AveragebygenesymbolPeterHRLR))){
  temp <- lm(AveragebygenesymbolPeterHRLR[i,]~locoData[,6])
  Peter_GeneBehaviorCorrelations[i,1]<-summary.lm(temp)$coefficients[2,1]
  Peter_GeneBehaviorCorrelations[i,2]<-summary.lm(temp)$r.squared 
  Peter_GeneBehaviorCorrelations[i,3]<-summary.lm(temp)$coefficients[2,4]
  rm(temp)
  temp <- lm(AveragebygenesymbolPeterHRLR[i,]~locoData[,7])
  Peter_GeneBehaviorCorrelations[i,4]<-summary.lm(temp)$coefficients[2,1]
  Peter_GeneBehaviorCorrelations[i,5]<-summary.lm(temp)$r.squared 
  Peter_GeneBehaviorCorrelations[i,6]<-summary.lm(temp)$coefficients[2,4]
  rm(temp)
  temp <- lm(AveragebygenesymbolPeterHRLR[i,]~locoData[,10])
  Peter_GeneBehaviorCorrelations[i,7]<-summary.lm(temp)$coefficients[2,1]
  Peter_GeneBehaviorCorrelations[i,8]<-summary.lm(temp)$r.squared 
  Peter_GeneBehaviorCorrelations[i,9]<-summary.lm(temp)$coefficients[2,4]
  rm(temp)
}

write.csv(Peter_GeneBehaviorCorrelations, "Peter_GeneBehaviorCorrelations.csv")




# Cigdem's data:

#Creating color vector for plots
locoData$Lineage
#[1] HR LR HR LR LR HR HR LR HR HR LR LR IR IR IR IR IR IR

colorPhen<-Cigdem_sampleInfo$Rat_Type
colorPhen <-gsub("HR", "green3", colorPhen) #Making HRs green
colorPhen <-gsub("LR", "red", colorPhen) #Making LRs red


geneT <- t(AveragebygenesymbolCigdemHRLR)
str(geneT)
#pulling out column numbers associated with GOIs
goiCols <- which(Cigdem_Annotation$GeneSymbol%in%genesOfInterest==TRUE)

#Subsetting GOIs data
goiData <- cbind.data.frame(geneT[,goiCols])
colnames(goiData)<-Cigdem_Annotation$GeneSymbol[goiCols]
Cigdem_locoData<-cbind.data.frame(Cigdem_sampleInfo, goiData, stringsAsFactors=F )

dim(Cigdem_locoData)
colnames(Cigdem_locoData)
# [1] "Sample_ID"   "Core_ID"     "Sample_Desc" "Rat_Type"    "Treatment"   "SI_Score"    "Group"       "Sox9"        "Trhr"       
# [10] "Bmp4"        "Trh"         "C1qb"        "C1qc"        "C1qa"        "Cd4"         "Mfge8"       "Ucp2"        "Fxyd7"      
# [19] "Ncan"        "Cav1"    

#Create folder for boxplots

setwd("~/Documents/Microarray Gen/HRLR/HRLR_Studies/Cigdem_HRLR_RNAseq/GeneBehavior_Plots")

#This is me cheating a little for the automatic xlab names:

for(i in c(8:20)){
    tiff(paste0(colnames(Cigdem_locoData[i]), 
                "vs", "SocialInteraction", "Plot.tiff", sep=" "), res=300, compression="lzw",
         width = 6, height = 6, units = 'in')
    par(mai=c(.9,1.05,0.95,0.4))
    
    plot(Cigdem_locoData[,i]~Cigdem_locoData[,6], data=Cigdem_locoData, col="black", bg=colorPhen, pch=21, cex=1.5,
         xlab="Social Interaction Score", cex.lab=2, cex.main=2.2, ylab=paste(colnames(Cigdem_locoData[i]),": Log(2) FPM", sep=""))
    
    temp <- lm(Cigdem_locoData[,i]~Cigdem_locoData[,6])
    abline(temp, lwd=2)
    mtext(paste("R-squared=", signif(summary.lm(temp)$r.squared, 2), ", P-value=", signif(summary.lm(temp)$coefficients[2,4], 2)), cex=1.5, line=1)
    if(temp$coefficients[2] > 0){
      legend("bottomright", legend=c("HR","LR"), col=c("green3", "red"), 
             pch=19, cex=1.3)
    }else{legend("topright", legend=c("HR", "LR"), col=c("green3", "red"), 
                 pch=19, cex=1.3)}
    dev.off()
}

#There is a discrepancy between the SI score that Cigdem has and Fan's SI score for rat 53816 (HR_1)

Cigdem_locoData$SI_ScoreCorrected<-Cigdem_locoData$SI_Score
Cigdem_locoData$SI_ScoreCorrected[1]<-28.3333

setwd("~/Documents/Microarray Gen/HRLR/HRLR_Studies/Cigdem_HRLR_RNAseq/GeneBehavior_Plots/If_53816_is28")

for(i in c(8:20)){
  tiff(paste0(colnames(Cigdem_locoData[i]), 
              "vs", "SocialInteraction", "Plot.tiff", sep=" "), res=300, compression="lzw",
       width = 6, height = 6, units = 'in')
  par(mai=c(.9,1.05,0.95,0.4))
  
  plot(Cigdem_locoData[,i]~Cigdem_locoData$SI_ScoreCorrected, data=Cigdem_locoData, col="black", bg=colorPhen, pch=21, cex=1.5,
       xlab="Social Interaction Score", cex.lab=2, cex.main=2.2, ylab=paste(colnames(Cigdem_locoData[i]),": Log(2) FPM", sep=""))
  
  temp <- lm(Cigdem_locoData[,i]~Cigdem_locoData$SI_ScoreCorrected)
  abline(temp, lwd=2)
  mtext(paste("R-squared=", signif(summary.lm(temp)$r.squared, 2), ", P-value=", signif(summary.lm(temp)$coefficients[2,4], 2)), cex=1.5, line=1)
  if(temp$coefficients[2] > 0){
    legend("bottomright", legend=c("HR","LR"), col=c("green3", "red"), 
           pch=19, cex=1.3)
  }else{legend("topright", legend=c("HR", "LR"), col=c("green3", "red"), 
               pch=19, cex=1.3)}
  dev.off()
}


Cigdem_GeneBehaviorCorrelations<-matrix(0, nrow(AveragebygenesymbolCigdemHRLR), ncol=6)
row.names(Cigdem_GeneBehaviorCorrelations)<-row.names(AveragebygenesymbolCigdemHRLR)
colnames(Cigdem_GeneBehaviorCorrelations)<-c("SocialInteraction_Beta", "SocialInteraction_Rsquared", "SocialInteraction_Pval", "SocialInteractionCorrected_Beta", "SocialInteractionCorrected_Rsquared", "SocialInteractionCorrected_Pval")

for(i in c(1:nrow(AveragebygenesymbolCigdemHRLR))){
  temp <- lm(AveragebygenesymbolCigdemHRLR[i,]~Cigdem_locoData$SI_Score)
  Cigdem_GeneBehaviorCorrelations[i,1]<-summary.lm(temp)$coefficients[2,1]
  Cigdem_GeneBehaviorCorrelations[i,2]<-summary.lm(temp)$r.squared 
  Cigdem_GeneBehaviorCorrelations[i,3]<-summary.lm(temp)$coefficients[2,4]
  rm(temp)
  temp <- lm(AveragebygenesymbolCigdemHRLR[i,]~Cigdem_locoData$SI_ScoreCorrected)
  Cigdem_GeneBehaviorCorrelations[i,4]<-summary.lm(temp)$coefficients[2,1]
  Cigdem_GeneBehaviorCorrelations[i,5]<-summary.lm(temp)$r.squared 
  Cigdem_GeneBehaviorCorrelations[i,6]<-summary.lm(temp)$coefficients[2,4]
  rm(temp)
}

write.csv(Cigdem_GeneBehaviorCorrelations, "Cigdem_GeneBehaviorCorrelations.csv")


tiff(paste0("SI vs Phenotype.tiff", sep=" "), 
     res=300, compression="lzw",
     width = 4, height = 4.5, units = 'in')
par(mai=c(0.8,1,0.95,0.4))

boxplot(Cigdem_locoData$SI_Score~Cigdem_locoData$Rat_Type, data=Cigdem_locoData, 
        col=c("green3","red"), cex=1, pch=19, ylab="Social Interaction Score", cex.lab=1.3)

dev.off()


summary.lm(lm(Cigdem_locoData$SI_Score~Cigdem_locoData$Rat_Type))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                  19.334      2.575   7.508 6.87e-05 ***
#   Cigdem_locoData$Rat_TypeLR   -6.136      3.642  -1.685     0.13    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 5.758 on 8 degrees of freedom
# Multiple R-squared:  0.2619,	Adjusted R-squared:  0.1697 
# F-statistic: 2.839 on 1 and 8 DF,  p-value: 0.1305

summary.lm(lm(Cigdem_locoData$SI_ScoreCorrected~Cigdem_locoData$Rat_Type))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                  22.201      2.630   8.441 2.96e-05 ***
#   Cigdem_locoData$Rat_TypeLR   -9.003      3.719  -2.421   0.0418 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 5.881 on 8 degrees of freedom
# Multiple R-squared:  0.4228,	Adjusted R-squared:  0.3506 
# F-statistic: 5.859 on 1 and 8 DF,  p-value: 0.04182


tiff(paste0("SICorrected vs Phenotype.tiff", sep=" "), 
     res=300, compression="lzw",
     width = 4, height = 4.5, units = 'in')
par(mai=c(0.8,1,0.95,0.4))

boxplot(Cigdem_locoData$SI_ScoreCorrected~Cigdem_locoData$Rat_Type, data=Cigdem_locoData, 
        col=c("green3","red"), cex=1, pch=19, ylab="Social Interaction Score", cex.lab=1.3)

dev.off()


#LM of SI effects


out<-c(
  capture.output(
    for (i in 8:20){
      print(paste(colnames(Cigdem_locoData)[i], "vs. SI", sep="  "))
      print(summary.lm(lm(Cigdem_locoData[,i]~Cigdem_locoData$SI_ScoreCorrected, data=Cigdem_locoData)))
      print(" ")
    }
  )
)
cat(out, file="lm of SI_F43Data.txt", sep="\n", append=FALSE)
rm(out)