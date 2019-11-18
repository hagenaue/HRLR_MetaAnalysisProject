#Snooping around at how Zeisel's single cell data (2018) compares to other expression data for the hippocampus:
#Megan Hagenauer, 4-9-2018

setwd("~/Documents/Microarray Gen/Zeisel_WholeBrainscRNAseq")

ClusterExpression<-read.delim("Zeisel_MouseBrainAllCellTypes_ClusterExpression_ForR.txt", sep="\t", header=T, stringsAsFactors = F)

str(ClusterExpression)
# 'data.frame':	27998 obs. of  266 variables:
#   $ Gene   : chr  "Cbln2" "Ptchd2" "P2rx2" "Ptger4" ...
# $ ENT9   : num  13.976 0.391 7.609 0.929 8.74 ...
# $ ENT8   : num  0.4227 0 1.3814 0.0928 9 ...
# $ ENT6   : num  0.127 0 0.186 0.195 3.602 ...
# $ ENT5   : num  0.0849 0 0.5174 0.1274 1.7529 ...

colnames(ClusterExpression)[1]<-"GeneSymbol_Mouse"

ClusterMetaData<-read.delim("CellCluster_MetaData_ForR.txt", sep="\t", header=T, stringsAsFactors = F)
str(ClusterMetaData)
# 'data.frame':	265 obs. of  45 variables:
#   $ Bucket                   : chr  "/Users/sten/build_20171205/L4_Enteric_Neurons.loom" "/Users/sten/build_20171205/L4_Enteric_Neurons.loom" "/Users/sten/build_20171205/L4_Enteric_Neurons.loom" "/Users/sten/build_20171205/L4_Enteric_Neurons.loom" ...
# $ Class                    : chr  "Neurons" "Neurons" "Neurons" "Neurons" ...
# $ ClusterName              : chr  "ENT9" "ENT8" "ENT6" "ENT5" ...
# $ ClusterScore             : num  10 7.99 8.33 6 6.4 ...

setwd("~/Documents/Microarray Gen/HRLR/SecondaryAnalyses/HC_CoexpressionModules/ComparedToZeisel")

HippoSeq_DendrogramMarkerGenes<-read.csv("HippoSeq_DendrogramMarkerGenes2_forR.csv", header=T, stringsAsFactors = F)
str(HippoSeq_DendrogramMarkerGenes)

# 'data.frame':	248 obs. of  7 variables:
#   $ enriched            : chr  "dg_d-dg_v" "dg_d-dg_v" "dg_d-dg_v" "dg_d-dg_v" ...
# $ depleted            : chr  "ca4-ca3_d-ca3_v-ca2-ca1_d-ca1_v" "ca4-ca3_d-ca3_v-ca2-ca1_d-ca1_v" "ca4-ca3_d-ca3_v-ca2-ca1_d-ca1_v" "ca4-ca3_d-ca3_v-ca2-ca1_d-ca1_v" ...
# $ GeneSymbol_Mouse    : chr  "Ablim3" "Akap7" "Arhgap20" "Btg2" ...
# $ coronalIshAvailable : int  1 0 0 1 1 0 1 1 1 1 ...
# $ coronalIshConsistent: int  1 0 0 0 1 0 1 1 1 1 ...
# $ GeneSymbol_Human    : chr  "ABLIM3" "AKAP7" "ARHGAP20" "BTG2" ...
# $ EnrichedVSDepleted  : chr  "dg_d-dg_v__VS__ca4-ca3_d-ca3_v-ca2-ca1_d-ca1_v" "dg_d-dg_v__VS__ca4-ca3_d-ca3_v-ca2-ca1_d-ca1_v" "dg_d-dg_v__VS__ca4-ca3_d-ca3_v-ca2-ca1_d-ca1_v" "dg_d-dg_v__VS__ca4-ca3_d-ca3_v-ca2-ca1_d-ca1_v" ...



library(plyr)

HippoSeq_vs_ClusterExpression<-join(HippoSeq_DendrogramMarkerGenes, ClusterExpression, by="GeneSymbol_Mouse", type="inner")
str(HippoSeq_vs_ClusterExpression)

# 'data.frame':	232 obs. of  272 variables:
#   $ enriched            : chr  "dg_d-dg_v" "dg_d-dg_v" "dg_d-dg_v" "dg_d-dg_v" ...
# $ depleted            : chr  "ca4-ca3_d-ca3_v-ca2-ca1_d-ca1_v" "ca4-ca3_d-ca3_v-ca2-ca1_d-ca1_v" "ca4-ca3_d-ca3_v-ca2-ca1_d-ca1_v" "ca4-ca3_d-ca3_v-ca2-ca1_d-ca1_v" ...
# $ GeneSymbol_Mouse    : chr  "Ablim3" "Akap7" "Arhgap20" "Btg2" ...
# $ coronalIshAvailable : int  1 0 0 1 1 0 1 1 1 1 ...
# $ coronalIshConsistent: int  1 0 0 0 1 0 1 1 1 1 ...
# $ GeneSymbol_Human    : chr  "ABLIM3" "AKAP7" "ARHGAP20" "BTG2" ...
# $ EnrichedVSDepleted  : chr  "dg_d-dg_v__VS__ca4-ca3_d-ca3_v-ca2-ca1_d-ca1_v" "dg_d-dg_v__VS__ca4-ca3_d-ca3_v-ca2-ca1_d-ca1_v" "dg_d-dg_v__VS__ca4-ca3_d-ca3_v-ca2-ca1_d-ca1_v" "dg_d-dg_v__VS__ca4-ca3_d-ca3_v-ca2-ca1_d-ca1_v" ...
# $ ENT9                : num  0.615 0.183 0.343 1.462 0 ...

head(as.matrix(HippoSeq_vs_ClusterExpression[,c(8:272)]))

heatmap(as.matrix(HippoSeq_vs_ClusterExpression[,c(8:272)]))
#That seems pretty overwhelming...

HippoSeq_vs_ClusterExpression_Matrix<-as.matrix(HippoSeq_vs_ClusterExpression[,c(8:272)])
str(HippoSeq_vs_ClusterExpression_Matrix)

table(ClusterMetaData$Developmental_compartment)

Neural crest,  Spinal cord

table(ClusterMetaData$Tissue_CA1==0, ClusterMetaData$Developmental_compartment)
#Hmm... some of those are Neural Crest, or Mesoderm
#No Diencephalon, Mesencephalon, Rhombencephalon, or Spinal Cord
table(ClusterMetaData$Tissue_DentGyr==0, ClusterMetaData$Developmental_compartment)
#There are a few that are mesencephalon, rhombencephalon, neural crest...
#Still no spinal cord or diencephalon.

table(ClusterMetaData$Tissue_HC==0, ClusterMetaData$Developmental_compartment)
#Still no spinal cord or diencephalon. 

table(ClusterMetaData$Tissue_HC==0, ClusterMetaData$Tissue_DentGyr==0)
#Not completely overlapping
table(ClusterMetaData$Tissue_HC==0, ClusterMetaData$Tissue_CA1==0)
#Also not completely overlapping.

sum(ClusterMetaData$Region%in%c("Enteric nervous system", "Pons", "Medulla", "Spinal cord", "Cerebellum", "Dorsal root ganglion",  "Sympathetic ganglion", "Olfactory bulb"))
#[1] 115

LessInterestingClusters<-ClusterMetaData$Region%in%c("Enteric nervous system", "Pons", "Medulla", "Spinal cord", "Cerebellum", "Dorsal root ganglion",  "Sympathetic ganglion", "Olfactory bulb")

table(ClusterMetaData$Tissue_HC==0, LessInterestingClusters)
# #        FALSE TRUE
# FALSE    65    4
# TRUE     85  111
#There's still overlap! I wonder where...

cbind(ClusterMetaData$Tissue_HC==0, LessInterestingClusters)
# [129,] FALSE                    TRUE
# [133,] FALSE                    TRUE
# [134,] FALSE                    TRUE
# [233,] FALSE                    TRUE

ClusterMetaData[c(129, 133, 134, 233),]
#From olfactory bulb and cerebellum - so I guess maybe I shouldn't cut out those tissues.

table(ClusterMetaData$Tissue_DentGyr==0, LessInterestingClusters)
LessInterestingClusters
#       FALSE TRUE
# FALSE    45    2
# TRUE    105  113
#The cerebellum neurons again

table(ClusterMetaData$Tissue_CA1==0, LessInterestingClusters)
#Just the olfactory bulb astrocytes.

#Revised:
LessInterestingClusters<-ClusterMetaData$Region%in%c("Enteric nervous system", "Pons", "Medulla", "Spinal cord", "Dorsal root ganglion",  "Sympathetic ganglion")

HippoSeq_vs_ClusterExpression_Matrix_Relevant<-HippoSeq_vs_ClusterExpression_Matrix[, LessInterestingClusters==F]
heatmap(HippoSeq_vs_ClusterExpression_Matrix_Relevant)

HippoSeq_vs_ClusterMetaData_Relevant<-ClusterMetaData[LessInterestingClusters==F,]

write.csv(data.frame(HippoSeq_vs_ClusterExpression[, c(1:7)], HippoSeq_vs_ClusterExpression_Matrix_Relevant), "HippoSeq_vs_ClusterExpression_Matrix_Relevant.csv")
write.csv(HippoSeq_vs_ClusterMetaData_Relevant, "HippoSeq_vs_ClusterMetaData_Relevant.csv")

#Hard to visualize because some genes have very low expression, some have much higher. Let's transform it:

str(HippoSeq_vs_ClusterExpression_Matrix_Relevant)

HippoSeq_vs_ClusterExpression_Matrix_Relevant_Zscore<-t(scale(t(HippoSeq_vs_ClusterExpression_Matrix_Relevant)))

write.csv(data.frame(HippoSeq_vs_ClusterExpression[, c(1:7)], HippoSeq_vs_ClusterExpression_Matrix_Relevant_Zscore), "HippoSeq_vs_ClusterExpression_Matrix_Relevant_Zscore.csv")

#Well... that definitely wasn't as pretty as I had hoped, but still a little interesting.  Most of the region specific genes were neuronal in nature.


#Let's try another one:

MouseCoexpressionModules<-read.csv("MouseCoexpressionModules_forR.csv", header=T, stringsAsFactors = F) 
str(MouseCoexpressionModules)
# 'data.frame':	25697 obs. of  3 variables:
#   $ GeneSymbol_Human: chr  "39878" "1190002N15RIK" "1200013B08RIK" "1600021P15RIK" ...
# $ moduleHippo     : chr  "black" "black" "black" "black" ...
# $ moduleStriatum  : chr  "white" "cyan" "grey" "red" ...


ClusterExpression_toJoin<-data.frame(toupper(ClusterExpression$GeneSymbol_Mouse), ClusterExpression)
colnames(ClusterExpression_toJoin)[1]<-"GeneSymbol_Human"

MouseCoexpressionModules_vs_ClusterExpression<-join(MouseCoexpressionModules, ClusterExpression_toJoin, by="GeneSymbol_Human", type="inner")
str(MouseCoexpressionModules_vs_ClusterExpression)
image(as.matrix(MouseCoexpressionModules_vs_ClusterExpression[, c(5:269)]))
#that worked well, LOL.

MouseCoexpressionModules_vs_ClusterExpression_Matrix<-as.matrix(MouseCoexpressionModules_vs_ClusterExpression[, c(5:269)])
str(MouseCoexpressionModules_vs_ClusterExpression_Matrix)

MouseCoexpressionModules_vs_ClusterExpression_Matrix_Relevant<-MouseCoexpressionModules_vs_ClusterExpression_Matrix[, LessInterestingClusters==F]

write.csv(data.frame(MouseCoexpressionModules_vs_ClusterExpression[,c(1:4)], MouseCoexpressionModules_vs_ClusterExpression_Matrix_Relevant), "MouseCoexpressionModules_vs_ClusterExpression_Matrix_Relevant.csv")


MouseCoexpressionModules_vs_ClusterExpression_Matrix_Relevant_Zscore<-t(scale(t(MouseCoexpressionModules_vs_ClusterExpression_Matrix_Relevant)))

write.csv(data.frame(MouseCoexpressionModules_vs_ClusterExpression[,c(1:4)], MouseCoexpressionModules_vs_ClusterExpression_Matrix_Relevant_Zscore), "MouseCoexpressionModules_vs_ClusterExpression_Matrix_Relevant_Zscore.csv")
#Interesting.


Mean_MouseCoexpressionModules_vs_ClusterExpression<-matrix(0, length(table(MouseCoexpressionModules_vs_ClusterExpression$moduleHippo)), ncol(MouseCoexpressionModules_vs_ClusterExpression_Matrix_Relevant_Zscore))
row.names(Mean_MouseCoexpressionModules_vs_ClusterExpression)<-names(table(MouseCoexpressionModules_vs_ClusterExpression$moduleHippo))
colnames(Mean_MouseCoexpressionModules_vs_ClusterExpression)<-colnames(MouseCoexpressionModules_vs_ClusterExpression_Matrix_Relevant_Zscore)

for(i in c(1:ncol(MouseCoexpressionModules_vs_ClusterExpression_Matrix_Relevant_Zscore))){
Mean_MouseCoexpressionModules_vs_ClusterExpression[,i]<-tapply(MouseCoexpressionModules_vs_ClusterExpression_Matrix_Relevant_Zscore[,i],MouseCoexpressionModules_vs_ClusterExpression$moduleHippo, function(y) mean(y, na.rm=T))
}       

write.csv(Mean_MouseCoexpressionModules_vs_ClusterExpression, "Mean_MouseCoexpressionModules_vs_ClusterExpression.csv")
#Interesting - that was informative, but only for a small subset of modules.


HumanCoexpressionModules<-read.csv("HumanCoexpressionModules_forR.csv", header=T, stringsAsFactors = F) 
str(HumanCoexpressionModules)
# 'data.frame':	25697 obs. of  3 variables:
#   $ GeneSymbol_Human: chr  "39878" "1190002N15RIK" "1200013B08RIK" "1600021P15RIK" ...
# $ moduleHippo     : chr  "black" "black" "black" "black" ...
# $ moduleStriatum  : chr  "white" "cyan" "grey" "red" ...

HumanCoexpressionModules_vs_ClusterExpression<-join(HumanCoexpressionModules, ClusterExpression_toJoin, by="GeneSymbol_Human", type="inner")
str(HumanCoexpressionModules_vs_ClusterExpression)
image(as.matrix(HumanCoexpressionModules_vs_ClusterExpression[, c(4:268)]))
#that worked well, LOL.



HumanCoexpressionModules_vs_ClusterExpression_Matrix<-as.matrix(HumanCoexpressionModules_vs_ClusterExpression[, c(4:268)])
str(HumanCoexpressionModules_vs_ClusterExpression_Matrix)

HumanCoexpressionModules_vs_ClusterExpression_Matrix_Relevant<-HumanCoexpressionModules_vs_ClusterExpression_Matrix[, LessInterestingClusters==F]

write.csv(data.frame(HumanCoexpressionModules_vs_ClusterExpression[,c(1:3)], HumanCoexpressionModules_vs_ClusterExpression_Matrix_Relevant), "HumanCoexpressionModules_vs_ClusterExpression_Matrix_Relevant.csv")


HumanCoexpressionModules_vs_ClusterExpression_Matrix_Relevant_Zscore<-t(scale(t(HumanCoexpressionModules_vs_ClusterExpression_Matrix_Relevant)))

write.csv(data.frame(HumanCoexpressionModules_vs_ClusterExpression[,c(1:3)], HumanCoexpressionModules_vs_ClusterExpression_Matrix_Relevant_Zscore), "HumanCoexpressionModules_vs_ClusterExpression_Matrix_Relevant_Zscore.csv")
#Interesting.


Mean_HumanCoexpressionModules_vs_ClusterExpression<-matrix(0, length(table(HumanCoexpressionModules_vs_ClusterExpression$Module)), ncol(HumanCoexpressionModules_vs_ClusterExpression_Matrix_Relevant_Zscore))
row.names(Mean_HumanCoexpressionModules_vs_ClusterExpression)<-names(table(HumanCoexpressionModules_vs_ClusterExpression$Module))
colnames(Mean_HumanCoexpressionModules_vs_ClusterExpression)<-colnames(HumanCoexpressionModules_vs_ClusterExpression_Matrix_Relevant_Zscore)

for(i in c(1:ncol(HumanCoexpressionModules_vs_ClusterExpression_Matrix_Relevant_Zscore))){
  Mean_HumanCoexpressionModules_vs_ClusterExpression[,i]<-tapply(HumanCoexpressionModules_vs_ClusterExpression_Matrix_Relevant_Zscore[,i],HumanCoexpressionModules_vs_ClusterExpression$Module, function(y) mean(y, na.rm=T))
}       

write.csv(Mean_HumanCoexpressionModules_vs_ClusterExpression, "Mean_HumanCoexpressionModules_vs_ClusterExpression.csv")

####

#Alright, switching gears to take a peek at the top genes from the adult meta-analysis (FDR<0.10):

setwd("~/Documents/Microarray Gen/HRLR/SecondaryAnalyses/Zeisel_scRNASeq")

AdultMetaAnalysis_FDR10<-read.csv("Adult_FDR10_forR.csv", header=T, stringsAsFactors = F)
str(AdultMetaAnalysis_FDR10)
#'data.frame':	191 obs. of  9 variables:
# $ GeneSymbol: chr  "Tmem144" "Asb15" "Kif15" "Pkhd1l1" ...
# $ rawp      : num  3.04e-08 4.05e-07 4.27e-07 6.13e-07 7.87e-07 1.04e-06 2.03e-06 2.12e-06 2.86e-06 4.27e-06 ...
# $ BH        : num  0.000495 0.002317 0.002317 0.002491 0.002562 ...
# $ BY        : num  0.00509 0.02381 0.02381 0.0256 0.02632 ...
# $ estimate  : num  -3.57 2.82 -2.21 -3.17 2.78 ...
# $ SE        : num  0.644 0.556 0.437 0.636 0.563 ...
# $ CI_lb     : num  -4.83 1.73 -3.07 -4.42 1.68 ...
# $ CI_ub     : num  -2.3 3.91 -1.35 -1.92 3.89 ...
# $ datasets  : int  3 3 4 3 3 4 2 3 4 5 ...

colnames(AdultMetaAnalysis_FDR10)[1]<-"GeneSymbol_Mouse"

AdultMetaAnalysis_FDR10_vs_ClusterExpression<-join(AdultMetaAnalysis_FDR10, ClusterExpression, by="GeneSymbol_Mouse", type="left")
str(AdultMetaAnalysis_FDR10_vs_ClusterExpression)

AdultMetaAnalysis_FDR10_vs_ClusterExpression_Matrix<-as.matrix(AdultMetaAnalysis_FDR10_vs_ClusterExpression[,c(10:274)])

AdultMetaAnalysis_FDR10_vs_ClusterExpression_Matrix_Relevant<-AdultMetaAnalysis_FDR10_vs_ClusterExpression_Matrix[,LessInterestingClusters==F]

write.csv(data.frame(AdultMetaAnalysis_FDR10_vs_ClusterExpression[,c(1:9)], AdultMetaAnalysis_FDR10_vs_ClusterExpression_Matrix_Relevant), "AdultMetaAnalysis_FDR10_vs_ClusterExpression_Matrix_Relevant.csv")

AdultMetaAnalysis_FDR10_vs_ClusterExpression_Matrix_Relevant_Zscore<-t(scale(t(AdultMetaAnalysis_FDR10_vs_ClusterExpression_Matrix_Relevant)))

write.csv(data.frame(AdultMetaAnalysis_FDR10_vs_ClusterExpression[,c(1:9)], AdultMetaAnalysis_FDR10_vs_ClusterExpression_Matrix_Relevant_Zscore), "AdultMetaAnalysis_FDR10_vs_ClusterExpression_Matrix_Relevant_Zscore.csv")


write.csv(HippoSeq_vs_ClusterMetaData_Relevant, "ClusterMetaData_Relevant.csv")




