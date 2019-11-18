#Some script to explore the human and mouse coexpression networks that I downloaded (Park et al. BMC Systems Biology 2011, 5:43 & ...)
## Megan Hagenauer 6/16/2017

setwd("~/Documents/pdfs of articles")
MouseCoexpressionModules<-read.csv("Park_MouseHippocampusStriatumCoexpressionModules.csv", header=T, stringsAsFactors = F)
HumanCoexpressionModules<-read.csv("HumanHippocampus_GeneCoexpressionModules_S2.csv", header=T, stringsAsFactors = F)
HippoSeq_DendrogramMarkerGenes<-read.csv("HippoSeq_S1_DendrogramMarkerGenes.csv", header=T, stringsAsFactors=F)
HippoSeq_DGMossyCellMarkerGenes<-read.csv("Hipposeq_DGMossyCellMarkers.csv", header=T, stringsAsFactors=F)
HippoSeq_DGGranuleCellMarkerGenes<-read.csv("Hipposeq_DGGranuleCellMarkers.csv", header=T, stringsAsFactors=F)

setwd("~/Documents/Microarray Gen/Cell Type Paper")
CellTypeGenes<-read.csv("SupplTable1_CellTypeSpecificGenes_Master3.csv", header=T, stringsAsFactors = F)

library(plyr)

colnames(MouseCoexpressionModules)
head(MouseCoexpressionModules)
colnames(CellTypeGenes)
colnames(CellTypeGenes)[5]<-"GeneSymbol_Mouse"
colnames(CellTypeGenes)[4]<-"GeneSymbol_Human"
colnames(MouseCoexpressionModules)[3]<-"GeneSymbol_Human"

colnames(HumanCoexpressionModules)[3]<-"GeneSymbol_Human"
head(MouseCoexpressionModules)

colnames(HippoSeq_DendrogramMarkerGenes)[4]<-"GeneSymbol_Mouse"
head(HippoSeq_DendrogramMarkerGenes)

colnames(HippoSeq_DGMossyCellMarkerGenes)[3]<-"GeneSymbol_Mouse"
str(HippoSeq_DGMossyCellMarkerGenes)
table(HippoSeq_DGMossyCellMarkerGenes$enriched)

HippoSeq_DGMossyCellMarkerGenes2<-data.frame(HippoSeq_DGMossyCellMarkerGenes, toupper(HippoSeq_DGMossyCellMarkerGenes[,3]), stringsAsFactors = F)
colnames(HippoSeq_DGMossyCellMarkerGenes2)[8]<-"GeneSymbol_Human"

colnames(HippoSeq_DGGranuleCellMarkerGenes)[2]<-"GeneSymbol_Mouse"
str(HippoSeq_DGGranuleCellMarkerGenes)
HippoSeq_DGGranuleCellMarkerGenes2<-data.frame(HippoSeq_DGGranuleCellMarkerGenes, toupper(HippoSeq_DGGranuleCellMarkerGenes[,2]), stringsAsFactors = F)
str(HippoSeq_DGGranuleCellMarkerGenes2)
colnames(HippoSeq_DGGranuleCellMarkerGenes2)[9]<-"GeneSymbol_Human"

HippoSeq_DendrogramMarkerGenes2<-data.frame(HippoSeq_DendrogramMarkerGenes, toupper(HippoSeq_DendrogramMarkerGenes[,4]), paste(HippoSeq_DendrogramMarkerGenes2$enriched, "__VS__", HippoSeq_DendrogramMarkerGenes2$depleted, sep=""), stringsAsFactors = F)
str(HippoSeq_DendrogramMarkerGenes2)
colnames(HippoSeq_DendrogramMarkerGenes2)[8]<-"GeneSymbol_Human"
colnames(HippoSeq_DendrogramMarkerGenes2)[9]<-"EnrichedVSDepleted"
  
write.csv(HippoSeq_DendrogramMarkerGenes2, "HippoSeq_DendrogramMarkerGenes2.csv")
write.csv(HippoSeq_DGGranuleCellMarkerGenes2, "HippoSeq_DGGranuleCellMarkerGenes.csv")
write.csv(HippoSeq_DGMossyCellMarkerGenes2, "HippoSeq_DGMossyCellMarkerGenes.csv")
write.csv(MouseCoexpressionModules, "MouseCoexpressionModules.csv")
write.csv(HumanCoexpressionModules, "HumanCoexpressionModules.csv")


CellTypeGenesVsMouseCoexpressionModules<-join(CellTypeGenes,MouseCoexpressionModules, by="GeneSymbol_Human", type="inner")
dim(CellTypeGenesVsMouseCoexpressionModules)
colnames(CellTypeGenesVsMouseCoexpressionModules)
table(CellTypeGenesVsMouseCoexpressionModules$CellType_Primary, CellTypeGenesVsMouseCoexpressionModules$moduleHippo)

CellTypeGenesVsMouseCoexpressionModules<-join(CellTypeGenes,MouseCoexpressionModules, by="GeneSymbol_Human", type="full")
table(CellTypeGenesVsMouseCoexpressionModules$CellType_Primary, CellTypeGenesVsMouseCoexpressionModules$moduleHippo)

HippoSeq_DendrogramMarkerGenesVsMouseCoexpressionModules<-join(HippoSeq_DendrogramMarkerGenes2, MouseCoexpressionModules, by="GeneSymbol_Human", type="full")
table(HippoSeq_DendrogramMarkerGenesVsMouseCoexpressionModules$EnrichedVSDepleted, HippoSeq_DendrogramMarkerGenesVsMouseCoexpressionModules$moduleHippo)

HippoSeq_DGMossyCellMarkerGenesVsMouseCoexpressionModules<-join(HippoSeq_DGMossyCellMarkerGenes2, MouseCoexpressionModules, by="GeneSymbol_Human", type="full")

HippoSeq_DGGranuleCellMarkerGenesVsMouseCoexpressionModules<-join(HippoSeq_DGGranuleCellMarkerGenes2, MouseCoexpressionModules, by="GeneSymbol_Human", type="full")
colnames(HippoSeq_DGGranuleCellMarkerGenesVsMouseCoexpressionModules)

setwd("~/Documents/Microarray Gen/HRLR/HC_CoexpressionModules")


write.csv(rbind(table(CellTypeGenesVsMouseCoexpressionModules$CellType_Primary, CellTypeGenesVsMouseCoexpressionModules$moduleHippo), table(CellTypeGenesVsMouseCoexpressionModules$moduleHippo)), "SummaryTable_CellTypeGenesVsHCMouseCoexpressionModule.csv")


write.csv(rbind(table(CellTypeGenesVsMouseCoexpressionModules$Tag, CellTypeGenesVsMouseCoexpressionModules$moduleHippo), table(CellTypeGenesVsMouseCoexpressionModules$moduleHippo)), "SummaryTable_CellTypeTagVsHCMouseCoexpressionModule.csv")

write.csv(rbind(table(HippoSeq_DendrogramMarkerGenesVsMouseCoexpressionModules$EnrichedVSDepleted, HippoSeq_DendrogramMarkerGenesVsMouseCoexpressionModules$moduleHippo), table(HippoSeq_DendrogramMarkerGenesVsMouseCoexpressionModules$moduleHippo)), "SummaryTable_HippoSeqEnrichedVsHCMouseCoexpressionModule.csv")

write.csv(rbind(table(HippoSeq_DGMossyCellMarkerGenesVsMouseCoexpressionModules$enriched, HippoSeq_DGMossyCellMarkerGenesVsMouseCoexpressionModules$moduleHippo), table(HippoSeq_DGMossyCellMarkerGenesVsMouseCoexpressionModules$moduleHippo)), "SummaryTable_HippoSeq_DGMossyCellMarkerGenesVsMouseCoexpressionModules.csv")

write.csv(rbind(table(is.na(HippoSeq_DGGranuleCellMarkerGenesVsMouseCoexpressionModules$gene_id)==F, HippoSeq_DGGranuleCellMarkerGenesVsMouseCoexpressionModules$moduleHippo), table(HippoSeq_DGGranuleCellMarkerGenesVsMouseCoexpressionModules$moduleHippo)), "SummaryTable_HippoSeq_DGGranuleCellMarkerGenesVsMouseCoexpressionModules.csv")


CellTypeGenesVsHumanCoexpressionModules<-join(CellTypeGenes,HumanCoexpressionModules, by="GeneSymbol_Human", type="full")
dim(CellTypeGenesVsHumanCoexpressionModules)
table(CellTypeGenesVsHumanCoexpressionModules$CellType_Primary, CellTypeGenesVsHumanCoexpressionModules$Module)

HippoSeq_DendrogramMarkerGenesVsHumanCoexpressionModules<-join(HippoSeq_DendrogramMarkerGenes2, HumanCoexpressionModules, by="GeneSymbol_Human", type="full")

HippoSeq_DGGranuleCellMarkerGenesVsHumanCoexpressionModules<-join(HippoSeq_DGGranuleCellMarkerGenes2, HumanCoexpressionModules, by="GeneSymbol_Human", type="full")

HippoSeq_DGMossyCellMarkerGenesVsHumanCoexpressionModules<-join(HippoSeq_DGMossyCellMarkerGenes2, HumanCoexpressionModules, by="GeneSymbol_Human", type="full")

write.csv(rbind(table(CellTypeGenesVsHumanCoexpressionModules$CellType_Primary, CellTypeGenesVsHumanCoexpressionModules$Module), table(CellTypeGenesVsHumanCoexpressionModules$Module)), "SummaryTable_CellTypeGenesVsHCHumanCoexpressionModule.csv")


write.csv(rbind(table(CellTypeGenesVsHumanCoexpressionModules$Tag, CellTypeGenesVsHumanCoexpressionModules$Module), table(CellTypeGenesVsHumanCoexpressionModules$Module)), "SummaryTable_CellTypeTagVsHCHumanCoexpressionModule.csv")


write.csv(rbind(table(HippoSeq_DendrogramMarkerGenesVsHumanCoexpressionModules$EnrichedVSDepleted, HippoSeq_DendrogramMarkerGenesVsHumanCoexpressionModules$Module), table(HippoSeq_DendrogramMarkerGenesVsHumanCoexpressionModules$Module)), "SummaryTable_HippoSeqEnrichedVsHCHumanCoexpressionModule.csv")

write.csv(rbind(table(is.na(HippoSeq_DGGranuleCellMarkerGenesVsHumanCoexpressionModules$gene_id)==F, HippoSeq_DGGranuleCellMarkerGenesVsHumanCoexpressionModules$Module), table(HippoSeq_DGGranuleCellMarkerGenesVsHumanCoexpressionModules$Module)), "SummaryTable_HippoSeqDGGranuleCellVsHCHumanCoexpressionModule.csv")

write.csv(rbind(table(HippoSeq_DGMossyCellMarkerGenesVsHumanCoexpressionModules$enriched, HippoSeq_DGMossyCellMarkerGenesVsHumanCoexpressionModules$Module), table(HippoSeq_DGMossyCellMarkerGenesVsHumanCoexpressionModules$Module)), "SummaryTable_HippoSeqDGMossyCellVsHCHumanCoexpressionModule.csv")


###
#Just out of curiousity:

Adult_GSEAresults<-read.delim("Adult_GSEA_100000perm.txt", sep="\t", header=T, stringsAsFactors = F)
str(Adult_GSEAresults)
#hmmm... I guess that will take a little more effort - it is reading in the full list of gene symbols as a character string. 

Adult_GSEAresults_PositiveRegCellProlif<-read.csv("Adult_GSEA_PositiveRegCellProlif.csv", stringsAsFactors = F, header=F)

Adult_GSEAresults_PositiveRegCellProlif[,1]%in%HumanCoexpressionModules$GeneSymbol_Human
colnames(Adult_GSEAresults_PositiveRegCellProlif)[1]<-"GeneSymbol_Human"
HumanCoexpressionModulesVsAdult_GSEAresults_PositiveRegCellProlif<-join(HumanCoexpressionModules, Adult_GSEAresults_PositiveRegCellProlif, by="GeneSymbol_Human", type="inner")
table(HumanCoexpressionModulesVsAdult_GSEAresults_PositiveRegCellProlif$Module)
# M1  M16  M18  M19  M20  M21  M24   M8   M9 none 
# 1    1    2    2    1    2    2    6    1   24 

table(HumanCoexpressionModules$Module)
# M1  M10  M11  M12  M13  M14  M15  M16  M17  M18  M19   M2  M20  M21  M22  M23  M24   M3   M4   M5   M6   M7   M8   M9 none 
# 1225   48  150   72  148  128   40  150   42   40  576   38  459  100   74  147  253  160   35  158   30   49  732  180 6547 

MouseCoexpressionModulesVsAdult_GSEAresults_PositiveRegCellProlif<-join(MouseCoexpressionModules, Adult_GSEAresults_PositiveRegCellProlif, by="GeneSymbol_Human", type="inner")
table(MouseCoexpressionModulesVsAdult_GSEAresults_PositiveRegCellProlif$moduleHippo)

# black          blue         brown      darkgrey darkturquoise         green   greenyellow          grey     lightcyan   lightyellow 
# 4             7            12             1             2             5             1            22             1             1 
# magenta     turquoise        yellow 
# 1            16             2 

table(MouseCoexpressionModules$moduleHippo)

# black           blue          brown           cyan      darkgreen       darkgrey    darkmagenta darkolivegreen        darkred 
# 520           2592           2820            229            131            492             43             44            154 
# darkturquoise   green    greenyellow           grey         grey60      lightcyan    lightyellow        magenta         orange 
# 128           1249            414           8445            184            696            604            469            100 
# paleturquoise royalblue    saddlebrown        sienna3        skyblue       skyblue3      steelblue      turquoise         violet 
# 65            163             73             39             73             34             67           4574             60 
# white         yellow    yellowgreen 
# 76           1120             39 


Adult_GSEAresults_PositiveRegCellProlif[,1]%in%HippoSeq_DGGranuleCellMarkerGenes2$GeneSymbol_Human
#None.
Adult_GSEAresults_PositiveRegCellProlif[,1]%in%HippoSeq_DGMossyCellMarkerGenes2$GeneSymbol_Human

Adult_GSEAresults_PositiveRegCellProlif[Adult_GSEAresults_PositiveRegCellProlif[,1]%in%HippoSeq_DendrogramMarkerGenes2$GeneSymbol_Human,1]
#Only 1 - [1] "S100B"

CellTypeGenes[CellTypeGenes$GeneSymbol_Human%in%Adult_GSEAresults_PositiveRegCellProlif[,1],]
#All different cell types - neurons, astrocytes, mural, endothelial, microglia.


 

CellTypeGenes[sample(1:3000, 10),]
