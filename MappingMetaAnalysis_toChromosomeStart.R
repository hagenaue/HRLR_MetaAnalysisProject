#Comparing the exploratory locomotion QTL results from the Zhou et al. 2019 (PNAS) study using our HRxLR F2 intercross to the HR/LR hippocampal transcriptomic meta-analysis results
#Megan Hagenauer, 08/2019
#R-studio version 1.0.153

###########################

#First: Accessing the chromosomal coordinates for all genes so that we can determine their distance to the QTLs 
library(org.Rn.eg.db)
library(plyr)

#Since the meta-analysis results are annotated by symbol, I need to map symbol to Entrez ID (which is what the rest of the annotation package uses)
x <- org.Rn.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

GeneSymbol<-unlist(xx, use.names=FALSE)
EntrezGeneID<-rep(names(xx), lengths(xx))
table(lengths(xx))

# 1 
# 42306
#A 1:1 mapping. Nice...

EntrezVsGeneSymbol<-data.frame(EntrezGeneID, GeneSymbol, stringsAsFactors=F)

#Note: the function CHRLOC retrieves the coordinates for the start and CHRLOCEND retrieves the end for each gene
x<-org.Rn.egCHRLOC
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

#Messing around to see how it works:
xx[1]
# $`24152`
# 3 
# 150492009
xx$"24152"
# 3 
# 150492009 

xx$"24152"[1]
# 3 
# 150492009 
names(xx$"24152"[1])
#[1] "3"

xx$"24152"[[1]]
#[1] 150492009

head(xx)
# $`24152`
# 3 
# 150492009 
# 
# $`24153`
# 4         4 
# 154309425 154423167 

ChromosomalLoc<-unlist(xx, use.names=FALSE)
head(ChromosomalLoc)
# [1]  150492009  154309425  154423167  128027879 -260124417  -88392247

hist(ChromosomalLoc)
#Why are half of the values negative?
#From the manual:
#"Chromosomal locations on both the sense and antisense strands are measured as the number of base pairs from the p (5’ end of the sense strand) to q (3’ end of the sense strand) arms. Chromosomal locations on the antisense strand have a leading "-" sign (e. g. -1234567).
#Since some genes have multiple start sites, this field can map to multiple locations."
str(ChromosomalLoc)
#int [1:18391] 150492009 154309425 154423167 128027879 -260124417 -88392247 -99120378 -49836685 -49836064 -49836801 ...
#How do I grab the chromosome as well as the position?
ChromosomalLoc<-unlist(xx, use.names=TRUE)
head(ChromosomalLoc)

#Interesting - includes the Entrez symbol as the name with some random decimal after it... but no chromosome.
24152.3    24153.4    24153.4    24157.8    24158.2   24159.10 
150492009  154309425  154423167  128027879 -260124417  -88392247 
#Oh - I think the random decimal is the chromosome. Weird. Maybe I can just use the decimal is delineation.
#Oh crap - if I just use the decimal as delineation, the .10 ends up being 1 instead of 10. Gah. That sucks. I didn't catch it the the first time around, and now I have to re-do everything. Grrr....
#And there isn't any way to use the order to clarify which .1's are actually 1 or 10. 

#What if I try to treat it as a character vector so that no assumptions are made?
strsplit((names(ChromosomalLoc)[6]), split="[.]")
# [[1]]
# [1] "24159" "10"

#O.k., that works.
length(ChromosomalLoc)
#[1] 18391
matrix(unlist(strsplit((names(ChromosomalLoc)[6:10]), split="[.]")), nrow=5, ncol=2, byrow=TRUE)
# [,1]    [,2]
# [1,] "24159" "10"
# [2,] "24161" "4" 
# [3,] "24161" "6" 
# [4,] "24161" "6" 
# [5,] "24161" "6" 
#That works

ChromosomalLoc_2<-matrix(unlist(strsplit((names(ChromosomalLoc)), split="[.]")), nrow=18391, ncol=2, byrow=TRUE)

EntrezID_Chromosome_ChromosomalLoc<-data.frame("EntrezID_Chromosome"=ChromosomalLoc_2[,1], "Chromosome"=ChromosomalLoc_2[,2],"ChromosomeLocation"=ChromosomalLoc)
head(EntrezID_Chromosome_ChromosomalLoc)
# EntrezID_Chromosome Chromosome ChromosomeLocation
# 1               24152          3          150492009
# 2               24153          4          154309425
# 3               24153          4          154423167
# 4               24157          8          128027879

setwd("~/Documents/Microarray Gen/HRLR/SecondaryAnalyses/RatGenomeDatabase_eQTL")
#I just want to take a look at this directly:
write.csv(EntrezID_Chromosome_ChromosomalLoc, "ChromosomalLoc_Rnor_6.csv")
#Looks good.

#################################################

EntrezGeneID<-rep(names(xx), lengths(xx))
head (EntrezGeneID)
#[1] "24152" "24153" "24153" "24157" "24158" "24159"

table(lengths(xx))
# 1     2     3     4     5     6     7     8    11 
# 16444   756    94    14     9     3     1     2     1
#Hmmm....Yep, Definitely not a 1 to 1 correspondence. Booger.

#What if I just grabbed the first listed?  They should all be within close proximity to each other, and I am interested in +/- 1 MB from the gene.

EntrezVsGeneSymbol[1,2]
[1] "Slc39a4l"
EntrezVsGeneSymbol[1,1]
[1] "100008565"

xx$"100008565"
#Huh. Crap.  That doesn't seem right - not all Entrez IDs will have gene symbols, but all genes should have chromosomal locations...
#Am I sure that xx is indexed by Entrez ID?

#Sanity Check:
#NCBI database: entrez gene Id 100008565 is Slc39a41, not in current annotation release
#So I guess that answers it. Sort of. Weird.

EntrezVsGeneSymbol[2000,1]

tail(EntrezVsGeneSymbol)
#Sanity Check: This matches the NCBI database:
#47850        94342       Bag6
xx$"47850"
#NULL

head(names(xx))
[1] "24152" "24153" "24157" "24158" "24159" "24161"
head(as.numeric(names(xx)))
#[1] 24152 24153 24157 24158 24159 24161
hist(as.numeric(names(xx)))
min(as.numeric(names(xx)))
#[1] 24152
max(as.numeric(names(xx)))
#[1] 108348453

#Sanity Check: 
#NCBI database Gene entrez id 24152 is Asip on Chr 3, Rnor_6: 150492010..150579870, Rnor_5: 56860395..156949277
#The start that we are getting in the output is: 150492009 (so the annotation is Rnor_6)


#################################
#The exploratory locomotion QTL results from the Zhou et al. 2019 study use an old genome build (Rnor5). 

#I ended up using the NCBI Genome Remapping Service to provide Rnor_5 analagous start information for each gene so that I could properly compare it to the output from the Zhou et al. study. I'm sure that this mapping should be provided somewhere else in a more official capacity, but I couldn't find it (?!).
#Link to the converter:
#https://www.ncbi.nlm.nih.gov/genome/tools/remap#tab=asm&src_org=Rattus%20norvegicus&src_asm=GCF_000001895.5&tgt_asm=GCF_000001895.4&min_ratio=0.5&max_ratio=2.0&allow_locations=true&merge_fragments=true&in_fmt=guess&out_fmt=guess&genome_workbench=true
#Accessed 8/8/2019
#Input was a bed file containing the Rnor_6 chromosome and start information for each Entrez ID. For simplicity sake, the end was set as the start + 100 bp.
#Default parameters were chosen: 
#Minimum ratio of bases that must be remapped: 0.5
#Maximum ratio for difference between source length and target length: 2.0
#Allow multiple locations to be returned: Y
#Merge Fragments: Y

#To make it easier to combine this annotation with the QTL results, I removed several columns that weren't necessary in Excel and also added a few columns: 
#A rounded version of the annotation for the Rnor5_start (because the QTL results use 1 MB bins)


NCBIconverter_GeneLocations_Rnor6_toRnor5<-read.csv("NCBIconverter_GeneLocations_Rnor6_toRnor5_forR.csv", header=T, stringsAsFactors = F)
head(NCBIconverter_GeneLocations_Rnor6_toRnor5)
feat_name Rnor6_Chr Rnor6_Start Rnor6_Stop Rnor5_Chr Rnor5_Start Rnor5_Stop Rnor5_Start_Round Feat_name_Rnor5_Start
# 1 3_150492010         3   150492010  150492109         3   156860395  156860494         157000000           3_157000000
# 2 4_154309426         4   154309426  154309525         4   221393233  221393332         221000000           4_221000000
# 3 4_154423168         4   154423168  154423267         4   221506287  221506386         222000000           4_222000000
# 4 8_128027880         8   128027880  128027979         8   127234641  127234740         127000000           8_127000000
# 5 2_260124418         2   260124418  260124517         2   278788485  278788584         279000000           2_279000000
# 6 10_88392248        10    88392248   88392347        10    88185495   88185594          88000000           10_88000000

#Then I had to recombine this with the EntrezIDs from the Rnor6 annotation, because the NCBI converter would only accept bed files (tab-delimited feat_name and coordinates - no extraneous information):
#Note that this is a version of the "ChromosomalLoc_Rnor_6.csv" output from earlier, but with a column called "feat_name" which is a combination of the chromosome and starting coordinates that were used in the bed file submitted to the NCBI converter.
GeneLocations_Rnor6<-read.csv("GeneLocations_Rnor6_forR.csv", header=T, stringsAsFactors = F)
head(GeneLocations_Rnor6)
# EntrezID_Chromosome Chromosome ChromosomeLocation   feat_name
# 1               24152          3          150492009 3_150492009
# 2               24153          4          154309425 4_154309425
# 3               24153          4          154423167 4_154423167
# 4               24157          8          128027879 8_128027879
# 5               24158          2          260124417 2_260124417
# 6               24159         10           88392247 10_88392247

#Oh weird - the feature names and start locations are always -1 bp shifted in the NCBI output. Huh. I'll have to correct for that if I'm going to join the two outputs. I'm feeling lazy, so I just did that quickly in Excel and read it in again:

GeneLocations_Rnor6<-read.csv("GeneLocations_Rnor6_forR.csv", header=T, stringsAsFactors = F)
head(GeneLocations_Rnor6)
# 1               24152          3          150492009 3_150492010
# 2               24153          4          154309425 4_154309426
# 3               24153          4          154423167 4_154423168

GeneLocations_Rnor6_Rnor5<-join(GeneLocations_Rnor6, NCBIconverter_GeneLocations_Rnor6_toRnor5, by="feat_name", type="left")
head(GeneLocations_Rnor6_Rnor5)

# EntrezID_Chromosome Chromosome ChromosomeLocation   feat_name Rnor6_Chr Rnor6_Start Rnor6_Stop Rnor5_Chr Rnor5_Start Rnor5_Stop
# 1               24152          3          150492009 3_150492010         3   150492010  150492109         3   156860395  156860494
# 2               24153          4          154309425 4_154309426         4   154309426  154309525         4   221393233  221393332
# 3               24153          4          154423167 4_154423168         4   154423168  154423267         4   221506287  221506386
# 4               24157          8          128027879 8_128027880         8   128027880  128027979         8   127234641  127234740
# 5               24158          2          260124417 2_260124418         2   260124418  260124517         2   278788485  278788584
# 6               24159         10           88392247 10_88392248        10    88392248   88392347        10    88185495   88185594
# Rnor5_Start_Round Feat_name_Rnor5_Start
# 1         157000000           3_157000000
# 2         221000000           4_221000000
# 3         222000000           4_222000000
# 4         127000000           8_127000000
# 5         279000000           2_279000000
# 6          88000000           10_88000000

#Cool, looks good. 

colnames(GeneLocations_Rnor6_Rnor5)[1]<-"EntrezGeneID"

#Final two joins:
EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5<-join(EntrezVsGeneSymbol, GeneLocations_Rnor6_Rnor5, by="EntrezGeneID", type="left")
head(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5)

# EntrezGeneID GeneSymbol Chromosome ChromosomeLocation   feat_name Rnor6_Chr Rnor6_Start Rnor6_Stop Rnor5_Chr Rnor5_Start Rnor5_Stop
# 1    100008565   Slc39a4l       <NA>                 NA        <NA>      <NA>          NA         NA      <NA>          NA         NA
# 2    100034253      Gnl3l          X           20037557  X_20037558         X    20037558   20037657         X    20785115   20785214
# 3    100036582    Olr1867          8           20565731  8_20565732         8    20565732   20565831         8    20639977   20640076
# 4    100036765     Ccdc92         12           37211315 12_37211316        12    37211316   37211415        12    39085314   39085413
# 5    100049583      Trex1          8          117796126 8_117796127         8   117796127  117796226         8   117147872  117147971
# 6    100101342 Vom2r-ps11       <NA>                 NA        <NA>      <NA>          NA         NA      <NA>          NA         NA
# Rnor5_Start_Round Feat_name_Rnor5_Start
# 1                NA                  <NA>
#   2          21000000            X_21000000
# 3          21000000            8_21000000
# 4          39000000           12_39000000
# 5         117000000           8_117000000
# 6                NA                  <NA>

str(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5)
#Looks good

write.csv(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5, "EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5.csv")
#Alright, I finally have a map of the gene locations in Rnor5 and Rnor6 coordinates. That took waaaay more effort than I had originally expected.

########################################################

#Alright, combining it with our genetic results:

#This analysis used the HRLR F2 R-qtl H-K output from the Zhou lab sent to us prior to submitting the manuscript (20180720). It provided LODs for each 1 MB bin across the entire genome as well as LODs for additional SNPs. However, there are 11 spots in the genome where there isn't an LOD for 2 MB... but there is often a SNP in between.  To simplify things for this analysis, I chose to treat the SNPS that are in between as the missing 1 MB bin. If there were multiple SNPs in between, I chose the SNP with the largest LOD.
#Note: The only situation where this actually matters for our analysis is the bin surrounding Mfge8.
# 
# c1.loc140
# S141117448
# c1.loc142
# 
# c1.loc142
# S142512409
# S143131279
# c1.loc144
# 
# c1.loc259
# S260116369
# S260340076
# c1.loc261
# 
# c3.loc66
# S66926470
# c3.loc68
# 
# c3.loc157
# S157390309
# S157929831
# c3.loc159
# 
# c4.loc176
# S177350346
# S177375498
# c4.loc178
# 
# c13.loc47
# S47850178
# c13.loc49
# 
# c14.loc36
# S36042863
# S36943390
# c14.loc38
# 
# c14.loc45
# S45942095
# c14.loc47
# 
# c15.loc98
# S98990779
# c15.loc100
# 
# c17.loc40
# S40712180
# S40717025
# S40900079
# S40959776
# c17.loc42

#Note: the file that I'm reading in is also adapted as follows:
#1. I removed the columns for other regression analyses (e.g., H-K controlling for kinship)
#2. I have added a column that calculated the maximum LOD +/-1 MB from each coordinate (1 MB bin)
#3. I have added a feature location column using the same notation as above (chromosome_coordinate)
HRLR_F2_R_qtl_Results<-read.csv("20180720_HRLR_F2_R-qtl_Results_forR_v2.csv", header=T, stringsAsFactors=F)
head(HRLR_F2_R_qtl_Results)
# Loci LOD..H.K. MaxLOD_1MBUpDown Chromosome Loc    Start feat_location
# 1  c1.loc5 0.7879744         0.945018          1   5  5000000     1_5000000
# 2  c1.loc6 0.9450180         1.094077          1   6  6000000     1_6000000
# 3  c1.loc7 1.0940767         1.214148          1   7  7000000     1_7000000
# 4  c1.loc8 1.2141475         1.261184          1   8  8000000     1_8000000
# 5  c1.loc9 1.2611842         1.300115          1   9  9000000     1_9000000
# 6 c1.loc10 1.3001146         1.331538          1  10 10000000    1_10000000

colnames(HRLR_F2_R_qtl_Results)[7]<-"Feat_name_Rnor5_Start"

HRLR_F2_R_qtl_Results_toJoin<-HRLR_F2_R_qtl_Results[, -(c(4:6))]

head(HRLR_F2_R_qtl_Results_toJoin)
str(HRLR_F2_R_qtl_Results_toJoin)

str(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5)

EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_wQTLresults<-join(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5, HRLR_F2_R_qtl_Results_toJoin, by="Feat_name_Rnor5_Start", type="left")

head(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_wQTLresults)

write.csv(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_wQTLresults, "EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_wQTLresults.csv")

########################################################

#...and comparing our results to the HR/LR meta-analysis output:

setwd("~/Documents/Microarray Gen/HRLR/ThesisMetaAnalysisOutput")

#Note: I used the version of the results that has the genename/dates confusion problem fixed... but if you open any of the files in Excel the problem simply recreates itself.
AdultMeta_CohDandPval<-read.csv("AdultMeta_CohDandPval_GeneDatesFixed.csv", header=T, stringsAsFactors = F)
head(AdultMeta_CohDandPval)

AdultMeta_CohDandPval_wQTLresults<-join(AdultMeta_CohDandPval, EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_wQTLresults, by="GeneSymbol", type="left")

head(AdultMeta_CohDandPval_wQTLresults)

write.csv(AdultMeta_CohDandPval_wQTLresults, "AdultMeta_CohDandPval_wQTLresults.csv")

#Just to peek... but clearly turned out to be pretty noisy to look at.
plot(log10(AdultMeta_CohDandPval_wQTLresults$rawp)~AdultMeta_CohDandPval_wQTLresults$MaxLOD_1MBUpDown)

#So let's do some basic overlap analysis:

unique(AdultMeta_CohDandPval_wQTLresults$GeneSymbol[which(AdultMeta_CohDandPval_wQTLresults$BH<0.05 & AdultMeta_CohDandPval_wQTLresults$MaxLOD_1MBUpDown>4)])
#[1] "C2cd3"   "Ezr"     "Hmgn5b"  "Pkhd1l1" "Trhr"    "Ucp2"  
# 6 genes

unique(AdultMeta_CohDandPval_wQTLresults$GeneSymbol[which(AdultMeta_CohDandPval_wQTLresults$BH<0.05 & AdultMeta_CohDandPval_wQTLresults$MaxLOD_1MBUpDown<4)])
# [1] "Aar2"      "Acox3"     "Acss2"     "Afmid"     "Blvra"     "C1qa"      "Car9"      "Cav1"      "Ccdc137"   "Ddx20"     "Dnaaf3"
# [12] "Egfem1"    "Epb41l4a"  "Etv4"      "Fmo5"      "Fn3k"      "Fxyd7"     "LOC363337" "Mfge8"     "Mki67"     "Mvb12b"    "Myh6" # [23] "Ntn4"      "Nudt4"     "Oard1"     "Pkib"      "Pld4"      "Rhpn2"     "Rnls"      "Robo3"     "Rpl17"     "Samd5"     "Slc19a3"  
# [34] "Slc27a1"   "Slc39a12"  "Slc4a11"   "Slc9a3r1"  "Tdg"       "Tek"       "Tmem144"   "Tmem2"     "Tnnt1"     "Trappc6a"  "Trmt10a"  
# [45] "Tuba8"     "Tubg1"     "Uhrf1"     "Zfp110"    "Zfp90"     "Selenop" 
# 50 genes

length(unique(AdultMeta_CohDandPval_wQTLresults$GeneSymbol[which(AdultMeta_CohDandPval_wQTLresults$BH>0.05 & AdultMeta_CohDandPval_wQTLresults$MaxLOD_1MBUpDown>4)]))
#[1] 303

length(unique(AdultMeta_CohDandPval_wQTLresults$GeneSymbol[which(AdultMeta_CohDandPval_wQTLresults$BH>0.05 & AdultMeta_CohDandPval_wQTLresults$MaxLOD_1MBUpDown<4)]))
#[1]  11291

11291+303+50+6
#[1] 11650

#Just double-checking that things look correct:
head(unique(AdultMeta_CohDandPval_wQTLresults$GeneSymbol[which(AdultMeta_CohDandPval_wQTLresults$BH>0.05 & AdultMeta_CohDandPval_wQTLresults$MaxLOD_1MBUpDown<4)]))

head(unique(AdultMeta_CohDandPval_wQTLresults[which(AdultMeta_CohDandPval_wQTLresults$BH>0.05 & AdultMeta_CohDandPval_wQTLresults$MaxLOD_1MBUpDown<4),]))

head(AdultMeta_CohDandPval_wQTLresults[which(AdultMeta_CohDandPval_wQTLresults$BH>0.05 & AdultMeta_CohDandPval_wQTLresults$MaxLOD_1MBUpDown<4),])

AdultMeta_CohDandPval_wQTLresults[which(AdultMeta_CohDandPval_wQTLresults$BH>0.05 & AdultMeta_CohDandPval_wQTLresults$MaxLOD_1MBUpDown<4),][c(1:100),]

#I double-checked in Excel too - yes, 11,650 does look like the number of non-duplicated gene symbols with accompanying adult meta-analysis results which mapped to Rnor5 chromosomal regions included in the Zhou 2019 QTL study (which is most of the genome - their study is only missing a few MB at the beginning of several of the chromosomes)
#More details:
##2178 gene symbols didn't map because they are located in regions of the genome not included in the Zhou 2019 QTL study.
##2496 gene symbols didn't map because they lack coordinates in the Rnor6 genome build (!)

#So the cross tab would be:

Intersection_FDR05Adult_QTL_4<-matrix(c(6, 50, 303, 11291), nrow=2, ncol=2)
colnames(Intersection_FDR05Adult_QTL_4)<-c("FDR<0.05", "FDR>0.05")
row.names(Intersection_FDR05Adult_QTL_4)<-c("QTL_LOD>4", "QTL_LOD<4")
Intersection_FDR05Adult_QTL_4

#          FDR<0.05 FDR>0.05
#QTL_LOD>4        6      303
#QTL_LOD<4       50    11291

#Using Fisher's exact test (since one cell has less than 10):
fisher.test(Intersection_FDR05Adult_QTL_4)

# Fisher's Exact Test for Count Data
# 
# data:  Intersection_FDR05Adult_QTL_4
# p-value = 0.003533
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 1.554292 10.543691
# sample estimates:
# odds ratio 
# 4.470441 

#Very pretty! 

#Venn Diagram:

library(limma)

pdf("VennDiagram_AdultMetaAnalysis05_VsQTL_LOD4.pdf", height=4, width=4)
library(eulerr)
v <- euler(c("A"=50, "B"=303, "A&B"=6))
plot(v)
dev.off()
#That seems... unenlightening. But people seem to like them. Why?


#Another option would be plotting the percentages:

#Percent of Adult Meta-Analysis Top Genes (FDR<0.05) that overlap Sig QTLs for Exploratory Locomotion
6/(50+6)
#[1] 0.1071429

#Percent of genes that are not top genes in the adult meta-analysis that overlap sig QTLs for Exploratory Locomotion
296/(11291+303)
#[1] 0.02553045

#That is reasonably striking. A lot more informative than the Venn Diagram.

pdf("BarChart_Overlap_AdultMetaAnalysisFDR05_vs_LocomotorQTL.pdf", height=6, width=3.5)
barplot(height=c((0.1071429)*100, (0.02553045)*100), names.arg=c("FDR<0.05", "FDR>0.05"), xlab="Adult Meta-Analysis", ylab="% Overlapping QTL for Exploratory Locomotion (+/- 1MB)", ylim=c(0,20), col=1)
dev.off()
#I ended up re-doing the axis labels while making the combined figure later vs. messing around in R with the fonts and margins.


#Trying some other versions using less strict cut-offs to see the results:

unique(AdultMeta_CohDandPval_wQTLresults$GeneSymbol[which(AdultMeta_CohDandPval_wQTLresults$BH<0.05 & AdultMeta_CohDandPval_wQTLresults$MaxLOD_1MBUpDown>3)])
#[1] "C2cd3"   "Ezr"     "Hmgn5b"  "Mfge8"   "Pkhd1l1" "Trhr"    "Ucp2"    "Selenop"
#8 genes 


#For comparison:

unique(AdultMeta_CohDandPval_wQTLresults$GeneSymbol[which(AdultMeta_CohDandPval_wQTLresults$MaxLOD_1MBUpDown>3)])
#762 genes


#To make it similar to the PGSEA comparison:
unique(AdultMeta_CohDandPval_wQTLresults$GeneSymbol[which(AdultMeta_CohDandPval_wQTLresults$BH<0.10 & AdultMeta_CohDandPval_wQTLresults$MaxLOD_1MBUpDown>4)])
#10 genes
# [1] "C2cd3"   "Ezr"     "Hmgn5b"  "Olr35"   "Pdzd2"   "Pkhd1l1" "Spcs2"   "Trhr"    "Ttc30a1" "Ucp2"   

unique(AdultMeta_CohDandPval_wQTLresults$GeneSymbol[which(AdultMeta_CohDandPval_wQTLresults$BH<0.10 & AdultMeta_CohDandPval_wQTLresults$MaxLOD_1MBUpDown>0)])
#147 genes

10/147
#[1] 0.06802721
#6.8%
#So definitely more enriched when we use a stricter cut-off for significance in the meta-analysis. Still quite a bit more enriched than the genes that are not sig.


unique(AdultMeta_CohDandPval_wQTLresults$GeneSymbol[which(AdultMeta_CohDandPval_wQTLresults$rawp<0.001 & AdultMeta_CohDandPval_wQTLresults$MaxLOD_1MBUpDown>4)])
#[1] "C2cd3"   "Ezr"     "Hmgn5b"  "Pdzd2"   "Pkhd1l1" "Spcs2"   "Trhr"    "Ucp2"   


#It would be nice to do this for the developmental data too - but if so, we need to use a different cut-off than FDR<0.05 or FDR<0.10 because the results are weaker.




#################
#Left off here and didn't update code/results

# 
# unique(AdultMeta_CohDandPval_wQTLresults$GeneSymbol[which(AdultMeta_CohDandPval_wQTLresults$rawp<0.01 & AdultMeta_CohDandPval_wQTLresults$MaxLOD_1MBUpDown>4)])
# # [1] "Atp6v1c1" "B9d1"     "C2cd3"    "Capn5"    "Cir1"     "Cpq"      "Emc2"     "Eny2"     "Ezr"      "Hmgn5"   
# # [11] "Hmgn5b"   "Olr35"    "Pdk1"     "Pdzd2"    "Pkhd1l1"  "Plekhb1"  "Ppme1"    "Rspo2"    "Serpinh1" "Slc1a3"  
# # [21] "Spcs2"    "Tagap"    "Trhr"     "Ttc30a1"  "Ucp2"  
# #25 total, no mind-blowing new additions (except maybe slc1a3)
# 
# #For comparison:
# unique(AdultMeta_CohDandPval_wQTLresults$GeneSymbol[which(AdultMeta_CohDandPval_wQTLresults$rawp<0.01 & AdultMeta_CohDandPval_wQTLresults$MaxLOD_1MBUpDown<4)])
# #464 genes
# 
# unique(AdultMeta_CohDandPval_wQTLresults$GeneSymbol[which(AdultMeta_CohDandPval_wQTLresults$rawp>0.01 & AdultMeta_CohDandPval_wQTLresults$MaxLOD_1MBUpDown>4)])
# #313 genes
# 
# length(unique(AdultMeta_CohDandPval_wQTLresults$GeneSymbol[AdultMeta_CohDandPval_wQTLresults$MaxLOD_1MBUpDown>0]))
# #[1] 11578
# #Total genes that could have been possibly included in the analysis
# 
# 25/(313+25)
# #[1] 0.0739645
# #7 % of all genes overlapping the >4 LOD areas are p<0.01 in the adult meta-analysis
# 464/(11578-(313+25))
# #[1] 0.04128114
# #whereas 4% of all genes not overlapping >4 LOD areas are p<0.01 in the adult meta-analysis
# 
# 
# 
# AdultMeta_CohDandPval_wQTLresults[which(AdultMeta_CohDandPval_wQTLresults$rawp<0.001 & AdultMeta_CohDandPval_wQTLresults$MaxLOD_1MBUpDown>3), ]
# 
# 
# 
# 
