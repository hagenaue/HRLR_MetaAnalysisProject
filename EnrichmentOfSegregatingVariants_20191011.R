#Quickly Checking whether there is an enrichment of highly-segregating exome single-nucleotide variants in close proximity to our top differentially expressed genes.
##Megan Hagenauer, 10-11-2019


#Setting the working directory:
setwd("~/Documents/Microarray Gen/HRLR/SecondaryAnalyses/WholeExomeSequencing")

#Reading in the data:
F0_Results<-read.csv("Zhifeng_F0Results_newDownload_JL1_forR.csv", header=T, stringsAsFactors=F)
str(F0_Results)

# 'data.frame':	129237 obs. of  15 variables:
# $ ID                       : chr  "1_248147" "1_248273" "1_248307" "1_248495" ...
# $ Chr                      : chr  "1" "1" "1" "1" ...
# $ Post                     : int  248147 248273 248307 248495 248523 248556 315197 316702 316726 392043 ...
# $ Ref                      : chr  "G" "T" "T" "C" ...
# $ Alt                      : chr  "A" "A" "G" "G" ...
# $ HR.AA                    : int  1 11 12 9 11 6 11 0 1 3 ...
# $ HR.AB                    : int  11 1 0 3 1 6 1 0 0 8 ...
# $ HR.BB                    : int  0 0 0 0 0 0 0 1 1 0 ...
# $ LR.AA                    : int  0 12 11 8 11 9 12 2 1 2 ...
# $ LR.AB                    : int  12 0 1 4 1 3 0 0 0 6 ...
# $ LR.BB                    : int  0 0 0 0 0 0 0 0 1 2 ...
# $ Fishers.Exact.Test.Pvalue: num  1 1 1 1 1 ...
# $ X..log10FishersP.        : num  0 0 0 0 0 ...
# $ funcLoc                  : chr  "Intergenic" "Intergenic" "Intergenic" "Intergenic" ...
# $ GeneSym                  : chr  "" "" "" "" ...

#In Excel, I already snooped to see what a genome-wide cut-off would be for alpha:
0.05/129237
#[1] 3.868861e-07

#In the PNAS paper, Zhifeng reported that there were results for only 110172 variants. 
129237-110172
#[1] 19065
#So 19065 of the results were considered unusable. Let's figure out why.

#Perhaps some of these variants are unvarying?  Let's check:

#First, do any of the variants lack reads altogether?

min(apply(F0_Results[,6:11], 1, sum))
#[1] 1

#No, but you can't really claim group differences with only one read, so maybe that is the difference in the totals? Let's check:

table(apply(F0_Results[,6:11], 1, sum))

# 1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22    23    24 
# 1295  1182  1206  1202  1323  1275  1375  1373  1327  1545  1628  1711  1741  1971  2085  2186  2365  2715  3245  4196  5618  9262 17141 60270 

#The number of variants with 13 reads or less:
1295+1182+1206+1202+1323+1275+1375+1373+1327+1545+1628+1711+1741
[1] 18183

#The number of variants with 14 reads or less:
1295+1182+1206+1202+1323+1275+1375+1373+1327+1545+1628+1711+1741+1971 
[1] 20154

#hmmm... that does not exactly add up to 19065. There must have been some other exclusion criteria. Maybe it was based on a minimum number of reads for each genotype?

Total_HR_reads<-apply(F0_Results[,6:8], 1, sum)
Total_LR_reads<-apply(F0_Results[,9:11], 1, sum)

table(Total_HR_reads, Total_LR_reads)

# 0     1     2     3     4     5     6     7     8     9    10    11    12
# 0      0   622   248   143    67    40    22    15     5    11     5     6     3
# 1    673   629   444   294   190   113    60    41    23    15    11    14     5
# 2    305   465   446   391   242   210   135    77    51    35    23    17    14
# 3    154   326   423   421   345   247   163   124    90    59    37    19    12
# 4     69   225   296   383   381   317   265   211   149    75    53    26    20
# 5     54   153   260   312   349   417   354   285   215   145    95    54    21
# 6     28    88   181   233   363   403   421   385   339   210   155    78    44
# 7     14    61   116   198   295   412   466   484   462   348   224   122    58
# 8     10    33    62   150   217   298   463   562   578   522   427   269   144
# 9      5    37    49    94   157   280   424   559   722   769   754   544   349
# 10     8    19    25    62   123   203   295   519   783  1104  1334  1077   742
# 11     5     5    18    41    67   126   206   404   785  1488  2419  2956  3061
# 12     4     6    10    24    51    73   166   275   686  1773  5564 14080 60270

sum(Total_HR_reads<6)
#[1] 13714
sum(Total_LR_reads<6)
#[1] 15308

sum((Total_HR_reads<6)|(Total_LR_reads<6))
#[1] 18368

sum((Total_HR_reads<7)|(Total_LR_reads<7))
#[1] 22020
#Hmmm.... fail. Still not explaining the 19065 difference.  Oh well, it is not completely essential for my analysis that I match what they did - this is going to be a pretty blunt estimate anyway.

#In the PNAS paper, Zhifeng reported that 1584 of the variants were fully segregating. From glancing at the spreadsheet, this is difficult to determine using the p-values, because some of the fully segregating variants have p-values that are weaker due to smaller number of samples sequenced. Let's see what the approximate p-value cut-off would be to get this number of variants:

#First, a bonferonni-corrected alpha of 0.05:
sum(F0_Results$Fishers.Exact.Test.Pvalue<3.868861e-07)
#[1] 3344
3344/129237
#[1] 0.02587494 - so approximately 2.6% of the genome is segregated by this definition.

sum(F0_Results$Fishers.Exact.Test.Pvalue<3.868861e-08)
#[1] 2334
sum(F0_Results$Fishers.Exact.Test.Pvalue<3.868861e-09)
#[1] 1326

sum(F0_Results$Fishers.Exact.Test.Pvalue<1e-08)
#[1] 1645
#Pretty darn close.  Looking at the results, this p-value seems to solidly catch almost all of the fully segregating variants, although there are still a few missed and a few non fully segregating variants mixed in.


*********************************
  
#Speaking of which, since Zhifeng's study used Rnor5 coordinates, I'm going to need the Rnor5 coordinates for all of the genes included in our differential expression meta-analysis. 
#I created this as part of my analysis looking for enrichment in the QTL results. Let me find it again:

setwd("~/Documents/Microarray Gen/HRLR/SecondaryAnalyses/RatGenomeDatabase_eQTL")

EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5<-read.csv("EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_wQTLresults.csv", header=T, stringsAsFactors = F)
str(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5)

# 'data.frame':	49084 obs. of  17 variables:
#   $ X                    : int  42891 27779 44775 32800 44776 27787 44837 27788 25561 45751 ...
# $ EntrezGeneID         : int  502213 292438 679691 308257 679691 292461 679825 292462 25604 683983 ...
# $ GeneSymbol           : chr  "Vom2r3" "Vom2r2" "Vom2r5" "Vom2r6" ...
# $ Chromosome           : chr  "1" "1" "1" "1" ...
# $ ChromosomeLocation   : int  396699 400255 699412 744615 804861 1101686 1101964 1702695 1734862 1771720 ...
# $ feat_name            : chr  "1_396700" "1_400256" "1_699413" "1_744616" ...
# $ Rnor6_Chr            : chr  "1" "1" "1" "1" ...
# $ Rnor6_Start          : int  396700 400256 699413 744616 804862 1101687 1101965 1702696 1734863 1771721 ...
# $ Rnor6_Stop           : int  396799 400355 699512 744715 804961 1101786 1102064 1702795 1734962 1771820 ...
# $ Rnor5_Chr            : chr  "1" "1" "1" "1" ...
# $ Rnor5_Start          : int  388173 391729 688298 733501 792899 2810605 2810883 3400376 3432543 3468471 ...
# $ Rnor5_Stop           : int  388272 391828 688397 733600 792998 2810704 2810982 3400475 3432642 3468570 ...
# $ Rnor5_Start_Round    : int  0 0 1000000 1000000 1000000 3000000 3000000 3000000 3000000 3000000 ...
# $ Feat_name_Rnor5_Start: chr  "1_0" "1_0" "1_1000000" "1_1000000" ...
# $ Loci                 : chr  NA NA NA NA ...
# $ LOD..H.K.            : num  NA NA NA NA NA NA NA NA NA NA ...
# $ MaxLOD_1MBUpDown     : num  NA NA NA NA NA NA NA NA NA NA ...

#Alright, let's make a new column for -1 MB from start and +1 MB from stop:

#Double-checking gene length first:
summary(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5$Rnor5_Stop-EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5$Rnor5_Start)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   49.00   99.00   99.00   98.94   99.00  164.00   30728 

#Yeah -  I artificially set the gene length to be 100 kb and defined the coordinates by their start when using the conversion tool (Rnor5 vs Rnor6) because I discovered that the tool was having some difficulty figuring out what was analagous between assemblies for larger stretches.
#This is definitely going to be a blunt/imperfect analysis.

EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5$Rnor5_StartMinus1MB<-EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5$Rnor5_Start-1000000

EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5$Rnor5_StartPlus1MB<-EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5$Rnor5_Start+1000000

head(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5)

# X EntrezGeneID GeneSymbol Chromosome ChromosomeLocation   feat_name Rnor6_Chr Rnor6_Start Rnor6_Stop Rnor5_Chr Rnor5_Start Rnor5_Stop Rnor5_Start_Round
# 1 1    100008565   Slc39a4l       <NA>                 NA        <NA>      <NA>          NA         NA      <NA>          NA         NA                NA
# 2 2    100034253      Gnl3l          X           20037557  X_20037558         X    20037558   20037657         X    20785115   20785214          21000000
# 3 3    100036582    Olr1867          8           20565731  8_20565732         8    20565732   20565831         8    20639977   20640076          21000000
# 4 4    100036765     Ccdc92         12           37211315 12_37211316        12    37211316   37211415        12    39085314   39085413          39000000
# 5 5    100049583      Trex1          8          117796126 8_117796127         8   117796127  117796226         8   117147872  117147971         117000000
# 6 6    100101342 Vom2r-ps11       <NA>                 NA        <NA>      <NA>          NA         NA      <NA>          NA         NA                NA
# Feat_name_Rnor5_Start Rnor5_StartMinus1MB Rnor5_StartPlus1MB
# 1                  <NA>                  NA                 NA
# 2            X_21000000            19785115           21785115
# 3            8_21000000            19639977           21639977
# 4           12_39000000            38085314           40085314
# 5           8_117000000           116147872          118147872
# 6                  <NA>                  NA                 NA

#Looks good.

#Alright, let's figure out the maximum segregation in any particular stretch:

Min_Pval<-vector(mode="numeric", length=nrow(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5))

# nrow(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5)
# 
# 
# for(i in c(1:50)){
#   
#   if((is.na(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5$Chromosome[i])==T)){Min_Pval[i]<-NA}else{
#   
#     if((sum(F0_Results$Chr==EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5$Chromosome[i]))<1){Min_Pval[i]<-NA}else{
#     
#     Relevant_F0Results<-F0_Results[F0_Results$Chr==EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5$Chromosome[i],]
#     
#     if(sum((Relevant_F0Results$Post>EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5$Rnor5_StartMinus1MB[i])&(Relevant_F0Results$Post<EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5$Rnor5_StartPlus1MB[i])==T)<1){Min_Pval[i]<-NA}else{
#     
#       Relevant_F0Results_PlusMinus1MB<-Relevant_F0Results[(Relevant_F0Results$Post>EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5$Rnor5_StartMinus1MB[i])&(Relevant_F0Results$Post<EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5$Rnor5_StartPlus1MB[i]),]
#       
#       if(is.na(Relevant_F0Results_PlusMinus1MB$Fishers.Exact.Test.Pvalue)==T){Min_Pval[i]<-NA}else{
#         
#         Min_Pval[i]<-min(Relevant_F0Results_PlusMinus1MB$Fishers.Exact.Test.Pvalue, na.rm=T)
#     }}}}
#   
#   rm(Relevant_F0Results, Relevant_F0Results_PlusMinus1MB)
# }
# 
# 
# head(Min_Pval)
# 
# summary(Min_Pval)
# 
# EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5[c(200:250),]


#I'm getting a ton of error messages due to missing data. I'm just going to clean up the reference database first instead of embedding so many if/else statements.
length(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5$Chromosome)
[1] 49084

sum(is.na(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5$Chromosome))
[1] 30526

table(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5$Chromosome)

# 1 1_AABR07046231v1_random     1_KL567891v1_random                      10    10_KL568017v1_random                      11 
# 2533                       1                       2                    1524                       1                     493 
# 12                      13                      14                      15                      16                      17 
# 512                     579                     559                     562                     490                     505 
# 18                      19                       2     2_KL567908v1_random                      20    20_KL568103v1_random 
# 416                     454                    1137                       1                     594                       1 
# 3                       4     4_KL567939v1_random                       5                       6                       7 
# 1557                    1162                      12                    1096                     744                    1177 
# 7_AABR07046136v1_random     7_KL567988v1_random                       8                       9       Un_AABR07024031v1       Un_AABR07024040v1 
# 1                       1                    1072                     677                       1                       1 
# Un_AABR07024041v1       Un_AABR07024044v1       Un_AABR07024104v1       Un_AABR07024106v1       Un_AABR07024120v1       Un_AABR07024203v1 
# 3                       1                       5                       1                       1                       1 
# Un_AABR07024291v1       Un_AABR07024382v1       Un_AABR07024421v1           Un_KL568398v1           Un_KL568405v1           Un_KL568409v1 
# 3                       1                       1                       1                       4                       7 
# Un_KL568414v1           Un_KL568418v1           Un_KL568438v1           Un_KL568439v1           Un_KL568447v1           Un_KL568458v1 
# 3                       1                       2                       1                       5                       2 
# Un_KL568460v1           Un_KL568463v1           Un_KL568468v1           Un_KL568472v1           Un_KL568473v1                       X 
# 1                       1                       1                       1                       2                     636 
# X_KL568122v1_random     X_KL568132v1_random                       Y     Y_KL568139v1_random 
# 1                       1                       5                       1 

sum(table(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5$Chromosome))
#[1] 18558

18558+30526
#[1] 49084

EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr<-EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5[EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5$Chromosome%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21", "22", "X","Y"),]

table(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr$Chromosome)
# 1   10   11   12   13   14   15   16   17   18   19    2   20    3    4    5    6    7    8    9    X    Y 
# 2533 1524  493  512  579  559  562  490  505  416  454 1137  594 1557 1162 1096  744 1177 1072  677  636    5 

sum(table(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr$Chromosome))
#[1] 18484

length(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr$Chromosome)
#[1] 18484

#...so that will get rid of some of the errors.


Min_Pval<-vector(mode="numeric", length=nrow(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2))

nrow(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr)

table(F0_Results$Chr)
# 1    10    11    12    13    14    15    16    17    18    19     2    20     3     4     5     6     7     8     9     M     X 
# 16164 10041  3222  4038  3629  4965  4034  4085  3453  3477  3415  8425  5902 10760  7281  7637  5829  7918  6888  5064   109  2901 


F0_Results_StandardChr<-F0_Results[F0_Results$Chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21", "22", "X","Y"),]

table(F0_Results_StandardChr$Chr)
#...that will get rid of some more errors.

#How many still don't have Rnor5 coordinates?
sum(is.na(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr$Rnor5_Start))
#[1] 148

EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2<-EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr[is.na(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr$Rnor5_Start)==F,]


for(i in c(1:nrow(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2))){
  
      Relevant_F0Results<-F0_Results_StandardChr[F0_Results_StandardChr$Chr==EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Chromosome[i],]
      
      if(sum((Relevant_F0Results$Post>EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_StartMinus1MB[i])&(Relevant_F0Results$Post<EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_StartPlus1MB[i])==T)<1){Min_Pval[i]<-NA}else{
        
        Relevant_F0Results_PlusMinus1MB<-Relevant_F0Results[(Relevant_F0Results$Post>EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_StartMinus1MB[i])&(Relevant_F0Results$Post<EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_StartPlus1MB[i]),]
        
        if(is.na(Relevant_F0Results_PlusMinus1MB$Fishers.Exact.Test.Pvalue)==T){Min_Pval[i]<-NA}else{
          
          Min_Pval[i]<-min(Relevant_F0Results_PlusMinus1MB$Fishers.Exact.Test.Pvalue, na.rm=T)
        }}
  
  rm(Relevant_F0Results, Relevant_F0Results_PlusMinus1MB)
}

summary(Min_Pval)

# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# 0.000000 0.000000 0.000013 0.018121 0.001559 1.000000       27 

#So the vast majority of genes are within a segregating variant, but not necessarily a significantly segregating variant.

sum(Min_Pval<3.868861e-07, na.rm=T)
#[1] 6556

sum(Min_Pval>3.868861e-07, na.rm=T)
#[1] 11901

6556/11901
#[1] 0.5508781
#More than half of the genes in the genome are within a significantly segregating variant.
#How can that be if the median is greater than that? I guess the median is somehow influenced by NAs? If so, interesting. 

length(Min_Pval)
18484
6556+11901
18457
18457+27
18484

#Maybe I should use a stricter cut-off?
sum(Min_Pval<3.868861e-08, na.rm=T)
#[1] 5416
#Still huge percentage.
sum(Min_Pval<3.868861e-09, na.rm=T)
#[1] 3960
sum(Min_Pval<3.868861e-10, na.rm=T)
#[1] 2629
sum(Min_Pval<3.868861e-11, na.rm=T)
#[1] 1870
sum(Min_Pval<3.868861e-12, na.rm=T)
#[1] 1256


#Interesting.

EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$MinPval<-Min_Pval

str(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2)

#Alright, let's determine how that aligns with the differential expression results:


setwd("~/Documents/Microarray Gen/HRLR/ThesisMetaAnalysisOutput")

#Note: I used the version of the results that has the genename/dates confusion problem fixed... but if you open any of the files in Excel the problem simply recreates itself.
AdultMeta_CohDandPval<-read.csv("AdultMeta_CohDandPval_GeneDatesFixed.csv", header=T, stringsAsFactors = F)
head(AdultMeta_CohDandPval)

# X      rawp        BH BY GeneSymbol   estimate        SE      pval      CI_lb     CI_ub datasets d.F43 var.d.F43 g.F43 var.g.F43 d.F37 var.d.F37 g.F37
# 1 1 0.2408056 0.6953615  1     March1  0.6534483 0.5570860 0.2408056 -0.4384202 1.7453168        2  0.97      0.45  0.88      0.36    NA        NA    NA
# 2 2 0.5138646 0.8588518  1      Sept1  0.3658065 0.5603378 0.5138646 -0.7324355 1.4640484        2 -0.06      0.40 -0.05      0.33    NA        NA    NA
# 3 3 0.1215975 0.5686317  1    March11  0.8700617 0.5620180 0.1215975 -0.2314732 1.9715967        2  0.74      0.43  0.67      0.35    NA        NA    NA
# 4 4 0.8140614 0.9569539  1     Sept11  0.1257143 0.5345225 0.8140614 -0.9219305 1.1733591        2  0.20      0.40  0.18      0.33    NA        NA    NA
# 5 5 0.6948233 0.9212827  1     Sept15 -0.2100000 0.5352801 0.6948233 -1.2591298 0.8391298        2 -0.21      0.40 -0.19      0.33    NA        NA    NA
# 6 6 0.2808528 0.7250387  1     March2  0.5961111 0.5527708 0.2808528 -0.4872997 1.6795220        2  0.92      0.44  0.83      0.36    NA        NA    NA
# var.g.F37 d.NC.F34 var.d.NC.F34 g.NC.F34 var.g.NC.F34 d.Affy.F4 var.d.Affy.F4 g.Affy.F4 var.g.Affy.F4 d.RNAseq.F29 var.d.RNAseq.F29 g.RNAseq.F29
# 1        NA       NA           NA       NA           NA        NA            NA        NA            NA        -0.05             1.00         0.03
# 2        NA       NA           NA       NA           NA        NA            NA        NA            NA         1.92             1.46        -1.10
# 3        NA       NA           NA       NA           NA        NA            NA        NA            NA         1.23             1.19        -0.70
# 4        NA       NA           NA       NA           NA        NA            NA        NA            NA        -0.06             1.00         0.03
# 5        NA       NA           NA       NA           NA        NA            NA        NA            NA        -0.21             1.01         0.12
# 6        NA       NA           NA       NA           NA        NA            NA        NA            NA        -0.14             1.00         0.08
# var.g.RNAseq.F29
# 1             0.33
# 2             0.48
# 3             0.39
# 4             0.33
# 5             0.33
# 6             0.33


AdultMeta_CohDandPval_wExomeSeq_results<-join(AdultMeta_CohDandPval, EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2, by="GeneSymbol", type="left")

head(AdultMeta_CohDandPval_wExomeSeq_results)

write.csv(AdultMeta_CohDandPval_wExomeSeq_results, "AdultMeta_CohDandPval_wExomeSeq_results")


unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-08)])

# [1] "Aar2"     "Acox3"    "Acss2"    "Blvra"    "C2cd3"    "Car9"     "Ddx20"    "Dnaaf3"   "Egfem1"   "Epb41l4a" "Etv4"     "Exosc7"   "Ezr"      "Hmgn5b"  
# [15] "Kif15"    "Mfge8"    "Mki67"    "Mvb12b"   "Myh6"     "Oard1"    "Pkhd1l1"  "Pkib"     "Robo3"    "Rpl17"    "Samd5"    "Slc27a1"  "Tdg"      "Tmem144"
# [29] "Tnnt1"    "Trappc6a" "Trhr"     "Trmt10a"  "Tubg1"    "Ucp2"     "Zfp90"  


length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-07)]))
#[1] 42

length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval>3.868861e-07)]))
#[1] 19

42/(42+19)
#[1] 0.6885246

(42+19)
#[1] 61
#Out of 74 in the full dataset - so 13 are unable to be aligned with the exome sequencing results for some reason. 


length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH>0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-07)]))
#[1] 4855

length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH>0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval>3.868861e-07)]))
#[1] 8779

4855/(4855+8779)
#[1] 0.3560951

(4855+8779)
#[1] 13634

#Wow - that actually is an enrichment. Interesting.
#I should probably run some stats to back that up.

Intersection_FDR05Adult_SigSNP<-matrix(c(42, 19, 4855, 8779), nrow=2, ncol=2)
colnames(Intersection_FDR05Adult_SigSNP)<-c("FDR<0.05", "FDR>0.05")
row.names(Intersection_FDR05Adult_SigSNP)<-c("SigSNP", "NonSigSNP")
Intersection_FDR05Adult_SigSNP

#             FDR<0.05 FDR>0.05
# SigSNP          42     4855
# NonSigSNP       19     8779

#Using Fisher's exact test (since one cell has less than 10):
fisher.test(Intersection_FDR05Adult_SigSNP)

# Fisher's Exact Test for Count Data
# 
# data:  Intersection_FDR05Adult_SigSNP
# p-value = 2.452e-07
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.270782 7.284713
# sample estimates:
# odds ratio 
#   3.996654 

#Yep, pretty darn sig.

#I would actually like to see the results for a more significant cut-off - something closer to fully segregating.

length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-10)]))
#[1] 21

#Hey look - it is some of our favorite people:
# [1] "Aar2"     "Acox3"    "Acss2"    "C2cd3"    "Car9"     "Egfem1"   "Epb41l4a" "Etv4"     "Ezr"      "Mfge8"    "Mki67"    "Mvb12b"   "Oard1"    "Pkhd1l1" 
# [15] "Pkib"     "Tdg"      "Tmem144"  "Trappc6a" "Trhr"     "Tubg1"    "Ucp2" 

length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval>3.868861e-10)]))
#[1] 40

21/(21+40)
#[1] 0.3442623

length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH>0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-10)]))
#[1] 1971

length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH>0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval>3.868861e-10)]))
#[1] 11632

1971/(1971+11632)
#[1] 0.1448945

#That's pretty.


Intersection_FDR05Adult_SigSNP<-matrix(c(21, 40, 1971, 11632), nrow=2, ncol=2)
colnames(Intersection_FDR05Adult_SigSNP)<-c("FDR<0.05", "FDR>0.05")
row.names(Intersection_FDR05Adult_SigSNP)<-c("SigSNP", "NonSigSNP")
Intersection_FDR05Adult_SigSNP
#           FDR<0.05 FDR>0.05
# SigSNP          21     1971
# NonSigSNP       40    11632

#Using Fisher's exact test (since one cell has less than 10):
fisher.test(Intersection_FDR05Adult_SigSNP)

# 
# Fisher's Exact Test for Count Data
# 
# data:  Intersection_FDR05Adult_SigSNP
# p-value = 8.397e-05
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 1.731466 5.394157
# sample estimates:
# odds ratio 
# 3.097964 

#Interesting - even though the enrichment is greater, the effect is less significant because of the subgroup sample sizes.

#Oooh... wouldn't it be cool to plot the enrichment based on distance from the DE gene and significance of effect?   Let's give it a try.


#first, let's save the results that we have as a new vector:


EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Min_Pval_PlusMinus1MB<-Min_Pval

Min_Pval<-vector(mode="numeric", length=nrow(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2))


EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_StartMinus05MB<-EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_Start-500000

EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_StartPlus05MB<-EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_Start+500000


for(i in c(1:nrow(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2))){
  
  Relevant_F0Results<-F0_Results_StandardChr[F0_Results_StandardChr$Chr==EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Chromosome[i],]
  
  if(sum((Relevant_F0Results$Post>EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_StartMinus05MB[i])&(Relevant_F0Results$Post<EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_StartPlus05MB[i])==T)<1){Min_Pval[i]<-NA}else{
    
    Relevant_F0Results_PlusMinus05MB<-Relevant_F0Results[(Relevant_F0Results$Post>EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_StartMinus05MB[i])&(Relevant_F0Results$Post<EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_StartPlus05MB[i]),]
    
    if(is.na(Relevant_F0Results_PlusMinus05MB$Fishers.Exact.Test.Pvalue)==T){Min_Pval[i]<-NA}else{
      
      Min_Pval[i]<-min(Relevant_F0Results_PlusMinus05MB$Fishers.Exact.Test.Pvalue, na.rm=T)
    }}
  
  rm(Relevant_F0Results, Relevant_F0Results_PlusMinus05MB)
}

summary(Min_Pval)
# 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.00000 0.00000 0.00022 0.03890 0.01302 1.00000     133 

#I'm going to overwrite the old MinPval so that I can re-use the code above.
EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$MinPval<-Min_Pval

AdultMeta_CohDandPval_wExomeSeq_results<-join(AdultMeta_CohDandPval, EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2, by="GeneSymbol", type="left")

head(AdultMeta_CohDandPval_wExomeSeq_results)

setwd("~/Documents/Microarray Gen/HRLR/SecondaryAnalyses/WholeExomeSequencing")

write.csv(AdultMeta_CohDandPval_wExomeSeq_results, "AdultMeta_CohDandPval_wExomeSeq_results")



unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-07)])

# [1] "Aar2"     "Acox3"    "Acss2"    "Blvra"    "C2cd3"    "Car9"     "Cav1"     "Ddx20"    "Dnaaf3"   "Egfem1"   "Epb41l4a" "Etv4"     "Exosc7"   "Ezr"     
# [15] "Fn3k"     "Kif15"    "Mfge8"    "Mki67"    "Mvb12b"   "Myh6"     "Oard1"    "Pkhd1l1"  "Pkib"     "R3hdm4"   "Rpl17"    "Samd5"    "Tdg"      "Tmem144" 
# [29] "Tnnt1"    "Trappc6a" "Trhr"     "Trmt10a"  "Tubg1"    "Ucp2"  


length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-07)]))
#[1] 34

length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval>3.868861e-07)]))
#[1] 27

34/(34+27)
#[1] 0.557377

(34+27)
#[1] 61
#Out of 74 in the full dataset - so 13 are unable to be aligned with the exome sequencing results for some reason. 


length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH>0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-07)]))
#[1] 3576

length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH>0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval>3.868861e-07)]))
#[1] 10008

3576/(3576+10008)
#[1] 0.2632509

(3576+10008)
#[1] 13584

#How pretty.

Intersection_FDR05Adult_SigSNP<-matrix(c(34, 27, 3576, 10008), nrow=2, ncol=2)
colnames(Intersection_FDR05Adult_SigSNP)<-c("FDR<0.05", "FDR>0.05")
row.names(Intersection_FDR05Adult_SigSNP)<-c("SigSNP", "NonSigSNP")
Intersection_FDR05Adult_SigSNP

#             FDR<0.05 FDR>0.05
# SigSNP          34     3576
# NonSigSNP       27    10008

#Using Fisher's exact test:
fisher.test(Intersection_FDR05Adult_SigSNP)

Fisher's Exact Test for Count Data

# Fisher's Exact Test for Count Data
# 
# data:  Intersection_FDR05Adult_SigSNP
# p-value = 1.342e-06
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   2.060490 6.081352
# sample estimates:
#   odds ratio 
# 3.523837 


#To parallel the last analysis, let's do it using a more stringent cut off for segregation:

unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-10)])

# [1] "Aar2"     "Acox3"    "Acss2"    "C2cd3"    "Car9"     "Epb41l4a" "Etv4"     "Ezr"      "Mki67"    "Mvb12b"   "Oard1"    "Pkhd1l1"  "Pkib"     "Tdg"     
#[15] "Tmem144"  "Trappc6a" "Trhr"     "Tubg1"    "Ucp2"


length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-10)]))
#[1] 19

length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval>3.868861e-10)]))
#[1] 42

19/(19+42)
#[1] 0.3114754


length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH>0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-10)]))
#[1] 1330

length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH>0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval>3.868861e-10)]))
#[1] 12234

1330/(1330+12234)
#[1] 0.09805367

#Wow - that enrichment is even prettier. The more strict we go, the more enriched the results for the DE genes. Neat.

Intersection_FDR05Adult_SigSNP<-matrix(c(19, 42, 1330, 12234), nrow=2, ncol=2)
colnames(Intersection_FDR05Adult_SigSNP)<-c("FDR<0.05", "FDR>0.05")
row.names(Intersection_FDR05Adult_SigSNP)<-c("SigSNP", "NonSigSNP")
Intersection_FDR05Adult_SigSNP

#             FDR<0.05 FDR>0.05
# SigSNP          19     1330
# NonSigSNP       42    12234

#Using Fisher's exact test:
fisher.test(Intersection_FDR05Adult_SigSNP)

# Fisher's Exact Test for Count Data
# 
# data:  Intersection_FDR05Adult_SigSNP
# p-value = 3.728e-06
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.278221 7.337667
# sample estimates:
# odds ratio 
#   4.160875 


#Fun fun.  Let's see what happens if we make it even tighter to the segregating SNPs:


EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Min_Pval_PlusMinus05MB<-Min_Pval

Min_Pval<-vector(mode="numeric", length=nrow(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2))


EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_StartMinus025MB<-EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_Start-250000

EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_StartPlus025MB<-EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_Start+250000


for(i in c(1:nrow(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2))){
  
  Relevant_F0Results<-F0_Results_StandardChr[F0_Results_StandardChr$Chr==EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Chromosome[i],]
  
  if(sum((Relevant_F0Results$Post>EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_StartMinus025MB[i])&(Relevant_F0Results$Post<EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_StartPlus025MB[i])==T)<1){Min_Pval[i]<-NA}else{
    
    Relevant_F0Results_PlusMinus025MB<-Relevant_F0Results[(Relevant_F0Results$Post>EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_StartMinus025MB[i])&(Relevant_F0Results$Post<EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_StartPlus025MB[i]),]
    
    if(is.na(Relevant_F0Results_PlusMinus025MB$Fishers.Exact.Test.Pvalue)==T){Min_Pval[i]<-NA}else{
      
      Min_Pval[i]<-min(Relevant_F0Results_PlusMinus025MB$Fishers.Exact.Test.Pvalue, na.rm=T)
    }}
  
  rm(Relevant_F0Results, Relevant_F0Results_PlusMinus025MB)
}

summary(Min_Pval)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.00000 0.00000 0.00133 0.07042 0.05000 1.00000     297 


#I'm going to overwrite the old MinPval so that I can re-use the code above.
EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$MinPval<-Min_Pval

AdultMeta_CohDandPval_wExomeSeq_results<-join(AdultMeta_CohDandPval, EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2, by="GeneSymbol", type="left")

head(AdultMeta_CohDandPval_wExomeSeq_results)

setwd("~/Documents/Microarray Gen/HRLR/SecondaryAnalyses/WholeExomeSequencing")

write.csv(AdultMeta_CohDandPval_wExomeSeq_results, "AdultMeta_CohDandPval_wExomeSeq_results")



unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-07)])

# [1] "Aar2"     "Acox3"    "Acss2"    "Blvra"    "C2cd3"    "Car9"     "Cav1"     "Dnaaf3"   "Epb41l4a" "Etv4"     "Exosc7"   "Ezr"      "Fn3k"     "Kif15"   
# [15] "Mfge8"    "Mki67"    "Mvb12b"   "Myh6"     "Oard1"    "Pkhd1l1"  "Pkib"     "Rpl17"    "Samd5"    "Tdg"      "Tmem144"  "Tnnt1"    "Trappc6a" "Trhr"    
# [29] "Trmt10a"  "Tubg1"    "Ucp2"  


length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-07)]))
#[1] 31

length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval>3.868861e-07)]))
#[1] 29

31/(31+29)
#[1] 0.5166667


length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH>0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-07)]))
#[1] 2679

length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH>0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval>3.868861e-07)]))
#[1] 10841

2679/(2679+10841)
#[1] 0.1981509

#Really pretty.

Intersection_FDR05Adult_SigSNP<-matrix(c(31, 29, 2679, 10841), nrow=2, ncol=2)
colnames(Intersection_FDR05Adult_SigSNP)<-c("FDR<0.05", "FDR>0.05")
row.names(Intersection_FDR05Adult_SigSNP)<-c("SigSNP", "NonSigSNP")
Intersection_FDR05Adult_SigSNP

#             FDR<0.05 FDR>0.05
# SigSNP          31     2679
# NonSigSNP       29    10841

#Using Fisher's exact test:
fisher.test(Intersection_FDR05Adult_SigSNP)

# data:  Intersection_FDR05Adult_SigSNP
# p-value = 4.275e-08
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.516636 7.452410
# sample estimates:
# odds ratio 
#   4.325104 



unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-10)])

#  [1] "Aar2"     "Acox3"    "Acss2"    "C2cd3"    "Car9"     "Epb41l4a" "Mki67"    "Mvb12b"   "Oard1"    "Pkhd1l1"  "Pkib"     "Tdg"      "Tmem144"  "Trappc6a"
#[15] "Trhr"     "Tubg1"    "Ucp2"    


length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-10)]))
#[1] 17

length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval>3.868861e-10)]))
#[1] 43

17/(17+43)
#[1] 0.2833333


length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH>0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-10)]))
#[1] 920

length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH>0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval>3.868861e-10)]))
#[1] 12560

920/(920+12560)
#[1] 0.06824926

#Really pretty.

Intersection_FDR05Adult_SigSNP<-matrix(c(17, 43, 920, 12560), nrow=2, ncol=2)
colnames(Intersection_FDR05Adult_SigSNP)<-c("FDR<0.05", "FDR>0.05")
row.names(Intersection_FDR05Adult_SigSNP)<-c("SigSNP", "NonSigSNP")
Intersection_FDR05Adult_SigSNP

#             FDR<0.05 FDR>0.05
# SigSNP          17      920
# NonSigSNP       43    12560

#Using Fisher's exact test:
fisher.test(Intersection_FDR05Adult_SigSNP)

# Fisher's Exact Test for Count Data
# 
# data:  Intersection_FDR05Adult_SigSNP
# p-value = 3.76e-07
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   2.873351 9.703315
# sample estimates:
#   odds ratio 
# 5.395989 


#Ok, one more time.... getting even closer to the segregating SNP!



EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Min_Pval_PlusMinus025MB<-Min_Pval

Min_Pval<-vector(mode="numeric", length=nrow(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2))


EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_StartMinus010MB<-EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_Start-100000

EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_StartPlus010MB<-EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_Start+100000


for(i in c(1:nrow(EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2))){
  
  Relevant_F0Results<-F0_Results_StandardChr[F0_Results_StandardChr$Chr==EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Chromosome[i],]
  
  if(sum((Relevant_F0Results$Post>EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_StartMinus010MB[i])&(Relevant_F0Results$Post<EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_StartPlus010MB[i])==T)<1){Min_Pval[i]<-NA}else{
    
    Relevant_F0Results_PlusMinus010MB<-Relevant_F0Results[(Relevant_F0Results$Post>EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_StartMinus010MB[i])&(Relevant_F0Results$Post<EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$Rnor5_StartPlus010MB[i]),]
    
    if(is.na(Relevant_F0Results_PlusMinus010MB$Fishers.Exact.Test.Pvalue)==T){Min_Pval[i]<-NA}else{
      
      Min_Pval[i]<-min(Relevant_F0Results_PlusMinus010MB$Fishers.Exact.Test.Pvalue, na.rm=T)
    }}
  
  rm(Relevant_F0Results, Relevant_F0Results_PlusMinus010MB)
}

summary(Min_Pval)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  0.0001  0.0094  0.1323  0.1538  1.0000     789 


#I'm going to overwrite the old MinPval so that I can re-use the code above.
EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2$MinPval<-Min_Pval

AdultMeta_CohDandPval_wExomeSeq_results<-join(AdultMeta_CohDandPval, EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_StandardChr2, by="GeneSymbol", type="left")

head(AdultMeta_CohDandPval_wExomeSeq_results)

setwd("~/Documents/Microarray Gen/HRLR/SecondaryAnalyses/WholeExomeSequencing")

write.csv(AdultMeta_CohDandPval_wExomeSeq_results, "AdultMeta_CohDandPval_wExomeSeq_results.csv")



unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-07)])

# [1] "Aar2"     "Acox3"    "Acss2"    "Blvra"    "C2cd3"    "Car9"     "Dnaaf3"   "Epb41l4a" "Etv4"     "Exosc7"   "Kif15"    "Mfge8"    "Mki67"    "Mvb12b"  
# [15] "Myh6"     "Oard1"    "Pkhd1l1"  "Rpl17"    "Tdg"      "Tnnt1"    "Trappc6a" "Trmt10a"  "Tubg1"    "Ucp2"  


length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-07)]))
#[1] 24

length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval>3.868861e-07)]))
#[1] 36

24/(24+36)
#[1] 0.4


length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH>0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-07)]))
#[1] 1975

length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH>0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval>3.868861e-07)]))
#[1]  11321

1975/(1975+11321)
#[1] 0.1485409

#Really pretty.

Intersection_FDR05Adult_SigSNP<-matrix(c(24, 36, 1975, 11321), nrow=2, ncol=2)
colnames(Intersection_FDR05Adult_SigSNP)<-c("FDR<0.05", "FDR>0.05")
row.names(Intersection_FDR05Adult_SigSNP)<-c("SigSNP", "NonSigSNP")
Intersection_FDR05Adult_SigSNP

#             FDR<0.05 FDR>0.05
# SigSNP          24     1975
# NonSigSNP       36    11321

#Using Fisher's exact test:
fisher.test(Intersection_FDR05Adult_SigSNP)

# Fisher's Exact Test for Count Data
# 
# data:  Intersection_FDR05Adult_SigSNP
# p-value = 2.096e-06
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 2.175627 6.599934
# sample estimates:
# odds ratio 
# 3.820938 






unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-10)])

#  [1] "Acox3"    "Acss2"    "C2cd3"    "Car9"     "Epb41l4a" "Mki67"    "Mvb12b"   "Oard1"    "Pkhd1l1"  "Tdg"      "Trappc6a" "Tubg1"    "Ucp2"    


length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-10)]))
#[1] 13

length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH<0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval>3.868861e-10)]))
#[1] 47

13/(13+47)
#[1] 0.2166667


length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH>0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval<3.868861e-10)]))
#[1] 618

length(unique(AdultMeta_CohDandPval_wExomeSeq_results$GeneSymbol[which(AdultMeta_CohDandPval_wExomeSeq_results$BH>0.05 & AdultMeta_CohDandPval_wExomeSeq_results$MinPval>3.868861e-10)]))
#[1] 12653

618/(618+12653)
#[1] 0.0465677

#Really pretty.

Intersection_FDR05Adult_SigSNP<-matrix(c(13, 47, 618, 12653), nrow=2, ncol=2)
colnames(Intersection_FDR05Adult_SigSNP)<-c("FDR<0.05", "FDR>0.05")
row.names(Intersection_FDR05Adult_SigSNP)<-c("SigSNP", "NonSigSNP")
Intersection_FDR05Adult_SigSNP

#             FDR<0.05 FDR>0.05
# SigSNP          13      618
# NonSigSNP       47    12653

#Using Fisher's exact test:
fisher.test(Intersection_FDR05Adult_SigSNP)

# Fisher's Exact Test for Count Data

# data:  Intersection_FDR05Adult_SigSNP
# p-value = 3.497e-06
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   2.794629 10.708419
# sample estimates:
#   odds ratio 
# 5.660698 


#Well that was pretty cool.  I'm willing to bet that there was a better way to analyze that than running a bunch of individual Fisher's Exact tests, but it definitely shows a clear pattern. 
#Now I need to figure out how to illustrate it.  I may move into Excel for that for simplicity's sake.


