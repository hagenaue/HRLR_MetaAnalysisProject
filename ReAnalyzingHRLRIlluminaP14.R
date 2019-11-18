library(illuminaio)
#Note - please cite
#http://f1000researchdata.s3.amazonaws.com/manuscripts/3137/f51ea0dd-a971-45ed-8495-f94faba02038_2909%20-%20kasper%20hansen.pdf

setwd("/Users/mhh/Documents/Microarray Gen/Cigdem_HRLR_Drug/Sarah Illumina P14/4194720019")

#Oops - this is just code for the practice file: 4194720019_A_Grn<-system.file("extdata", "idat", "4194720019_A_Grn.idat"

#instead:
#http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/read.idat.html
#read.idat(idatfiles, bgxfile, dateinfo = FALSE, annotation = "Symbol", tolerance = 0L, verbose = TRUE)

4194720019_A_Grn<-read.idat(c("4194720019_A_Grn.idat", "4194720019_B_Grn.idat"), "RatRef-12_V1_0_R5_11222119_A.bgx", dateinfo = FALSE, annotation = "Symbol", tolerance = 0L, verbose = TRUE)
#Says that there is "unexpected input"

#I double checked that the underscore wasn't the problem.
#4194720019_A_Grn<-read.idat(c("4194720019AGrn.idat"), "RatRef-12_V1_0_R5_11222119_A.bgx", dateinfo = FALSE, annotation = "Symbol", tolerance = 0L, verbose = TRUE)

#Oh - it is probably unhappy because of the numeric name:
SampleA<-read.idat(c("4194720019AGrn.idat"), "RatRef-12_V1_0_R5_11222119_A.bgx", dateinfo = FALSE, annotation = "Symbol", tolerance = 0L, verbose = TRUE)
#still unhappy - apparently the read.idat file comes from a different package! limma!


setwd("/Users/mhh/Documents/Microarray Gen/Cigdem_HRLR_Drug/Sarah Illumina P14")

SampleA<-readIDAT(c("4194720019AGrn.idat"))
      str(SampleA)  
      
      # $ Barcode : chr "4194720019"
      # $ Section : chr "A"
      # $ ChipType: chr "BeadChip 12x1"
      # $ Quants  :'data.frame':	23401 obs. of  10 variables:
      #   ..$ MeanBinData         : num [1:23401] 1063 114 108 113 3627 ...
      # ..$ TrimmedMeanBinData  : num [1:23401] 1073 114 110 114 3667 ...
      # ..$ DevBinData          : num [1:23401] 164.8 32.6 46.2 32.6 648.6 ...
      # ..$ MedianBinData       : num [1:23401] 1065 110 109 106 3791 ...
      # ..$ BackgroundBinData   : num [1:23401] 700 699 699 699 701 ...
      # ..$ BackgroundDevBinData: num [1:23401] 3.07 2.82 3.62 1.83 3.8 ...
      # ..$ CodesBinData        : int [1:23401] 10411 10575 20224 20356 50010 50014 50017 50020 50021 50022 ...
      # ..$ NumBeadsBinData     : int [1:23401] 51 53 49 36 47 44 40 49 45 39 ...
      # ..$ NumGoodBeadsBinData : int [1:23401] 49 50 43 36 46 42 38 45 42 36 ...
      # ..$ IllumicodeBinData   : int [1:23401] 10411 10575 20224 20356 50010 50014 50017 50020 50021 50022 ...
      # $ RunInfo : chr [1:4, 1:5] "Decoding" "Scan" "Register" "Extract" ...
      # ..- attr(*, "dimnames")=List of 2
      # .. ..$ : chr [1:4] "" "" "" ""
      # .. ..$ : chr [1:5] "Name" "SoftwareApp" "Version" "Date" ...
      
      # Note: CodeBinData=ProbeID, MeanBinData=AVG_Signal, NumGoodBeadsBinData=NBeads in GenomeStudio
      #Note: DevBinData gives the standard deviation for the bead type and can generate the BEAD_STDERR values in GenomeStudio.
      
      #Let's try the limma version too:
    
      setwd("/Users/mhh/Documents/Microarray Gen/Cigdem_HRLR_Drug/Sarah Illumina P14/4194720019")
      
      library(limma)
      
      SampleA_limma<-read.idat(c("4194720019_A_Grn.idat"), "RatRef-12_V1_0_R5_11222119_A.bgx", dateinfo = FALSE, annotation = "Symbol", tolerance = 0L, verbose = TRUE)
      str(SampleA_limma)
      
      str(SampleA_limma[[3]])
      
      table(SampleA_limma[[3]]$Status)
      
      # biotin             cy3_hyb high_stringency_hyb        housekeeping            labeling 
      # 2                   6                   1                   6                   4 
      # low_stringency_hyb            negative             regular 
      # 8                 825               22523 
      
      #Interesting - that must be the info that comes from the BGX file.
      #Seems useful.  I wonder what else the limma package does for Illumina files.
      Probe_Annotation<-SampleA_limma[[3]]
      
      str(SampleA_limma[[4]])
      #This is the AVG_Signal data
      str(SampleA_limma[[5]][1])
      #This is the Number of Beads data
      str(SampleA_limma[[5]][2])
      #This is the StDev data
      
      #Let's see if all samples have the same dimensions:
      
      SampleB_limma<-read.idat(c("4194720019_B_Grn.idat"), "RatRef-12_V1_0_R5_11222119_A.bgx", dateinfo = FALSE, annotation = "Symbol", tolerance = 0L, verbose = TRUE)
      
      str(SampleB_limma)
      #yep, looks like they have the same dimensions.
  
      SampleC_limma<-read.idat(c("4194720019_C_Grn.idat"), "RatRef-12_V1_0_R5_11222119_A.bgx", dateinfo = FALSE, annotation = "Symbol", tolerance = 0L, verbose = TRUE)    
      
      SampleD_limma<-read.idat(c("4194720019_D_Grn.idat"), "RatRef-12_V1_0_R5_11222119_A.bgx", dateinfo = FALSE, annotation = "Symbol", tolerance = 0L, verbose = TRUE)  
  
      SampleE_limma<-read.idat(c("4194720019_E_Grn.idat"), "RatRef-12_V1_0_R5_11222119_A.bgx", dateinfo = FALSE, annotation = "Symbol", tolerance = 0L, verbose = TRUE) 
      
      SampleF_limma<-read.idat(c("4194720019_F_Grn.idat"), "RatRef-12_V1_0_R5_11222119_A.bgx", dateinfo = FALSE, annotation = "Symbol", tolerance = 0L, verbose = TRUE) 
  
      SampleG_limma<-read.idat(c("4194720019_G_Grn.idat"), "RatRef-12_V1_0_R5_11222119_A.bgx", dateinfo = FALSE, annotation = "Symbol", tolerance = 0L, verbose = TRUE)   
      
      SampleH_limma<-read.idat(c("4194720019_H_Grn.idat"), "RatRef-12_V1_0_R5_11222119_A.bgx", dateinfo = FALSE, annotation = "Symbol", tolerance = 0L, verbose = TRUE)   
    
      SampleI_limma<-read.idat(c("4194720019_I_Grn.idat"), "RatRef-12_V1_0_R5_11222119_A.bgx", dateinfo = FALSE, annotation = "Symbol", tolerance = 0L, verbose = TRUE)  
    
      SampleJ_limma<-read.idat(c("4194720019_J_Grn.idat"), "RatRef-12_V1_0_R5_11222119_A.bgx", dateinfo = FALSE, annotation = "Symbol", tolerance = 0L, verbose = TRUE)  
   
      SampleK_limma<-read.idat(c("4194720019_K_Grn.idat"), "RatRef-12_V1_0_R5_11222119_A.bgx", dateinfo = FALSE, annotation = "Symbol", tolerance = 0L, verbose = TRUE)     
      
      SampleL_limma<-read.idat(c("4194720019_L_Grn.idat"), "RatRef-12_V1_0_R5_11222119_A.bgx", dateinfo = FALSE, annotation = "Symbol", tolerance = 0L, verbose = TRUE)     
      
      
      AVG_Signal<-cbind(SampleA_limma[[4]], SampleB_limma[[4]], SampleC_limma[[4]], SampleD_limma[[4]],SampleE_limma[[4]],SampleF_limma[[4]],SampleG_limma[[4]],SampleH_limma[[4]],SampleI_limma[[4]], SampleJ_limma[[4]], SampleK_limma[[4]], SampleL_limma[[4]] )
     
      write.csv(data.frame(Probe_Annotation, AVG_Signal), "Illumina_HRLR_P14_F15_AVG_Signal.csv")
        
      str(AVG_Signal)
      log2AVG_Signal<-log2(AVG_Signal)
      
      #Note - this doesn't have any background subtraction - If I remember correctly, the background subtraction didn't help much when looking at the human data (Pritzker960).
      
      
      #Sample IDs from "Sara HR_LR Sample Sheet.xls"
      SampleID<-c("LR_26", "HR_27", "LR_28", "HR_29", "LR_30", "HR_33", "LR_38", "HR_39", "HR_41", "LR_42", "LR_46", "HR_47")
      
      Phenotype<-c("LR", "HR", "LR", "HR", "LR", "HR", "LR", "HR", "HR", "LR", "LR", "HR")
      
      colnames(log2AVG_Signal)<-SampleID
      
      getwd()
      setwd("/Users/mhh/Documents/Microarray Gen/Cigdem_HRLR_Drug/Sarah Illumina P14")
      
      #Time to run some QC:
      
      #Outputs a boxplot illustrating the distribution of log-transformed signal for each subject:
      png("Boxplot logAVGSignal.png")
      boxplot(data.frame(log2AVG_Signal), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of log signal values per sample (1 box=all probes)", xlab="Sample ID", ylab="Log Signal")
      dev.off()
      
      #3 of the HRs have less signal than everyone else.  Could be meaningful, but hard to say.
      
      heatmap(cor(log2AVG_Signal))
      #partially segregates by phenotype but not perfectly.
      
      library(preprocessCore)
      #Quantile normalize the log-transformed signal data (note: this is called "filtered" because I'm recycling old code):
      NormFiltered<-normalize.quantiles(log2AVG_Signal)
      row.names(NormFiltered)<-row.names(log2AVG_Signal)
      colnames(NormFiltered)<-colnames(log2AVG_Signal)
      
      #Output the quantile normalized data:
      write.table(NormFiltered, "Quantile Normalized filtered data.txt", sep="\t")
      
      write.csv(data.frame(Probe_Annotation, NormFiltered), "Illumina_HRLR_P14_F15_NormData.csv")
      
      png("Heatmap_samplesampleCorMatrix_NormalizedData.png")
      heatmap(cor(NormFiltered))
      dev.off()
      #Still doesn't segregate by phenotype.
      #check for outliers?
      
      write.csv(cor(NormFiltered), "SampleSampleCorMatrix_NormalizedData.csv" )
      
      
      png("09 Boxplot Sample Sample Correlations.png")
      boxplot(data.frame(cor(NormFiltered)), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of sample-sample correlations", xlab="Subject", ylab="Sample-Sample Correlations")
      Median10thQuantile<-median(apply((cor(NormFiltered)), 1, quantile, 0.1))
      MedianQuantile<-median(apply((cor(NormFiltered)), 1, quantile, 0.5))
      abline(a=Median10thQuantile, b=0, col=2)
      abline(a=MedianQuantile, b=0, col=3)
      mtext(paste("Median Sample-Sample Correlation=", round(MedianQuantile, digits=3), sep=" ")) 
      dev.off()
      
      
      #HR_27 might qualify as an outlier, although the sample-sample correlations are still high >0.98, everyone else is >0.99 , if the data were filtered it might be more extreme.
      
    #I'm going to leave him in for right now.
      
      pcaNormFiltered<-prcomp(t(NormFiltered))
      tmp<-pcaNormFiltered$x[,1:4]
      write.table(tmp, "PCA_1_4.txt", sep="\t")
      
      PC1<-pcaNormFiltered$x[,1]
      PC2<-pcaNormFiltered$x[,2]
      
      PC3<-pcaNormFiltered$x[,3]
      PC4<-pcaNormFiltered$x[,4]
      
      #Output a scree plot for the PCA:
      png("09 PCA Scree Plot1.png")
      plot(summary(pcaNormFiltered)$importance[2,]~(c(1:12)), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
      dev.off()
      
      png("09 PCA Scree Plot2.png")
      plot(summary(pcaNormFiltered)$importance[3,]~(c(1:12)), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Cumulative Proportion of Variance Explained", col=3)
      dev.off()
      
      #Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
      png("PC1 vs PC2.png")
      plot(PC1~PC2, main="Principal Components Analysis of Normalized Data", col=as.factor(Phenotype))
      #legend(min(PC1), max(PC2)+10, c("Control", "MDD"), text.col=c(3, 2), pch=19, col=c(3, 2))
      dev.off()
      
      
      #Output a scatterplot illustrating the relationship between Principal components 2 & 3 (PC3 & PC4):
      png("PC2 vs PC3.png")
      plot(PC2~PC3, main="Principal Components Analysis of Normalized Data", col=as.factor(Phenotype))
      #legend(min(PC3), max(PC4)+10, c("Control", "MDD"), text.col=c(3, 2), pch=19, col=c(3, 2))
      dev.off()
      
      #Output a scatterplot illustrating the relationship between Principal components 3 & 4 (PC3 & PC4):
      png("PC3 vs PC4.png")
      plot(PC3~PC4, main="Principal Components Analysis of Normalized Data", col=as.factor(Phenotype))
      #legend(min(PC3), max(PC4)+10, c("Control", "MDD"), text.col=c(3, 2), pch=19, col=c(3, 2))
      dev.off()
      
      #The data segregates in PC2/PC3
      
      png("Boxplot_PC1vsPhenotype.png")
      boxplot(PC1~Phenotype, col=2)
      dev.off()
      
      png("Boxplot_PC2vsPhenotype.png")
      boxplot(PC2~Phenotype, col=2)
      dev.off() 
      
      png("Boxplot_PC3vsPhenotype.png")
      boxplot(PC3~Phenotype, col=2)
      dev.off()
      
      png("Boxplot_PC4vsPhenotype.png")
      boxplot(PC4~Phenotype, col=2)
      dev.off()
      
      colnames(NormFiltered)
      
      png("Heatmap_samplesampleCorMatrix_NormalizedData_NoHR27.png")
      heatmap(cor(NormFiltered[,-2]))
      dev.off()
      #Still doesn't segregate by phenotype.
      #check for outliers?
      
    
      
      pcaNormFiltered<-prcomp(t(NormFiltered[,-2]))
      tmp<-pcaNormFiltered$x[,1:4]
      # write.table(tmp, "PCA_1_4.txt", sep="\t")
      # 
       PC1<-pcaNormFiltered$x[,1]
      PC2<-pcaNormFiltered$x[,2]
      # 
      PC3<-pcaNormFiltered$x[,3]
      PC4<-pcaNormFiltered$x[,4]
      # 
      #Output a scree plot for the PCA:
      # png("09 PCA Scree Plot1.png")
      # plot(summary(pcaNormFiltered)$importance[2,]~(c(1:12)), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
      # dev.off()
      # 
      # png("09 PCA Scree Plot2.png")
      # plot(summary(pcaNormFiltered)$importance[3,]~(c(1:12)), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Cumulative Proportion of Variance Explained", col=3)
      # dev.off()
      
      #Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
      png("PC1 vs PC2_noHR27.png")
      plot(PC1~PC2, main="Principal Components Analysis of Normalized Data", col=as.factor(Phenotype[-2]))
      #legend(min(PC1), max(PC2)+10, c("Control", "MDD"), text.col=c(3, 2), pch=19, col=c(3, 2))
      dev.off()
      
      
      #Output a scatterplot illustrating the relationship between Principal components 2 & 3 (PC3 & PC4):
      png("PC2 vs PC3_noHR27.png")
      plot(PC2~PC3, main="Principal Components Analysis of Normalized Data", col=as.factor(Phenotype[-2]))
      #legend(min(PC3), max(PC4)+10, c("Control", "MDD"), text.col=c(3, 2), pch=19, col=c(3, 2))
      dev.off()
      
      #Output a scatterplot illustrating the relationship between Principal components 3 & 4 (PC3 & PC4):
      png("PC3 vs PC4_noHR27.png")
      plot(PC3~PC4, main="Principal Components Analysis of Normalized Data", col=as.factor(Phenotype[-2]))
      #legend(min(PC3), max(PC4)+10, c("Control", "MDD"), text.col=c(3, 2), pch=19, col=c(3, 2))
      dev.off()
      
  