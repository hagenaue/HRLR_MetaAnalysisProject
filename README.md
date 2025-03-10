# HRLR_MetaAnalysisProject

This repository is part of the code release to accompany the manuscript "Genetic liability for internalizing versus externalizing behavior manifests in the developing and adult hippocampus: Insight from a meta-analysis of transcriptional profiling studies in a selectively-bred rat model" by Birt, Hagenauer, et al. 
https://www.sciencedirect.com/science/article/abs/pii/S0006322320316309 (https://doi.org/10.1016/j.biopsych.2020.05.024)

## Please note: The main body of code is available in Isabelle Birt's repository:
 https://github.com/isabellie4/PhenotypeProject
The code included in that repository focuses on the analysis and meta-analysis of the transcriptional profiling data (microarray, RNA-Seq). 

## General description of the current repository: 
https://github.com/hagenaue/HRLR_MetaAnalysisProject This repository was created by Megan Hagenauer, but also includes some code written by other contributors (Dr. Fan Meng, Isabelle Birt). The code included in this repository focuses on the basic preprocessing of the some of the transcriptional profiling data (microarray, RNA-Seq), as well as downstream analyses that were performed using the results from the meta-analysis, including the evaluation of overlap with genetic findings.

## Description of the individual files in the current respository:
Code used to analyze and produce figures for Dr. Kathryn Hilde and Alex Stefanov Bmp4 qPCR studies:
* PCR_ReanalysisCode_P14_P90.R
* Data: Figshare: https://doi.org/10.6084/m9.figshare.10321658.v1

Code used to perform some of the initial preprocessing of Dr. Sarah Clinton’s developmental Affymetrix microarray study from the F6 generation: MBNI_AffymetrixRgU34A_F6
* HRLR_ClintonOldDevelopment_Affy_RMA_ReAnnotation.R
* Dataset: Gene Expression Omnibus ID GSE29552: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29552

Code used to perform some of the initial preprocessing of Dr. Sarah Clinton’s P14 Illumina microarray study: MBNI_IlluminaRatRef12v1_F15
* ReAnalyzingHRLRIlluminaP14.R
* Dataset: Gene Expression Omnibus ID GSE140594: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140594

Code used to perform some of the initial preprocessing of Dr. Sarah Clinton’s P14 Affymetrix microarray study: MBNI_AﬀymetrixRae230_ F15
* HRLR_P14_Affy230_RMA_ReAnnotation.R
* Dataset: Gene Expression Omnibus ID GSE140595: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140595

Code used to perform the initial preprocessing of Dr. Sarah Clinton’s small P14 and adult RNA-Seq study: MBNI_RNASeq_F29 dataset (written by Dr. Fan Meng):
* Old RNASeq code.docx
* AnalysisCodeFromFan.txt
* Dataset: Gene Expression Omnibus ID GSE140597: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140597

Code used to analyze and plot the relationship between the behavioral data and pre-processed log(2) FPM count data from Dr. Cigdem Aydin and Dr. Peter Blandino’s adult RNA-Seq studies (MBNI_RNASeq_F37 and MBNI_RNASeq_43, respectively) 
* Blandino_BehaviorAnalysis_MHedit.R
* BMP4_Pathway_vsBehavior.R
* Dataset: Gene Expression Omnibus ID GSE140598: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140598
* Dataset: Gene Expression Omnibus ID GSE140287: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140287

Isabelle Birt’s developmental meta-analysis code (see Repository: https://github.com/isabellie4/PhenotypeProject), tweaked to change the formatting of figures for the manuscript:
* MetaAnalysis_Development.R

Isabelle Birt’s adult meta-analysis code (see Repository: https://github.com/isabellie4/PhenotypeProject), tweaked for four purposes:
1) To evaluate the correlation between the results from individual datasets
2) To compare the findings from the full meta-analysis to findings produced by a meta-analysis of just the later generation (F37 & F43) RNA-Seq studies
3) To perform a post-hoc evaluation of the importance of model type (fixed effects vs. random effects, including generation as a co-variate)
4) To change the formatting of figures for the manuscript
* Adult Meta-Analysis_MH_2018_04_20.R

Code used to compare the transcriptional profiling meta-analysis results to the exome sequencing results (quantitative trait loci coordinates from Zhou et al. 2019 Proc Natl Acad Sci USA. 116: 13107–13115, https://www.pnas.org/doi/pdf/10.1073/pnas.1820410116):
* MappingMetaAnalysis_toChromosomeStart.R

Code used to compare the transcriptional profiling meta-analysis results to the exome sequencing results (HR/LR segregating variants from Zhou et al. 2019 Proc Natl Acad Sci USA. 116: 13107–13115, https://www.pnas.org/doi/pdf/10.1073/pnas.1820410116):
* EnrichmentOfSegregatingVariants_20191011.R

Code used to extract the top meta-analysis results in the genomic regions highlighted by positional gene expression (PGE) analysis:
* JoiningPGEresultsWOurMetaanalysisResults.R

Isabelle Birt’s fgsea code (see Repository: https://github.com/isabellie4/PhenotypeProject), tweaked for three purposes:
1)  To compare the findings from the full meta-analysis to findings produced by a meta-analysis of just the later generation (F37 & F43) RNA-Seq studies. (note: Within this code, I also outputted the EntrezIDs necessary for running a PGE analysis of the late generation meta-analysis output too).
2) To perform a post-hoc evaluation of the importance of different fGSEA parameters (.gmt structure, maxSize, # of permutations).
3) To change the formatting of figures for the manuscript
* fgsea_MH_2018_04_12.R

Code used to explore the hippocampal specific co-expression networks (e.g. for enrichment within particular functional categories/cell types) – not really used in the manuscript:
 * CoexpressionNetworkMouseHumanHC.R

Code used to explore the overlap between hippocampal specific co-expression networks and top meta-analysis genes and the Mousebrain.org single-cell hippocampal expression database:
* Code_ZeiselvsHippocampalExpression_HC.R

Code used to examine the overlap between the top differentially expressed genes identified in the hippocampus of bHR/bLR rats versus the top differentially-expressed genes identified in the hippocampus of other bred rat models of internalizing/externalizing behavior. This code begins with a re-analysis of a publicly-available microarray dataset comparing Flinders Sensitive Line vs. Flinders Resistant Line (Blaveri et al. 2010), and ends with a broader analysis comparing across multiple models, using a database of differentially-expressed genes that I compiled from the re-analysis of Blaveri et al. (2010) and the published results from 8 other studies:
* Flinders_Code_HC_20190124.R
* Dataset: Gene Expression Omnibus ID GSE20388: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20388
