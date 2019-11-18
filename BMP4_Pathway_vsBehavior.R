
BMP4_pathway_vs_behavior_forR<-read.csv("BMP4_pathway_vs_behavior_forR.csv", header=T, stringsAsFactors = F)
str(BMP4_pathway_vs_behavior_forR)
colnames(BMP4_pathway_vs_behavior_forR)

Numeric_BMP4_pathway_vs_behavior_forR<-as.matrix(BMP4_pathway_vs_behavior_forR[, c(6:10, 13:31)])
str(Numeric_BMP4_pathway_vs_behavior_forR)
write.csv(cor(Numeric_BMP4_pathway_vs_behavior_forR), "Cor_BMP4Pathway_Vs_Behavior.csv")

pdf("Adult_BMP4_VS_Phenotype.pdf")
boxplot(BMP4_pathway_vs_behavior_forR$BMP4~BMP4_pathway_vs_behavior_forR$Lineage, col=c(3, 0, 2), ylab="BMP4 Expression (log2 CPM)", xlab="Phenotype", cex.axis=1.3, cex.lab=1.3)
stripchart(BMP4_pathway_vs_behavior_forR$BMP4~BMP4_pathway_vs_behavior_forR$Lineage, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, col = 'black', cex=1.7)
mtext("F(2,15)=2.41, p=0.1237", cex=1.3)
dev.off()


anova(lm(BMP4_pathway_vs_behavior_forR$BMP4~BMP4_pathway_vs_behavior_forR$Lineage))
# Analysis of Variance Table
# 
# Response: BMP4_pathway_vs_behavior_forR$BMP4
# Df  Sum Sq Mean Sq F value Pr(>F)
# BMP4_pathway_vs_behavior_forR$Lineage  2  4.9761  2.4880  2.4103 0.1237
# Residuals                             15 15.4839  1.0323    
# 


#Note that the colors here are not the same as those typically used in our lab to represent HR/LR
plot(BMP4_pathway_vs_behavior_forR$BMP4~BMP4_pathway_vs_behavior_forR$Boli, col=as.factor(BMP4_pathway_vs_behavior_forR$Lineage), pch=18)

summary.lm(lm(BMP4_pathway_vs_behavior_forR$BMP4~BMP4_pathway_vs_behavior_forR$Boli))

# Call:
#   lm(formula = BMP4_pathway_vs_behavior_forR$BMP4 ~ BMP4_pathway_vs_behavior_forR$Boli)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.7517 -0.6917  0.2078  0.6020  2.0209 
# 
# Coefficients:
#                                       Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)                          2.2427     0.2739   8.188 4.09e-07 ***
#   BMP4_pathway_vs_behavior_forR$Boli   0.3170     0.1232   2.573   0.0204 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.951 on 16 degrees of freedom
# Multiple R-squared:  0.2927,	Adjusted R-squared:  0.2485 
# F-statistic: 6.622 on 1 and 16 DF,  p-value: 0.02041

summary.lm(lm(BMP4_pathway_vs_behavior_forR$BMP4~BMP4_pathway_vs_behavior_forR$TotalLocoScore))

# Call:
#   lm(formula = BMP4_pathway_vs_behavior_forR$BMP4 ~ BMP4_pathway_vs_behavior_forR$TotalLocoScore)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.8880 -0.5456 -0.1646  0.5146  1.9236 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                   3.0969333  0.3286348   9.424 6.23e-08 ***
#   BMP4_pathway_vs_behavior_forR$TotalLocoScore -0.0004685  0.0002356  -1.988   0.0642 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.013 on 16 degrees of freedom
# Multiple R-squared:  0.1981,	Adjusted R-squared:  0.148 
# F-statistic: 3.953 on 1 and 16 DF,  p-value: 0.06416

summary.lm(lm(BMP4_pathway_vs_behavior_forR$BMP4~BMP4_pathway_vs_behavior_forR$DistanceTraveled))
# Call:
#   lm(formula = BMP4_pathway_vs_behavior_forR$BMP4 ~ BMP4_pathway_vs_behavior_forR$DistanceTraveled)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.1545 -0.7522  0.1872  0.5753  1.6190 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                    2.642e+00  3.376e-01   7.825 7.39e-07 ***
#   BMP4_pathway_vs_behavior_forR$DistanceTraveled 1.281e-06  4.421e-05   0.029    0.977    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.131 on 16 degrees of freedom
# Multiple R-squared:  5.242e-05,	Adjusted R-squared:  -0.06244 
# F-statistic: 0.0008388 on 1 and 16 DF,  p-value: 0.9773

summary.lm(lm(BMP4_pathway_vs_behavior_forR$BMP4~BMP4_pathway_vs_behavior_forR$TimeImmobile))

# Call:
#   lm(formula = BMP4_pathway_vs_behavior_forR$BMP4 ~ BMP4_pathway_vs_behavior_forR$TimeImmobile)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.1651 -0.7453  0.2241  0.5612  1.6043 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                 2.6762785  0.4129575   6.481 7.58e-06 ***
#   BMP4_pathway_vs_behavior_forR$TimeImmobile -0.0004874  0.0053906  -0.090    0.929    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.131 on 16 degrees of freedom
# Multiple R-squared:  0.0005108,	Adjusted R-squared:  -0.06196 
# F-statistic: 0.008176 on 1 and 16 DF,  p-value: 0.9291

summary.lm(lm(BMP4_pathway_vs_behavior_forR$BMP4~BMP4_pathway_vs_behavior_forR$PercentTimeOpenArm))

# 
# Call:
#   lm(formula = BMP4_pathway_vs_behavior_forR$BMP4 ~ BMP4_pathway_vs_behavior_forR$PercentTimeOpenArm)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.9697 -0.3695  0.3648  0.5693  1.0377 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                       3.32875    0.35363   9.413 6.33e-08 ***
#   BMP4_pathway_vs_behavior_forR$PercentTimeOpenArm -0.03406    0.01361  -2.503   0.0235 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.9586 on 16 degrees of freedom
# Multiple R-squared:  0.2814,	Adjusted R-squared:  0.2365 
# F-statistic: 6.267 on 1 and 16 DF,  p-value: 0.02352


pdf("Adult_BMP4_VS_OpenArms.pdf")
plot(BMP4_pathway_vs_behavior_forR$BMP4~BMP4_pathway_vs_behavior_forR$PercentTimeOpenArm, ylab="BMP4 Expression (log2 CPM)", xlab="EPM: % Time In Open Arms", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7)
points(BMP4_pathway_vs_behavior_forR$BMP4[BMP4_pathway_vs_behavior_forR$Lineage=="HR"]~BMP4_pathway_vs_behavior_forR$PercentTimeOpenArm[BMP4_pathway_vs_behavior_forR$Lineage=="HR"],col=3, pch=20, cex=1.7)
points(BMP4_pathway_vs_behavior_forR$BMP4[BMP4_pathway_vs_behavior_forR$Lineage=="LR"]~BMP4_pathway_vs_behavior_forR$PercentTimeOpenArm[BMP4_pathway_vs_behavior_forR$Lineage=="LR"],col=2, pch=20, cex=1.7)
temp<-lm(BMP4_pathway_vs_behavior_forR$BMP4~BMP4_pathway_vs_behavior_forR$PercentTimeOpenArm)
abline(temp, lwd=2)
mtext("R-squared=0.28, Beta=-0.034, T(16)=-2.503, p=0.0235", cex=1.3)
dev.off()


pdf("Adult_BMP4_VS_FecalBoli.pdf")
plot(BMP4_pathway_vs_behavior_forR$BMP4~BMP4_pathway_vs_behavior_forR$Boli, ylab="BMP4 Expression (log2 CPM)", xlab="EPM: Fecal Boli", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7)
points(BMP4_pathway_vs_behavior_forR$BMP4[BMP4_pathway_vs_behavior_forR$Lineage=="HR"]~BMP4_pathway_vs_behavior_forR$Boli[BMP4_pathway_vs_behavior_forR$Lineage=="HR"],col=3, pch=20, cex=1.7)
points(BMP4_pathway_vs_behavior_forR$BMP4[BMP4_pathway_vs_behavior_forR$Lineage=="LR"]~BMP4_pathway_vs_behavior_forR$Boli[BMP4_pathway_vs_behavior_forR$Lineage=="LR"],col=2, pch=20, cex=1.7)
temp<-lm(BMP4_pathway_vs_behavior_forR$BMP4~BMP4_pathway_vs_behavior_forR$Boli)
abline(temp, lwd=2)
mtext("R-squared=0.29, Beta=0.317, T(16)=2.573, p=0.0204", cex=1.3)
dev.off()

   