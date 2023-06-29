library("ggplot2")
library("dplyr")
library("reshape")
library("tidyverse")

#####################################
######PCA PLOTS########
setwd("~/Desktop/Pan_Disease/PAPER_VERSION_4/")

data <- read.csv("~/Desktop/Pan_Disease/AllSamples_Cleaned_OverlappingWithDb_PanDisease_3544_PCs_Sex_Cohort_DiseaseStatus_Origin_SubjectID_AgeAtDiagnosis_AllSubPhenotypes_WithAllAAb", sep="\t", header=T)
data$Group <- ifelse(data$DiseaseStatus == 2, "SIAD Case", "Control")
data$Group <- factor(data$Group, levels=c("SIAD Case","Control"))

pdf("Figure_PC1_PC2_PanDisease_BigFont.pdf")
ggplot(data, aes(x=PC1,y=PC2,color=Group)) + labs(col="Disease status") + labs(x = "PC1") + labs(y = "PC2") + scale_color_manual(values = c("darkorange3","darkblue")) + geom_point() +
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),
        legend.title=element_text(size=15), legend.text=element_text(size=15)              
  )
dev.off()

pdf("Figure_PC2_PC3_PanDisease_BigFont.pdf")
ggplot(data, aes(x=PC2,y=PC3,color=Group)) + labs(col="Disease status") + labs(x = "PC2") + labs(y = "PC3") + scale_color_manual(values = c("darkorange3","darkblue")) + geom_point() +
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),
        legend.title=element_text(size=15), legend.text=element_text(size=15)              
  )
dev.off()

pdf("Figure_PC3_PC4_PanDisease_BigFont.pdf")
ggplot(data, aes(x=PC3,y=PC4,color=Group)) + labs(col="Disease status") + labs(x = "PC3") + labs(y = "PC4") + scale_color_manual(values = c("darkorange3","darkblue")) + geom_point() +
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),
        legend.title=element_text(size=15), legend.text=element_text(size=15)              
  )
dev.off()


pdf("Figure_PC4_PC5_PanDisease_BigFont.pdf")
ggplot(data, aes(x=PC4,y=PC5,color=Group)) + labs(col="Disease status") + labs(x = "PC4") + labs(y = "PC5") + scale_color_manual(values = c("darkorange3","darkblue")) + geom_point() +
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),
        legend.title=element_text(size=15), legend.text=element_text(size=15)              
  )
dev.off()
#####################################



#####################################
####BAR PLOT Subphenotype#####
setwd("~/Desktop/Pan_Disease/PAPER_VERSION_4/")
data <- read.csv("File_SIADsubphenoStrat_4BarPlot", head=T, sep="\t")
head(data)
data$Proportion <- data$N / data$Tot

data$Phenotype <- factor(data$Phenotype, levels=c("ANA","dsDNA","SSA","SSA-Ro52","SSA-Ro60","SSB","RNP","SM","RF",
                                                  "ARTHRITIS","CANCER","CARDIAC INVOLVEMENT","LUNG INVOLVEMENT","RAYNAUD","RENAL DISORDER","SKIN INVOLVEMENT"))
data$SIAD <- factor(data$SIAD, levels=c("SLE","pSS","Myositis"))


####Absolute count#####
pdf("SIAD_Subphenotypes_PatientNumber.pdf")
ggplot(data=data, aes(x=Phenotype, y=N, fill=SIAD)) +
  geom_bar(stat="identity") +
  ylab("# patients") +
  xlab("") +
  theme(axis.text.x = element_text(size=10, angle=75, hjust=1), axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=10), axis.title.y = element_text(size=15)) +
  scale_fill_manual(values = c("#736F6E", "#98AFC7","#153E7E"))  
dev.off()

####Proportion#####
pdf("SIAD_Subphenotypes_PatientFraction.pdf")
ggplot(data=data, aes(x=Phenotype, y=Proportion, fill=SIAD)) +
  geom_bar(stat="identity") +
  ylab("Fraction of patients") +
  xlab("") +
  theme(axis.text.x = element_text(size=10, angle=75, hjust=1), axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=10), axis.title.y = element_text(size=15)) +
  scale_fill_manual(values = c("#736F6E", "#98AFC7","#153E7E"))
  dev.off()
  #####################################

  


#####################################
#################################CADD COMPARISON######
#setwd("~/Desktop/Pan_Disease/PAPER_FIRST_VERSION/")
#cadd_comp <- read.csv("Vars_CrossDisease_cadd13_and_cadd16_AllChromosomes_250KRandomVars", head=F, sep="\t")
#cor.test(cadd_comp$V5, cadd_comp$V7, method="spearman") #####Spearman rho0.81     #####Pearson 0.90
#cadd_comp_100K <- read.csv("Vars_CrossDisease_cadd13_and_cadd16_AllChromosomes_100KRandomVars", head=F, sep="\t")
#cor.test(cadd_comp_100K$V5, cadd_comp_100K$V7, method="spearman") #####Spearman rho0.81
#pdf("Plot_Cadd_DifferentVersion_100K_Spearman.pdf")
#plot(cadd_comp_100K$V5, cadd_comp_100K$V7, col="slateblue", xlab="CADD v1.3 (Raw)", ylab ="CADD v1.6 (Raw)", pch=16)
#abline(lm(cadd_comp_100K$V7 ~ cadd_comp_100K$V5), lty=6)
#legend(x = "bottomright", legend = "Spearman rho = 0.81")
#dev.off()


######VENN SIMILAR CLINICAL GROUPS#######
#library(venn)
####Create tables with genes per group and clinical subgroup#####
#setwd("~/Desktop/Pan_Disease/PAPER_VERSION_4/Venn_SimilarClinicalSubgroups_OverlappingGenes_GRSbased/")

#Arthritis <- read.csv("Arthritis.txt", sep="\t", header = F)
#Cancer <- read.csv("Cancer.txt", sep="\t", header = F)
#RenalDisorder <- read.csv("RenalDisorder.txt", sep="\t", header = F)
#SkinInvolvement <- read.csv("SkinInvolvement.txt", sep="\t", header = F)

#tot <- list(
#  "Arthritis" = Arthritis$V1,
#  "Cancer" = Cancer$V1,
#  "Renal\nDisorder" = RenalDisorder$V1,
#  "Skin\nInvolvement" = SkinInvolvement$V1
#)

#pdf("Venn_ClinicalSubgroups_OverlappingGenes_GRSbased.pdf")
#venn(tot, ilab=TRUE, zcolor = "style", sncs=0.8, ilcs=1)
#dev.off()


#######POWER CALCULATION#######
#library(genpwr)

#OR_1_2 <- power.calc(N=3544, Case.Rate=c(0.6468), k=NULL,
#                     MAF=seq(0.01, 0.30, 0.01),Alpha=c(0.000000125),
#                     True.Model='Additive', Test.Model='Additive', OR=1.2)

#OR_1_3 <- power.calc(N=3544, Case.Rate=c(0.6468), k=NULL,
#                     MAF=seq(0.01, 0.30, 0.01),Alpha=c(0.000000125),
#                     True.Model='Additive', Test.Model='Additive', OR=1.3)

#OR_1_4 <- power.calc(N=3544, Case.Rate=c(0.6468), k=NULL,
#                     MAF=seq(0.01, 0.30, 0.01),Alpha=c(0.000000125),
#                     True.Model='Additive', Test.Model='Additive', OR=1.4)

#OR_1_5 <- power.calc(N=3544, Case.Rate=c(0.6468), k=NULL,
#                 MAF=seq(0.01, 0.30, 0.01),Alpha=c(0.000000125),
#                 True.Model='Additive', Test.Model='Additive', OR=1.5)

#OR_2_0 <- power.calc(N=3544, Case.Rate=c(0.6468), k=NULL,
#                     MAF=seq(0.01, 0.30, 0.01),Alpha=c(0.000000125),
#                     True.Model='Additive', Test.Model='Additive', OR=2)

#OR_2_5 <- power.calc(N=3544, Case.Rate=c(0.6468), k=NULL,
#                     MAF=seq(0.01, 0.30, 0.01),Alpha=c(0.000000125),
#                     True.Model='Additive', Test.Model='Additive', OR=2.5)

#OR_3_0 <- power.calc(N=3544, Case.Rate=c(0.6468), k=NULL,
#                     MAF=seq(0.01, 0.30, 0.01),Alpha=c(0.000000125),
#                     True.Model='Additive', Test.Model='Additive', OR=3)

#OR_3_5 <- power.calc(N=3544, Case.Rate=c(0.6468), k=NULL,
#                     MAF=seq(0.01, 0.30, 0.01),Alpha=c(0.000000125),
#                     True.Model='Additive', Test.Model='Additive', OR=3.5)

#OR_4_0 <- power.calc(N=3544, Case.Rate=c(0.6468), k=NULL,
#                     MAF=seq(0.01, 0.30, 0.01),Alpha=c(0.000000125),
#                     True.Model='Additive', Test.Model='Additive', OR=4)

#total <- rbind(OR_1_2, OR_1_3, OR_1_4, OR_1_5, OR_2_0, OR_2_5, OR_3_0, OR_3_5, OR_4_0)
#colnames(total) <- c("TestModel","TrueModel","MAF","OR","NTotal","NCases","NControls","CaseRate","Power")

#pdf("Desktop/Pan_Disease/PAPER_VERSION_4//PowerPlot.pdf")
#ggplot(total, aes(x=MAF, y=Power, col=as.factor(OR))) + 
#  geom_line() +
#  scale_x_continuous(breaks=c(0.01,0.05,0.10,0.15,0.20,0.25,0.30)) +
#  scale_y_continuous(breaks=seq(0,1,0.1)) +
#  ylab("Power") +
#  xlab("MAF") +
#  labs(colour = "OR")
#dev.off()
#####################################
  

  
#####################################
####################################################
#Bonferroni threshold is 0.05/106484 
#Check variants (GC corrected) exceeding Bonferroni
#awk '$1 == chrom && $4 < 0.00000047 {print}' AllSamples_Cleaned_OverlappingWithDb_PanDisease_3544_NoFalsePositives_maf001_LogReg4PCASex.assoc.logistic.adjusted
setwd("Desktop/Pan_Disease/PAPER_VERSION_4/SingleVariant_Association/")
a <- read.csv("PanDisease_SingleVarAssoc_PValsGC_4Plot.txt", sep="\t", head=T)
new.names <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X")
pdf("SIAD_Cases_vs_Controls_ManhattanPlot.pdf")
qq(a$P)
manhattan(a, genomewideline = -log10(0.05/nrow(a)), suggestiveline = F, cex = 0.5, cex.axis = 0.8, col = c("darkblue", "darkorange3"), chrlabs = new.names, ylim=c(0,50))
dev.off()
#####################################



##########
#setwd("~/Desktop/Pan_Disease/PAPER_VERSION_4/Venn_ClinicalSubpheno_OverlappingPatients/")
#data <- read.csv("Pairwise_ClinicalOverlap_4Heatmap", sep="\t", head=T)
#head(data)
#mat <- matrix(NA, 7, 7)
#mat[upper.tri(mat)] <- data$Overlap
#colnames(mat) <- c("Arthritis","Cancer","Cardiac\nInvolvement","Lung\nInvolvement","Raynaud","Renal\nDisorder","Skin\nInvolvement")
#rownames(mat) <- c("Arthritis","Cancer","Cardiac\nInvolvement","Lung\nInvolvement","Raynaud","Renal\nDisorder","Skin\nInvolvement")

#pdf("HeatMap_ClinicalSubphenotypes_PairwiseFractionOverlap.pdf")
#melt(mat) %>%
#  ggplot(aes(Var1, Var2, fill = value)) + geom_tile() +
#  geom_text(aes(label = value)) +
#  scale_fill_gradient2(mid ="tomato", limit = c(0,1),name="Fraction overlap", na.value = "white") +
#  theme_classic() +
#  theme(axis.title = element_blank(),
#        axis.text.x = element_text(size=8))
#dev.off()




#####################################

##########
####BAR PLOT gene significance#####
setwd("~/Desktop/Pan_Disease/PAPER_VERSION_4/SUMMARY_FILES_PVals_Aggregate/")
data <- read.csv("File_SignificantAndNonSignificantGeneNumber_4BarPlot_v4", head=T, sep="\t")
head(data)
data$Proportion <- data$CorrectNumber / data$Tot

data$Significance <- factor(data$Significance, levels=c("Bonferroni","FDR","Non significant"))
data$Group <- factor(data$Group, levels=c("GS","SKAT-O"))


####Absolute count#####
pdf("SignificanceBreakdown_GeneNumber_v4.pdf")
ggplot(data=data, aes(x=Group, y=CorrectNumber, fill=Significance)) +
  geom_bar(stat="identity") +
  ylab("# genes") +
  xlab("") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=15), axis.title.y = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15)) +
  scale_fill_manual(values = c("darkred", "#fdae61","#153E7E"))  
dev.off()

####Proportion#####
#pdf("SignificanceBreakdown_GeneFraction.pdf")
#ggplot(data=data, aes(x=Group, y=Proportion, fill=Significance)) +
#  geom_bar(stat="identity") +
#  ylab("Fraction of genes") +
#  xlab("") +
#  theme(axis.text.x = element_text(size=20),
#        axis.text.y = element_text(size=15), axis.title.y = element_text(20),
#        legend.title = element_text(size=15),
#        legend.text = element_text(size=15)) +
#  scale_fill_manual(values = c("darkred", "#fdae61","#153E7E"))
#dev.off()

########################################################
########CONSIDERING SKATO CADD WEIGHTED#################
####COMPARISON BASED ON GRS SIGNIFICANT#################
########################################################

setwd("~/Desktop/Pan_Disease/PAPER_VERSION_4/SUMMARY_FILES_PVals_Aggregate/")
GS <- read.csv("GeneRiskScore_AllVars_AllGenes/CADD13_RawScore_LogitModels_PanDisease_CASES_vs_CONTROLS_AllVariants_Bonferroni_FDR.txt_sorted",head=T, sep="\t")
head(GS)
skato <- read.csv("SKATO_CADDWeighted_CADDWithLinearConversion_AllGenes/Results_SKATO_COV_CaddWeighted_CaddWithLinearConversion_PanDisease_CASES_VS_CONTROLS_Complete_sorted.txt",head=T, sep="\t")
head(skato)
#burden <- read.csv("BURDEN_CADDWeighted_CADDWithLinearConversion_AllGenes/Results_BURDEN_COV_CaddWeighted_CaddWithLinearConversion_PanDisease_CASES_VS_CONTROLS_Complete_sorted.txt",head=T, sep="\t")

colnames(skato)[1] <- "GENE"
skato <- skato %>% select("GENE","P.value","Bonferroni","FDR")

tot <- merge(GS, skato, by="GENE", all=T)

head(tot)

colnames(tot) <- c("GENE","Praw_GS","Pbonferroni_GS","GS","Praw_SKATO","Pbonferroni_SKATO","SKAT-O")

head(tot)
cor.test(-log10(tot$GS), -log10(tot$"SKAT-O"), method="spearman") ####Spearman rho=0.53; P<2.2e-16
cor.test(-log10(tot$GS), -log10(tot$"SKAT-O"), method="pearson") ####pearson rho=0.76; P<2.2e-16

plot(-log10(tot$GS), -log10(tot$"SKAT-O"))

mygene <- subset(tot, tot$GS <0.05)
nrow(mygene)

plot(-log10(mygene$GS), -log10(mygene$"SKAT-O"))
cor.test(-log10(mygene$GS), -log10(mygene$"SKAT-O"), method="spearman")  ####Spearman rho=0.69; P=9.059e-08
cor.test(-log10(mygene$GS), -log10(mygene$"SKAT-O"), method="pearson")  ####pearson rho=0.74; P=3.251e-09


plot(mygene$GS, mygene$"SKAT-O")
cor.test(mygene$GS, mygene$"SKAT-O", method="spearman")


#####heatmap######
library(scales)
mygene_ok <- mygene %>% select("GENE","GS","SKAT-O")
mygene_ok
mygene_4heatmap <- melt(mygene_ok)
head(mygene_4heatmap)
mygene_4heatmap$log10P <- -log10(mygene_4heatmap$value)

#mygene_4heatmap$variable <- factor(mygene_4heatmap$variable, levels=c("GRS","GenePy","SKAT-O"))

#library(RColorBrewer)
#col <- rev(brewer.pal(7, "YlGnBu"))
#col
pdf("Heatmap_SignificantGenesInGS_CompWithSKATO_MHCBlueColor_Italics_GOOD.pdf")

ggplot(mygene_4heatmap, aes(variable, reorder(GENE,log10P), fill=log10P)) + 
  geom_tile() +
  #scale_fill_gradientn(colors = c("#abdda4","#e6f598",
  #                                "#fee08b","#fdae61","#f46d43",
  #                                "#d53e4f")) +
  #scale_fill_gradientn(colors = c( "khaki1","#fee08b","#fdae61",
  #                                 "#f46d43","#d53e4f","darkred")) +
  #scale_fill_gradientn(colors=col) +             
  #scale_fill_viridis_c() +
  scale_fill_gradient(low = "khaki1", high = "darkred", limits=c(0, 55), oob=squish)+
  labs(y="", x="") +
  labs(fill='-log10P(FDR)') +
  theme_classic() +
  scale_x_discrete(expand = c(0, 0)) +
  theme(axis.text.y = element_text(colour = c("black", "black","black","black","black",
                                              "black","black","black","black","black",
                                              "black","black","black","black","black",
                                              "black","black","black","black","dodgerblue3",
                                              "dodgerblue3","black","black","dodgerblue3","black",
                                              "dodgerblue3","dodgerblue3","dodgerblue3","black","dodgerblue3",
                                              "dodgerblue3","dodgerblue3","dodgerblue3","dodgerblue3","dodgerblue3",
                                              "dodgerblue3","dodgerblue3","dodgerblue3","dodgerblue3","dodgerblue3",
                                              "dodgerblue3","dodgerblue3","dodgerblue3","dodgerblue3","dodgerblue3",
                                              "dodgerblue3"),face = "italic")) +
  theme(axis.text.x = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15)) 

dev.off()

########################################################
########################################################
########################################################


#####Venn FDR significant between GRS and SKATO#####
#head(mygene)
#SKATO_4Venn <- subset(mygene, mygene$"SKAT-O" < 0.05)[,1]
#GS_4Venn <- subset(mygene, mygene$GS < 0.05)[,1]
#ForVenn <- list(
#  "SKAT-O" = SKATO_4Venn,
#  "GS" = GS_4Venn
#)
#pdf("Venn_FDRSignificantGenes_46SignificantInGS.pdf")
#venn(ForVenn, ilab=TRUE, zcolor = "style", sncs=2, ilcs=2)
#dev.off()

########################################################
########################################################
########################################################




#########VOI frequency and phyloP score#######
#####DUSP1####
setwd("~/Desktop/Pan_Disease/PAPER_VERSION_4/DUSP1_v4/")
DUSP1 <- read.csv("DUSP1_4Plot", sep="\t", head=T)
DUSP1 <- read.csv("DUSP1_4Plot_WithVarIDs", sep="\t", head=T)

head(DUSP1)

DUSP1$DeltaPosMinusNeg <- DUSP1$SkinInvolvementPos - DUSP1$SkinInvolvementNeg
DUSP1$Group <- ifelse(DUSP1$SkinInvolvementPos >= 0.01 & DUSP1$SkinInvolvementNeg >= 0.01, "Common", "Rare")
DUSP1$Effect <- ifelse(DUSP1$DeltaPosMinusNeg < 0, "Protective", "Risk")
DUSP1$Effect[DUSP1$SkinInvolvementPos == 0 & DUSP1$SkinInvolvementNeg == 0] <- "Neutral"

#pdf("DUSP1_BarPot.pdf", width = 15, height = 7)
#ggplot(DUSP1,aes(Coord,DeltaPosMinusNeg, fill=Effect, color=Effect)) +
#  geom_bar(stat="identity") +
#  scale_fill_manual(values=c("gray41", "royalblue2", "darkred")) +
#  scale_color_manual(values=c("gray41", "royalblue2", "darkred")) +
#  theme(text = element_text(size=10),
#        axis.text.x = element_text(angle=75, hjust=1, size=11),
#        axis.text.y = element_text(size=13), axis.title.y = element_text(size=13),
#        legend.title = element_text(size=13),
#        legend.text = element_text(size=13)) +
#  xlab("") +
#  ylab("Delta allele frequency")
#dev.off()


pdf("DUSP1_BarPot_WithVarIDs.pdf", width = 15, height = 7)
ggplot(DUSP1,aes(factor(Coord),DeltaPosMinusNeg, fill=Effect, color=Effect)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("gray41", "royalblue2", "darkred")) +
  scale_color_manual(values=c("gray41", "royalblue2", "darkred")) +
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5, size=13),
        axis.text.y = element_text(size=13), axis.title.y = element_text(size=13),
        legend.title = element_text(size=13),
        legend.text = element_text(size=13)) +
  xlab("") +
  ylab("Delta allele frequency")
dev.off()



pdf("DUSP1_phyloP.pdf", width = 15, height = 3)
ggplot(DUSP1, aes(x=as.factor(Coord), y=phyloPam, group=1)) +
  geom_line(linetype = "solid")+
  geom_point() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size=11), axis.title.y = element_text(size=15)) +
  scale_y_continuous(breaks=seq(-9,9,1)) +
  geom_hline(yintercept=2.27, linetype="solid", 
             color = "forestgreen", size=0.5) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "gray36", size=0.5) +
  ylab("phyloP score")
dev.off()





#####DUSP LD correlation common variants######
setwd("~/Desktop/Pan_Disease/PAPER_VERSION_4/DUSP1_v4/LD/")

####Dprime#####
#data <- read.csv("Pairwise_DUSP1_VarsInLD_4Heatmap_Dprime", sep="\t", head=T)
#head(data)
#mat <- matrix(NA, 10, 10)
#mat[upper.tri(mat)] <- data$Dprime
#colnames(mat) <- c("chr5:172195092","chr5:172196711","chr5:172196752","chr5:172196997","chr5:172197039",
#                   "chr5:172198905","chr5:172198952","chr5:172199489","chr5:172199592","chr5:172199623")
#rownames(mat) <- c("chr5:172195092","chr5:172196711","chr5:172196752","chr5:172196997","chr5:172197039",
#                   "chr5:172198905","chr5:172198952","chr5:172199489","chr5:172199592","chr5:172199623")

#pdf("HeatMap_DUSP1_LD_DPrime.pdf")
#melt(mat) %>%
#  ggplot(aes(Var1, Var2, fill = value)) + geom_tile() +
#  geom_text(aes(label = value), size=3) +
#  scale_fill_gradient2(high ="darkred", limit = c(0,1),name="D'", na.value = "white") +
#  theme_classic() +
#  theme(axis.title = element_blank(),
#        axis.text.x = element_text(size=8, angle=75, hjust=1))
#dev.off()




####Rsquared with VarsID#####
setwd("~/Desktop/Pan_Disease/PAPER_VERSION_4/DUSP1_v4/LD/")
data <- read.csv("Pairwise_DUSP1_VarsInLD_4Heatmap_RSquared_WithVarsIDs", sep="\t", head=T)
head(data)
mat <- matrix(NA, 10, 10)
mat[upper.tri(mat)] <- data$Rsquared

#colnames(mat) <- c("6","11","12","15","16",
#                   "22","23","25","27","28")
#rownames(mat) <- c("6","11","12","15","16",
#                   "22","23","25","27","28")

colnames(mat) <- c("6\n(rs3805476)","11\n(rs2431663)","12\n(rs34471628)","15\n(rs7702178)","16\n(rs35084382)",
                   "22\n(rs881150)","23\n(rs322382)","25\n(rs322381)","27\n(rs322380)","28\n(rs3763067)")
rownames(mat) <- c("6\n(rs3805476)","11\n(rs2431663)","12\n(rs34471628)","15\n(rs7702178)","16\n(rs35084382)",
                   "22\n(rs881150)","23\n(rs322382)","25\n(rs322381)","27\n(rs322380)","28\n(rs3763067)")

#colnames(mat) <- c("6\n(chr5:172195092:G:A)","11\n(chr5:172196711:G:T)","12\n(chr5:172196752:A:G)","15\n(chr5:172196997:T:C)","16\n(chr5:172197039:T:C)",
#                   "22\n(chr5:172198905:T:A)","23\n(chr5:172198952:T:C)","25\n(chr5:172199489:G:A)","27\n(chr5:172199592:T:C)","28\n(chr5:172199623:G:A)")
#rownames(mat) <- c("6\n(chr5:172195092:G:A)","11\n(chr5:172196711:G:T)","12\n(chr5:172196752:A:G)","15\n(chr5:172196997:T:C)","16\n(chr5:172197039:T:C)",
#                   "22\n(chr5:172198905:T:A)","23\n(chr5:172198952:T:C)","25\n(chr5:172199489:G:A)","27\n(chr5:172199592:T:C)","28\n(chr5:172199623:G:A)")

pdf("HeatMap_DUSP1_LD_Rsquared_WithVarsIDs_v4.pdf")
melt(as.data.frame.table(mat)) %>%
  ggplot(aes(as.factor(Var1), as.factor(Var2), fill = value)) + geom_tile() +
  geom_text(aes(label = value), size=2) +
  scale_fill_gradient2(high ="darkred", limit = c(0,1), name="R2", na.value = "white") +
  theme_classic() +
  theme(axis.text.x = element_text(size=8, angle=30, hjust=0.5)) +    ###element_text(size=8, angle=30, hjust=0.5##
  xlab("Variant ID") + ylab("Variant ID")
dev.off()







######################################
######################################

#######BAR PLOT#######
setwd("~/Desktop/Pan_Disease/PAPER_VERSION_4/Check_DUSP1_Variants_Reporter_Assay_FromSergey/Files_and_Plots/")
data <- read.csv("Reporter_Dataframe_4BarPlot_noUTR", sep="\t", head=T)
head(data)

data_OnlyNonStimulated <- subset(data, data$STIMULATION_GROUP == "NON_STIM")
data_OnlyNonStimulated$ALLELE <- factor(data_OnlyNonStimulated$ALLELE, levels=c("Non-Effect", "Effect"))
data_OnlyNonStimulated$REGION_GROUP <- factor(data_OnlyNonStimulated$REGION_GROUP, levels=c("Region#1", "Region#2","Region#3","Region#4","Region#6",
                                                                                                  "Region#7", "Region#8","Region#9","Region#10","Region#11","Region#12",
                                                                                                  "Region#13", "Region#14","Region#15","Region#16","Region#17","Region#18"))


pdf("FigureBarPlot_NonStimulated.pdf", width = 12, height = 10)
ggplot(data_OnlyNonStimulated, aes(x=reorder(HAPLOTYPE, ORDER), y = MEAN, fill = ALLELE)) + 
  geom_col() +
  geom_errorbar(aes(ymin=MEAN-SEM, ymax=MEAN+SEM), width=.2,
                position=position_dodge(.9)) +
#  scale_fill_manual(values=c("grey69","grey30")) +
  scale_fill_manual(values=c("skyblue","slateblue")) +
    facet_grid(~ REGION_GROUP, 
             scales = "free_x",
             space = "free_x",
             switch = "x") +
  theme_classic() +
  theme(panel.spacing = unit(0.20, units = "cm"),
      strip.placement = "outside", # moves the states down
      strip.background = element_rect(fill = "white"),
      strip.text.x = element_text(angle=75, hjust=0.5, size=10),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Reporter fold-change expression") +
  xlab("")
dev.off()


###########UTR############
data <- read.csv("Reporter_Dataframe_4BarPlot_OnlyUTR", sep="\t", head=T)
head(data)

data_OnlyNonStimulated <- subset(data, data$STIMULATION_GROUP == "NON_STIM")
data_OnlyNonStimulated$ALLELE <- factor(data_OnlyNonStimulated$ALLELE, levels=c("Non-Effect", "Effect"))


pdf("FigureBarPlot_NonStimulated_OnlyUTR.pdf")
ggplot(data_OnlyNonStimulated, aes(x=reorder(HAPLOTYPE, ORDER), y = MEAN, fill = ALLELE)) + 
  geom_col() +
  geom_errorbar(aes(ymin=MEAN-SEM, ymax=MEAN+SEM), width=.2,
                position=position_dodge(.9)) +
  
#  geom_signif(data = data_OnlyNonStimulated,
#              aes(y_position=-0.01, xmin=2, xmax=2,
#                  annotations="***"), tip_length=0, manual = T, lwd = 0) +
#  
#                geom_signif(data = data.frame(Group = c("S1","S2")),
#                            aes(y_position=c(5.3, 8.3), xmin=c(0.8, 0.8), xmax=c(1.2, 1.2),
#                                annotations=c("**", "NS")), tip_length=0, manual = T)              
                
                
    #  scale_fill_manual(values=c("grey69","grey30")) +
  scale_fill_manual(values=c("skyblue","slateblue")) +
  facet_grid(~ REGION_GROUP, 
             scales = "free_x",
             space = "free_x",
             switch = "x") +
  theme_classic() +
  theme(panel.spacing = unit(0.20, units = "cm"),
        strip.placement = "outside", # moves the states down
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(angle=75, hjust=0.5, size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("RLU") +
  xlab("")
dev.off()


#########ALL NON_STIMULATED AND STIMULATED
setwd("~/Desktop/Pan_Disease/PAPER_VERSION_4/Check_DUSP1_Variants_Reporter_Assay_FromSergey/Files_and_Plots/")
data <- read.csv("Reporter_Dataframe_4BarPlot_noUTR_BothNonStimAndStim", sep="\t", head=T)


head(data)

data$ALLELE <- factor(data$ALLELE, levels=c("Non-Effect", "Effect"))
data$REGION_GROUP <- factor(data$REGION_GROUP, levels=c("Region#1", "Region#2","Region#3","Region#4","Region#6",
                                                        "Region#7", "Region#8","Region#9","Region#10","Region#11","Region#12",
                                                        "Region#13", "Region#14","Region#15","Region#16","Region#17","Region#18"))

a <- c("black","firebrick","black","firebrick","black","firebrick","black","firebrick","black","firebrick",
       "black","firebrick","black","firebrick","black","firebrick","black","firebrick","black","firebrick",
       "black","firebrick","black","firebrick","black","firebrick","black","firebrick","black","firebrick",
       "black","firebrick","black","firebrick","black","firebrick","black","firebrick","black","firebrick",
       "black","firebrick","black","firebrick","black","firebrick","black","firebrick","black","firebrick",
       "black","firebrick","black","firebrick","black","firebrick","black","firebrick","black","firebrick",
       "black","firebrick","black","firebrick","black","firebrick","black","firebrick","black","firebrick",
       "black","firebrick","black","firebrick","black","firebrick","black","firebrick","black","firebrick",
       "black","firebrick","black","firebrick","black","firebrick","black","firebrick","black","firebrick",
       "black","firebrick")


pdf("FigureBarPlot_NonStimulated_Stimulated.pdf", width = 15, height = 10)
ggplot(data, aes(x=reorder(HAPLOTYPE, ORDER), y = MEAN, fill = ALLELE)) + 
  geom_col() +
  geom_errorbar(aes(ymin=MEAN-SEM, ymax=MEAN+SEM), width=.2,
                position=position_dodge(.9)) +
  #  scale_fill_manual(values=c("grey69","grey30")) +
  scale_fill_manual(values=c("skyblue","slateblue")) +
  facet_grid(~ REGION_GROUP, 
             scales = "free_x",
             space = "free_x",
             switch = "x") +
  theme_classic() +
  theme(panel.spacing = unit(0.20, units = "cm"),
        strip.placement = "outside", # moves the states down
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(angle=75, hjust=0.5, size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour = a)) +
  ylab("Reporter fold-change expression") +
  xlab("")
dev.off()



###########UTR - NON STIMULATED AND STIMULATED############
data <- read.csv("Reporter_Dataframe_4BarPlot_OnlyUTR_BothNonStimAndStim", sep="\t", head=T)
head(data)

data$ALLELE <- factor(data$ALLELE, levels=c("Non-Effect", "Effect"))

a <- c("black","firebrick","black","firebrick","black","firebrick")

pdf("FigureBarPlot_OnlyUTR_NonStimulated_Stimulated.pdf")
ggplot(data, aes(x=reorder(HAPLOTYPE, ORDER), y = MEAN, fill = ALLELE)) + 
  geom_col() +
  geom_errorbar(aes(ymin=MEAN-SEM, ymax=MEAN+SEM), width=.2,
                position=position_dodge(.9)) +
  #  scale_fill_manual(values=c("grey69","grey30")) +
  scale_fill_manual(values=c("skyblue","slateblue")) +
  facet_grid(~ REGION_GROUP, 
             scales = "free_x",
             space = "free_x",
             switch = "x") +
  theme_classic() +
  theme(panel.spacing = unit(0.20, units = "cm"),
        strip.placement = "outside", # moves the states down
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(angle=75, hjust=0.5, size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour = a)) +
  ylab("RLU") +
  xlab("")
dev.off()


######HaCat and HeLa cells expression levels###########
setwd("~/Desktop/Pan_Disease/PAPER_VERSION_4/Check_DUSP1_Variants_Reporter_Assay_FromSergey/Files_and_Plots/")

data <- read.csv("Reporter_Dataframe_4BarPlot_HaCat_HeLa_Cells", head=T, sep="\t")
head(data)

data$Cell <- factor(data$Cell, levels=c("HaCat", "HeLa"))
data$Group <- factor(data$Group, levels=c("Non stimulated", "Stimulated"))

pdf("FigureBarPlot_HaCat_HeLa.pdf")

ggplot(data, aes(Cell, Value, fill = Group)) + geom_bar(stat = "identity", position = 'dodge') + labs(y = "Relative transcript level") + labs(x = "") +
#  scale_fill_manual(values=c("skyblue","slateblue")) +
  scale_fill_manual(values=c("black","firebrick")) +
  theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16), axis.text.y = element_text(size=16),
        legend.title=element_text(size=16), legend.text=element_text(size=16)) +
  labs(fill="")
dev.off()







