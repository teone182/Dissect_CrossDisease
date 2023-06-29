#!/bin/bash

while read line; do

#Set variables
chromosome=$line
PhenoFolder='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/007_AllSamples_Cleaned_OverlappingWith2021DatabaseInfo'
SourceFolder='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/040_Cleaned_OverlappingWith2021DatabaseInfo_GenePy_+-2Kb_GeneRiskScore'
WorkFolder='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/040_Cleaned_OverlappingWith2021DatabaseInfo_GenePy_+-2Kb_GeneRiskScore/Results'
GLMFolder='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/040_Cleaned_OverlappingWith2021DatabaseInfo_GenePy_+-2Kb_GeneRiskScore/Results/GLM'
Scripts='/proj/sens2017142/nobackup/matteob/SCRIPTS/02_PAN_DISEASE_PROJECT/40_GeneRiskScore_Scripts'

echo "#!/bin/bash -l" > 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "#SBATCH -A sens2017142" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "#SBATCH -p core -n 1" >> 02_PanDisease_GeneRiskScore_Downstream.sc;  ###when you generate Annovar files, use 4 cores#####
echo "#SBATCH -J PanDisease_Dowmstream_AllSubPhenos_${chromosome}" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "#SBATCH -t 0-100:00:00" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "#SBATCH -o /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/040_Cleaned_OverlappingWith2021DatabaseInfo_GenePy_+-2Kb_GeneRiskScore/02_PanDisease_GeneRiskScore_2Kb_Downstream_${chromosome}.output" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "#SBATCH -e /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/040_Cleaned_OverlappingWith2021DatabaseInfo_GenePy_+-2Kb_GeneRiskScore/02_PanDisease_GeneRiskScore_2Kb_Downstream_${chromosome}.error" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "#SBATCH --mail-user matteo.bianchi@imbim.uu.se" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "#SBATCH --mail-type=ALL" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

echo "module load bioinfo-tools" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "module load vcftools" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "module load bcftools" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "module load snpEff" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "module load python/2.7.9" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "module load R_packages/4.0.0" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "module load tabix" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "module load BEDTools" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

echo "chromosome=$chromosome" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "PhenoFolder=$PhenoFolder" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "SourceFolder=$SourceFolder" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "WorkFolder=$WorkFolder" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "GLMFolder=$GLMFolder" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "Scripts=$Scripts" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

#####Create gene subsets files######

echo "awk '\$1 == \"$chromosome\" {print}' $Scripts/TargetGenesSpaceCoordinates_PlusMinus2Kb_NEW_1847.bed > $Scripts/List_of_Genes_${chromosome}" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

echo "cd $WorkFolder" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

#######Fix unwanted spaces########
echo "sed 's/ /_/g' $SourceFolder/PanDisease_${chromosome}.meta  > $WorkFolder/PanDisease_${chromosome}.meta_SpacesFixed" >> 02_PanDisease_GeneRiskScore_Downstream.sc

echo "$Scripts/GeneRiskScore_WithLinearConversion.pl $WorkFolder/PanDisease_${chromosome}.meta_SpacesFixed $WorkFolder/PanDisease_${chromosome}_scores.meta" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

echo "rm $WorkFolder/PanDisease_${chromosome}.meta_SpacesFixed" >> 02_PanDisease_GeneRiskScore_Downstream.sc

echo "while read chr start end gene" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "do"  >> 02_PanDisease_GeneRiskScore_Downstream.sc;

echo "cd $WorkFolder" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

####IF WE WANT FOR EVERY GENE ONLY THE VARIANTS IN UTRs, INTRONS AND INTERGENIC REGIONS (100Kb) THEN WE USE THE FOLLOWING:
echo "grep -e \"^Chr\" -e \"\s\${gene}\s\" -e \"\s\${gene};\" -e \";\${gene}\s\" -e \";\${gene};\" $WorkFolder/PanDisease_${chromosome}_scores.meta | grep -v \"intergenic\" > PanDisease_${chromosome}_scores.meta_\${gene}" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

###get only the variant scores#####
echo "cut -f 13- PanDisease_${chromosome}_scores.meta_\${gene} > VariantScores_AllVariants_${chromosome}_\${gene}.txt" >> 02_PanDisease_GeneRiskScore_Downstream.sc;


##############
######We need to have lists of samples: for example cases, controls, subgroups, etc.
#List_PanDisease_CASES_ALL
#List_PanDisease_CONTROLS_ALL

echo "while read cases controls" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "do" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

echo "cd $WorkFolder" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

#######
echo "R --vanilla <<EOF" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "library(dplyr)" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

echo "data <- read.csv(\"VariantScores_AllVariants_${chromosome}_\${gene}.txt\", sep=\"\\t\", head=T, stringsAsFactors = FALSE)" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "\${cases} <- read.csv(\"$Scripts/List_PanDisease_\${cases}\", head=F)" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "\${controls} <- read.csv(\"$Scripts/List_PanDisease_\${controls}\", head=F)" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "\${cases}_names <- as.vector(\${cases}\\\$V1)" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "\${controls}_names <- as.vector(\${controls}\\\$V1)" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "\${cases}_names_ok <- gsub(\"-\", \".\", \${cases}_names)" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "\${controls}_names_ok <- gsub(\"-\", \".\", \${controls}_names)" >> 02_PanDisease_GeneRiskScore_Downstream.sc;


echo "data_\${cases} <- data %>% select(all_of(\${cases}_names_ok))" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "data_\${cases}_fixed <- data.frame(apply(data_\${cases}, 2, function(x) as.numeric(as.character(x))))" >> 02_PanDisease_GeneRiskScore_Downstream.sc
echo "data_\${cases}_transposed <- as.data.frame(t(data_\${cases}_fixed))" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "data_\${cases}_transposed\\\$TotalScore <- rowSums(data_\${cases}_transposed, na.rm=T)" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "data_\${cases}_transposed\\\$IID <- rownames(data_\${cases}_transposed)" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "data_\${cases}_final <- data.frame(data_\${cases}_transposed\\\$IID, data_\${cases}_transposed\\\$TotalScore)" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "colnames(data_\${cases}_final) <- c(\"IID\", \"TotalScore\")" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

echo "data_\${controls} <- data %>% select(all_of(\${controls}_names_ok))" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "data_\${controls}_fixed <- data.frame(apply(data_\${controls}, 2, function(x) as.numeric(as.character(x))))" >> 02_PanDisease_GeneRiskScore_Downstream.sc
echo "data_\${controls}_transposed <- as.data.frame(t(data_\${controls}_fixed))" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "data_\${controls}_transposed\\\$TotalScore <- rowSums(data_\${controls}_transposed, na.rm=T)" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "data_\${controls}_transposed\\\$IID <- rownames(data_\${controls}_transposed)" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "data_\${controls}_final <- data.frame(data_\${controls}_transposed\\\$IID, data_\${controls}_transposed\\\$TotalScore)" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "colnames(data_\${controls}_final) <- c(\"IID\", \"TotalScore\")" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

echo "total <- rbind(data_\${cases}_final, data_\${controls}_final)" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

echo "write.table(total, \"Total_${chromosome}_\${gene}.txt\", quote=F, sep=\"\t\", row.names=F)" >> 02_PanDisease_GeneRiskScore_Downstream.sc;  ################JUST TO CHECK####

#######Stats test - GLM######

##################Append cohort Info####################
echo "IdsInfo <- read.csv(\"$PhenoFolder/AllSamples_Cleaned_OverlappingWithDb_PanDisease_3544_PCs_Sex_Cohort_DiseaseStatus\", sep=\"\\t\", header=T)" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "IdsInfo\\\$IID <- gsub('-', '.', IdsInfo\\\$IID)" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

echo "data_merged <- merge(total, IdsInfo, by=\"IID\")" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "data_merged\\\$BinaryStatus <- as.numeric(ifelse(data_merged\\\$DiseaseStatus == \"1\", \"0\", \"1\"))" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

########GLM ALL - cases vs controls#####

echo "logit_model <- glm(data_merged\\\$BinaryStatus ~ data_merged\\\$TotalScore+data_merged\\\$PC1+data_merged\\\$PC2+data_merged\\\$PC3+data_merged\\\$PC4+data_merged\\\$Sex, data = data_merged, family=\"binomial\")" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

echo "options(max.print=5000)" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

echo "sink(\"$GLMFolder/LogitModel_PanDisease_\${cases}_\${controls}_${chromosome}_\${gene}.txt\")" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "summary(logit_model)" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "sink()" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

echo "EOF" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

echo "done<$Scripts/List_ALL_SUBPHENOTYPES" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "done<$Scripts/List_of_Genes_${chromosome}" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

echo "rm $Scripts/List_of_Genes_${chromosome}" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "rm $WorkFolder/VariantScores_AllVariants_${chromosome}_*" >> 02_PanDisease_GeneRiskScore_Downstream.sc;
echo "rm $WorkFolder/PanDisease_${chromosome}_*" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

#Let me know it is done
echo "echo "Job is done"" >> 02_PanDisease_GeneRiskScore_Downstream.sc;

#sbatch the source file
sbatch 02_PanDisease_GeneRiskScore_Downstream.sc;

#Must be the same as the file input when calling the script
done <List_chromosomes.list

