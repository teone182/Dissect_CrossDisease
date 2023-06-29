#!/bin/bash -l
#SBATCH -A sens2017142
#SBATCH -p core -n 1
#SBATCH -J GLM
#SBATCH -t 2-00:00:00
#SBATCH -o /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/040_Cleaned_OverlappingWith2021DatabaseInfo_GenePy_+-2Kb_GeneRiskScore/03_PanDisease_GeneRiskScore_2Kb_Downstream_GLMSummarise_ALL_SUBPHENOTYPES_WithAllAAb_BETTER_221126.output
#SBATCH -e /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/040_Cleaned_OverlappingWith2021DatabaseInfo_GenePy_+-2Kb_GeneRiskScore/03_PanDisease_GeneRiskScore_2Kb_Downstream_GLMSummarise_ALL_SUBPHENOTYPES_WithAllAAb_BETTER_221126.error
#SBATCH --mail-user matteo.bianchi@imbim.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load R_packages/4.0.0

WorkFolder='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/040_Cleaned_OverlappingWith2021DatabaseInfo_GenePy_+-2Kb_GeneRiskScore/Results/GLM' &&
SourceFolder='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/040_Cleaned_OverlappingWith2021DatabaseInfo_GenePy_+-2Kb_GeneRiskScore/Results/GLM' &&
Scripts='/proj/sens2017142/nobackup/matteob/SCRIPTS/02_PAN_DISEASE_PROJECT/40_GeneRiskScore_Scripts' &&

cd $WorkFolder

######Replace TotalScore with Gene name
while read cases controls
do
while read chrom start end gene
do

cd $WorkFolder

sed "s/TotalScore/${gene}/g" $SourceFolder/LogitModel_PanDisease_${cases}_${controls}_${chrom}_${gene}.txt > $WorkFolder/LogitModel_PanDisease_${cases}_${controls}_${chrom}_${gene}_OK.txt

done<$Scripts/TargetGenesSpaceCoordinates_PlusMinus2Kb_NEW_1847.bed
done<$Scripts/List_ALL_SUBPHENOTYPES


######concatenate all GLM results

cd $WorkFolder

while read cases controls
do

cd $WorkFolder

cat LogitModel_PanDisease_${cases}_${controls}_chr*_OK.txt > LogitModels_PanDisease_${cases}_vs_${controls}.txt

rm LogitModel_PanDisease_${cases}_${controls}_chr*_OK.txt

done<$Scripts/List_ALL_SUBPHENOTYPES


cd $WorkFolder

while read case_group control_group
do

cd $WorkFolder

#####GLM RESULTS#####
grep -wF "Pr(>|z|)" -A 2 LogitModels_PanDisease_${case_group}_vs_${control_group}.txt | grep -vF "(Intercept)" | grep -v "^--" | sed "s/data_merged\\$//g" > LogitModels_PanDisease_${case_group}_vs_${control_group}_AllVariantsFinal.txt &&

grep -v 'Estimate' LogitModels_PanDisease_${case_group}_vs_${control_group}_AllVariantsFinal.txt | grep -v "\sNA\s" > LogitModels_PanDisease_${case_group}_vs_${control_group}_AllVariantsFinal_OK.txt &&


sed 's/< //g' LogitModels_PanDisease_${case_group}_vs_${control_group}_AllVariantsFinal_OK.txt | awk -v OFS="\t" '{print $1,$5}' > LogitModels_PanDisease_${case_group}_vs_${control_group}_4FDR &&

R --vanilla <<EOF
data <- read.csv("LogitModels_PanDisease_${case_group}_vs_${control_group}_4FDR", head=F, sep="\t")
data\$Bonferroni <- p.adjust(data\$V2, method="bonferroni")
data\$FDR <- p.adjust(data\$V2, method="BH")
write.table(data, "LogitModels_PanDisease_${case_group}_vs_${control_group}_AllVariants_Bonferroni_FDR.txt", quote=F, sep="\t", row.names=F, col.names=F)
EOF

sed '1iGENE\tP\tBonferroni\tFDR' LogitModels_PanDisease_${case_group}_vs_${control_group}_AllVariants_Bonferroni_FDR.txt | sort -gk4,4 > LogitModels_PanDisease_${case_group}_vs_${control_group}_AllVariants_Bonferroni_FDR.txt_sorted

rm LogitModels_PanDisease_${case_group}_vs_${control_group}_4FDR

done<$Scripts/List_ALL_SUBPHENOTYPES

