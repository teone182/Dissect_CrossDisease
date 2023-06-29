#!/bin/bash -l
#SBATCH -A sens2017142
#SBATCH -p core -n 1
#SBATCH -J table
#SBATCH -t 0-01:00:00
#SBATCH -o /proj/sens2017142/nobackup/matteob/SCRIPTS/02_PAN_DISEASE_PROJECT/40_GeneRiskScore_Scripts/CreateLists.output
#SBATCH -e /proj/sens2017142/nobackup/matteob/SCRIPTS/02_PAN_DISEASE_PROJECT/40_GeneRiskScore_Scripts/CreateLists.error
#SBATCH --mail-user matteo.bianchi@imbim.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load R_packages/4.0.0


WorkFolder='/proj/sens2017142/nobackup/matteob/SCRIPTS/02_PAN_DISEASE_PROJECT/40_GeneRiskScore_Scripts' &&
SourceFolder='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/007_AllSamples_Cleaned_OverlappingWith2021DatabaseInfo' &&
Scripts='/proj/sens2017142/nobackup/matteob/SCRIPTS/02_PAN_DISEASE_PROJECT/40_GeneRiskScore_Scripts' &&

#####Create the input file for GLM###############

cd $WorkFolder

R --vanilla <<EOF

a <- read.csv("$SourceFolder/AllSamples_Cleaned_OverlappingWithDb_PanDisease_3544_PCs_Sex_Cohort_DiseaseStatus_Origin_SubjectID_AgeAtDiagnosis_ALL_SEROLOGICAL_CLINICAL_SUBPHENOTYPES_BETTER", sep="\t", head=T)

controls <- subset(a, a\$DiseaseStatus == 1)[2]

write.table(controls, "List_PanDisease_CONTROLS", quote=F, row.names=F, col.names=F)

EOF

cd $WorkFolder

while read pheno
do

cd $WorkFolder

R --vanilla <<EOF

a <- read.csv("$SourceFolder/AllSamples_Cleaned_OverlappingWithDb_PanDisease_3544_PCs_Sex_Cohort_DiseaseStatus_Origin_SubjectID_AgeAtDiagnosis_ALL_SEROLOGICAL_CLINICAL_SUBPHENOTYPES_BETTER", sep="\t", head=T)

${pheno}_yes <- subset(a, a\$${pheno} == "YES")[2]
${pheno}_no <- subset(a, a\$${pheno} == "NO")[2]
${pheno}_no_and_miss <- subset(a, a\$${pheno} == "NO" | (is.na(a\$${pheno}) & a\$DiseaseStatus == 2))[2] 

write.table(${pheno}_yes, "List_PanDisease_CASES_${pheno}_POSITIVE", quote=F, row.names=F, col.names=F)
write.table(${pheno}_no, "List_PanDisease_CASES_${pheno}_NEGATIVE", quote=F, row.names=F, col.names=F)
write.table(${pheno}_no_and_miss, "List_PanDisease_CASES_${pheno}_NEGATIVE_and_MISSING", quote=F, row.names=F, col.names=F)

EOF

done<$Scripts/List_Subphenos_Groups


