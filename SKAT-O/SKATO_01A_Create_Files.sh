#!/bin/bash
#SBATCH -A sens2017142
#SBATCH -p core -n 1
#SBATCH -J CreateFiles
#SBATCH -t 0-06:00:00
#SBATCH -o /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/060_SKATO/Output_Error/PrepareFiles.output
#SBATCH -e /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/060_SKATO/Output_Error/PrepareFiles.error
#SBATCH --mail-user matteo.bianchi@imbim.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load plink/1.90b4.9
module load vcftools
module load bcftools
module load R_packages/4.0.0
module load tabix

SourceFolderPhenos='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/007_AllSamples_Cleaned_OverlappingWith2021DatabaseInfo'
SourceFolderVars='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/040_Cleaned_OverlappingWith2021DatabaseInfo_+-2Kb_GeneRiskScore'
WorkFolder='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/060_SKATO/PrepareFiles'
SetID_File_Folder='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/060_SKATO/PrepareFiles/SetID_File_Folder'
Scripts='/proj/sens2017142/nobackup/matteob/SCRIPTS/02_PAN_DISEASE_PROJECT/41_SKATO_Scripts'


cd $WorkFolder

######Create plink bim,bed,fam

plink -vcf $SourceFolderVars/PanDisease_Annotated.vcf.gz --double-id --allow-no-sex --make-bed --out PanDisease 

#######Create covariate file#####
cut -f 2-7,14 $SourceFolderPhenos/AllSamples_Cleaned_OverlappingWithDb_PanDisease_3544_PCs_Sex_Cohort_DiseaseStatus_Origin_SubjectID_AgeAtDiagnosis_ALL_SEROLOGICAL_CLINICAL_SUBPHENOTYPES_BETTER > PanDisease_COVARIATES_File   ####NOTE THE HEADER###

######Create setID file with all the gene/variants to be tested######
while read chrom start end gene
do

cd $SetID_File_Folder

grep -e "\s$gene\s" -e "\s$gene;" -e ";$gene\s" -e ";$gene;" $SourceFolderVars/PanDisease_${chrom}.hg19_multianno.txt | grep -v "intergenic" | awk -v FS="\t" -v OFS="" -v var="$gene" '{print var,"\t",$1,":",$2}' > Vars_${gene}.list

done<$Scripts/TargetGenesSpaceCoordinates_PlusMinus2Kb_NEW_1847.bed

cd $WorkFolder

cat $SetID_File_Folder/Vars_*.list > PanDisease_SetID_4SKATO
#rm $SetID_File_Folder/*


########Create Matrix with subphenotype comparisons


cd $WorkFolder

R --vanilla <<EOF

a <- read.csv("$SourceFolderPhenos/AllSamples_Cleaned_OverlappingWithDb_PanDisease_3544_PCs_Sex_Cohort_DiseaseStatus_Origin_SubjectID_AgeAtDiagnosis_ALL_SEROLOGICAL_CLINICAL_SUBPHENOTYPES_BETTER", sep="\t", head=T)

controls <- subset(a, a\$DiseaseStatus == 1)[2]
write.table(controls, "List_PanDisease_CONTROLS", quote=F, row.names=F, col.names=F)

cases <- subset(a, a\$DiseaseStatus == 2)[2]
write.table(cases, "List_PanDisease_CASES", quote=F, row.names=F, col.names=F)

EOF

cd $WorkFolder

while read pheno
do

cd $WorkFolder

R --vanilla <<EOF

a <- read.csv("$SourceFolderPhenos/AllSamples_Cleaned_OverlappingWithDb_PanDisease_3544_PCs_Sex_Cohort_DiseaseStatus_Origin_SubjectID_AgeAtDiagnosis_ALL_SEROLOGICAL_CLINICAL_SUBPHENOTYPES_BETTER", sep="\t", head=T)

${pheno}_yes <- subset(a, a\$${pheno} == "YES")[2]
${pheno}_no <- subset(a, a\$${pheno} == "NO")[2]

write.table(${pheno}_yes, "List_PanDisease_CASES_${pheno}_POSITIVE", quote=F, row.names=F, col.names=F)
write.table(${pheno}_no, "List_PanDisease_CASES_${pheno}_NEGATIVE", quote=F, row.names=F, col.names=F)

EOF

done<$Scripts/List_Subphenos_Groups


cd $WorkFolder

while read cases controls
do

cd $WorkFolder

awk -v OFS="\t" '{print $1,$1,"2"}' List_PanDisease_${cases} > List_PanDisease_${cases}_TMP
awk -v OFS="\t" '{print $1,$1,"1"}' List_PanDisease_${controls} > List_PanDisease_${controls}_TMP
cat List_PanDisease_${cases}_TMP List_PanDisease_${controls}_TMP | sed "1iFID\tIID\t${cases}_VS_${controls}" > Data_${cases}_VS_${controls}

rm List_PanDisease_${cases}_TMP List_PanDisease_${controls}_TMP

done<$Scripts/List_SUBPHENOTYPES_Prep4SkatO

cd $WorkFolder

R --vanilla <<EOF

multMerge_data = function(mypath){
  filenames = list.files(path = mypath, pattern=glob2rx("Data_CASES_*"), full.names = TRUE)
  datalist = lapply(filenames, 
      function(x){read.csv(file = x, header = TRUE, sep="\t", stringsAsFactors = FALSE)})

  Reduce(function(x,y) {merge(x, y, by=c("FID","IID"), all=T)}, datalist)
}

data_merged <- multMerge_data("$WorkFolder/")

write.table(data_merged, "PanDisease_ALLSubPhenotypesMatrix.txt", quote=F, row.names=F, sep="\t")

EOF




