#!/bin/bash
#SBATCH -A sens2017142
#SBATCH -p core -n 1
#SBATCH -J CreateFiles
#SBATCH -t 0-01:00:00
#SBATCH -o /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/060_SKATO/Output_Error/PrepareFiles_CaddWeights.output
#SBATCH -e /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/060_SKATO/Output_Error/PrepareFiles_caddWeights.error
#SBATCH --mail-user matteo.bianchi@imbim.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load plink/1.90b4.9
module load vcftools
module load bcftools
module load R_packages/4.0.0
module load tabix

SourceFolderVars='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/040_Cleaned_OverlappingWith2021DatabaseInfo_+-2Kb_GeneRiskScore'
WorkFolder='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/060_SKATO/PrepareFiles'
CaddWeights_File_Folder='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/060_SKATO/PrepareFiles/CaddWeights_File_Folder'
Scripts='/proj/sens2017142/nobackup/matteob/SCRIPTS/02_PAN_DISEASE_PROJECT/41_SKATO_Scripts'


cd $WorkFolder

######Create setID file with all the gene/variants to be tested######
while read chrom start end gene
do

cd $CaddWeights_File_Folder

grep -e "\s$gene\s" -e "\s$gene;" -e ";$gene\s" -e ";$gene;" $SourceFolderVars/PanDisease_${chrom}.hg19_multianno.txt | grep -v "intergenic" | awk -v FS="\t" -v OFS="" '{print $1,":",$2,"\t",$20}' > Vars_${gene}_CaddWeights.list

done<$Scripts/TargetGenesSpaceCoordinates_PlusMinus2Kb_NEW_1847.bed

cd $WorkFolder

cat $CaddWeights_File_Folder/Vars_*_CaddWeights.list > PanDisease_CaddWeights_4SKAT
rm $CaddWeights_File_Folder/*

$Scripts/CADDRaw_LinearConversion.pl PanDisease_CaddWeights_4SKAT PanDisease_CaddWeights_WithLinearConversion_4SKAT


