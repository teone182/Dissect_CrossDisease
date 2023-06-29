#!/bin/bash -l
#SBATCH -A sens2017142
#SBATCH -p core -n 1
#SBATCH -J Assoc
#SBATCH -t 0-02:00:00
#SBATCH -o /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/008_AllSamples_Cleaned_OverlappingWith2021DatabaseInfo_ASSOCIATION/LogReg_210215.output
#SBATCH -e /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/008_AllSamples_Cleaned_OverlappingWith2021DatabaseInfo_ASSOCIATION/LogReg_210215.error
#SBATCH --mail-user matteo.bianchi@imbim.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load plink/1.90b4.9
module load R
module load vcftools
module load tabix

SourceFolder='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/007_AllSamples_Cleaned_OverlappingWith2021DatabaseInfo'
WorkFolder='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/008_AllSamples_Cleaned_OverlappingWith2021DatabaseInfo_ASSOCIATION'

cd $WorkFolder &&


#####ALL CASES v CONTROLS######maf 0.01#####
plink --vcf $SourceFolder/PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM_Cleaned_OverlappingWithDB_NoFalsePositives.vcf.gz --double-id --allow-no-sex --maf 0.01 --logistic --covar $SourceFolder/AllSamples_Cleaned_OverlappingWithDb_PanDisease_3544_PCs_Sex_Cohort_DiseaseStatus --covar-number 1,2,3,4,11 --adjust -ci 0.95 --pheno $SourceFolder/AllSamples_Cleaned_OverlappingWithDb_PanDisease_3544_CasesControls_PhenoFile.pheno --out AllSamples_Cleaned_OverlappingWithDb_PanDisease_3544_NoFalsePositives_maf001_LogReg4PCASex &&

grep 'ADD' AllSamples_Cleaned_OverlappingWithDb_PanDisease_3544_NoFalsePositives_maf001_LogReg4PCASex.assoc.logistic | sort -g -k12,12 > AllSamples_Cleaned_OverlappingWithDb_PanDisease_3544_NoFalsePositives_maf001_LogReg4PCASex.assoc.logistic_sorted &&


#######Calculate frequencies##########
plink --vcf $SourceFolder/PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM_Cleaned_OverlappingWithDB_NoFalsePositives.vcf.gz --double-id --allow-no-sex --pheno $SourceFolder/AllSamples_Cleaned_OverlappingWithDb_PanDisease_3544_CasesControls_PhenoFile.pheno --out PanDiseaseCleaned_3544_DiffFreq --freq case-control



