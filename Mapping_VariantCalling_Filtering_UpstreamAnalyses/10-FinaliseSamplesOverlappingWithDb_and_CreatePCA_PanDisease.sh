#!/bin/bash -l
#SBATCH -A sens2017142
#SBATCH -p core -n 1
#SBATCH -J MakeFiles
#SBATCH -t 3-00:00:00
#SBATCH -o /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/output/CreateFinalFile_and_PCA_210205.output
#SBATCH -e /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/output/CreateFinalFile_and_PCA_210205.error
#SBATCH --mail-user matteo.bianchi@imbim.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load plink/1.90b4.9
module load R
module load vcftools
module load tabix

WorkFolder='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/007_AllSamples_Cleaned_OverlappingWith2021DatabaseInfo'
OriginalFolder='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/006_AllSamples_VQSR'

cd $WorkFolder &&

#######Create a final file with only samples overlapping the Db##############
vcftools --gzvcf $OriginalFolder/PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM.vcf.gz --keep AllSamples_Cleaned_OverlappingWithDb_PanDisease_3544.list --recode --recode-INFO-all --out PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM_Cleaned_OverlappingWithDB &&

mv PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM_Cleaned_OverlappingWithDB.recode.vcf PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM_Cleaned_OverlappingWithDB.vcf &&

bgzip PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM_Cleaned_OverlappingWithDB.vcf &&
tabix -p vcf PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM_Cleaned_OverlappingWithDB.vcf.gz &&


###########To Find FalsePositve SNPs#######
#plink --vcf PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM_Cleaned_OverlappingWithDB.vcf.gz --double-id --allow-no-sex --freq case-control --pheno AllSamples_Cleaned_OverlappingWithDb_PanDisease_3544_CasesControls_PhenoFile.pheno --out PanDiseaseCleaned_3544_DiffFreq

####R######
#data <- read.table("PanDiseaseCleaned_3544_DiffFreq.frq.cc", head=T)
#rem <- subset(data, (data$MAF_A == 0 & data$MAF_U > 0.01) | data$MAF_U == 0 & data$MAF_A > 0.01)   ####remove those variants with 0 maf in one group and greater than 0.01 in the other group#######
#rem
#     CHR           SNP A1 A2 MAF_A   MAF_U NCHROBS_A NCHROBS_U
#6623   1 chr1:23141540  T  G     0 0.02004      4314      2246
#9829   1 chr1:40091280  T  G     0 0.03196      4466      2378
#####################
###########+++++++++++++++++++++++#######


#######REMOVE FALSE POSITIVES SNPs##############
vcftools --gzvcf PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM_Cleaned_OverlappingWithDB.vcf.gz --exclude FalsePositiveSNPs.list --recode --recode-INFO-all --out PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM_Cleaned_OverlappingWithDB_NoFalsePositives &&

mv PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM_Cleaned_OverlappingWithDB_NoFalsePositives.recode.vcf PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM_Cleaned_OverlappingWithDB_NoFalsePositives.vcf &&

bgzip PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM_Cleaned_OverlappingWithDB_NoFalsePositives.vcf &&
tabix -p vcf PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM_Cleaned_OverlappingWithDB_NoFalsePositives.vcf.gz &&


#####Create PCA########
plink --vcf PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM_Cleaned_OverlappingWithDB_NoFalsePositives.vcf.gz --allow-no-sex --double-id --exclude range price_high_ld_regions.txt --make-bed --out PanDiseaseCleaned_3544_high_ld_excl_noFP &&

plink --bfile PanDiseaseCleaned_3544_high_ld_excl_noFP --allow-no-sex --maf 0.05 --double-id --indep-pairwise 1500 150 0.2 --out PanDiseaseCleaned_3544_high_ld_excl_indep_0.2_noFP &&

plink --bfile PanDiseaseCleaned_3544_high_ld_excl_noFP --allow-no-sex --double-id --exclude PanDiseaseCleaned_3544_high_ld_excl_indep_0.2_noFP.prune.out --make-bed --out PanDiseaseCleaned_3544_high_ld_excl_indep_0.2_noFP &&

######get the PCA######
plink --bfile PanDiseaseCleaned_3544_high_ld_excl_indep_0.2_noFP --allow-no-sex --maf 0.05 --double-id --pca 10 --out PanDiseaseCleaned_3544_MAF5R202Pruned_PCA_NoFP

#######Get various stats#######
plink --vcf PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM_Cleaned_OverlappingWithDB_NoFalsePositives.vcf.gz --missing --double-id --allow-no-sex --out PanDiseaseCleaned_3544_Missingness
