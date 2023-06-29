#!/bin/bash -l
#SBATCH -A sens2017142
#SBATCH -p core -n 4
#SBATCH -J filter2_stats_scand_myo
#SBATCH -t 4-00:00:00
#SBATCH -o /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/output/PanDisease_SecondFilter_and_GetStats-210119.output
#SBATCH -e /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/output/PanDisease_SecondFilter_and_GetStats-210119.error
#SBATCH --mail-user matteo.bianchi@imbim.uu.se
#SBATCH --mail-type=ALL

# Load modules
module load bioinfo-tools
module load vcftools
module load bcftools
module load tabix
module load samtools
module load R
module load plink
module load plinkseq

ref='/proj/sens2017142/sens2017107/DISSECT/Reference/homo_sapiens_bwa.0.7.12.fa'

hapmap='/proj/sens2017142/sens2017107/DISSECT/Reference/hapmap.vcf'
omni='/proj/sens2017142/sens2017107/DISSECT/Reference/omni.vcf'
snphc='/proj/sens2017142/sens2017107/DISSECT/Reference/1000Gsnp.vcf'
dbsnp='/proj/sens2017142/sens2017107/DISSECT/Reference/dbSNP_147.vcf'
mills='/proj/sens2017142/sens2017107/DISSECT/Reference/mills.vcf'

InDir='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/005_AllSamples_JointGenotyping'
OutDir='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/006_AllSamples_VQSR'

cd $OutDir

# Calculate HWE p-values based on only controls
vcftools --gzvcf PanDisease_SNP_PASS_DP_GQ_P_AB.vcf.gz --keep PanDisease_AllControlsBeforeSampleFiltering.list --hardy --out PanDisease_SNP_PASS_DP_GQ_P_AB_controls &&

# Get positions that have HWE p-value>0.000001 (include those)
awk '($6 + 0) >= 0.000001 {if (NR!=1) print $1 "\t" $2}' PanDisease_SNP_PASS_DP_GQ_P_AB_controls.hwe > PanDisease_SNP_PASS_DP_GQ_P_AB_controls_hwe_pass.txt &&

# Filter on HWE (include positions that pass HWE)
vcftools --gzvcf PanDisease_SNP_PASS_DP_GQ_P_AB.vcf.gz --positions PanDisease_SNP_PASS_DP_GQ_P_AB_controls_hwe_pass.txt  --recode --recode-INFO-all --out PanDisease_SNP_PASS_DP_GQ_P_AB_HWE &&

# Zip vcf files
mv PanDisease_SNP_PASS_DP_GQ_P_AB_HWE.recode.vcf PanDisease_SNP_PASS_DP_GQ_P_AB_HWE.vcf &&
bgzip PanDisease_SNP_PASS_DP_GQ_P_AB_HWE.vcf &&

# Index vcf file
tabix -p vcf PanDisease_SNP_PASS_DP_GQ_P_AB_HWE.vcf.gz &&


# Generate plinkseq i-stat
pseq PanDisease_SNP_PASS_DP_GQ_P_AB_HWE.vcf.gz i-stats > PanDisease_snp_pass_dp_gq_p_ab_hwe.istat &&

#  For chr X

# Extract chr X only
vcftools --gzvcf PanDisease_SNP_PASS_DP_GQ_P_AB_HWE.vcf.gz --chr "chrX" --recode --recode-INFO-all --out PanDisease_SNP_PASS_DP_GQ_P_AB_HWE_chrX &&

# Zip vcf files
mv PanDisease_SNP_PASS_DP_GQ_P_AB_HWE_chrX.recode.vcf PanDisease_SNP_PASS_DP_GQ_P_AB_HWE_chrX.vcf &&
bgzip PanDisease_SNP_PASS_DP_GQ_P_AB_HWE_chrX.vcf &&

# Index vcf file
tabix -p vcf PanDisease_SNP_PASS_DP_GQ_P_AB_HWE_chrX.vcf.gz &&

# Generate plinkseq i-stat
pseq PanDisease_SNP_PASS_DP_GQ_P_AB_HWE_chrX.vcf.gz i-stats > PanDisease_snp_pass_dp_gq_p_ab_hwe_chrx.istat

########CHECK THE STATS FILES FOR PINPOINTING SAMPLES THAT NEED TO BE REMOVED####################

# Let me know it is done
echo "Job is done"

date
