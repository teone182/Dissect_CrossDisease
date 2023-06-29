#!/bin/bash -l
#SBATCH -A sens2017142
#SBATCH -p core -n 1
#SBATCH -J ThirdFilter
#SBATCH -t 7-00:00:00
#SBATCH -o /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/output/PanDisease_ThirdFilter-210201.output
#SBATCH -e /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/output/PanDisease_ThirdFilter-210201.error
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

###########STEP	2###############


# Exclude individuals based on laser, king, QC outliers, discordant sex and genotype and clinical
vcftools --gzvcf PanDisease_SNP_PASS_DP_GQ_P_AB_HWE.vcf.gz --remove ALL_SAMPLES_ToBeRemoved_PanDisease_210201.list --recode --recode-INFO-all --out SNP_PASS_QC1 &&

# Zip vcf files
mv SNP_PASS_QC1.recode.vcf SNP_PASS_QC1.vcf &&
bgzip SNP_PASS_QC1.vcf &&

# Index vcf file
tabix -p vcf SNP_PASS_QC1.vcf.gz &&



# Filter variants on max-missing and generate individual missingness file
vcftools --gzvcf SNP_PASS_QC1.vcf.gz --out SNP_PASS_QC1_0.9 --missing-indv --remove-filtered-geno-all --remove-filtered-all --max-missing 0.9 &&

# Filter variants on max-missing and generate new VCF file
vcftools --gzvcf SNP_PASS_QC1.vcf.gz --out SNP_PASS_QC1_0.9 --remove-filtered-geno-all --remove-filtered-all --max-missing 0.9 --recode --recode-INFO-all &&

# Zip vcf files
mv SNP_PASS_QC1_0.9.recode.vcf SNP_PASS_QC1_0.9.vcf &&
bgzip SNP_PASS_QC1_0.9.vcf &&

# Index vcf file
tabix -p vcf SNP_PASS_QC1_0.9.vcf.gz &&



# Get samples with < 80% call rate (to exclude)
awk '$5 >=0.2 {if (NR!=1) print $1 }' SNP_PASS_QC1_0.9.imiss > SNP_PASS_QC1_0.9_imiss_0.2.txt &&

# Remove samples with < 80% call and filter variants on polymorphism
vcftools --gzvcf SNP_PASS_QC1_0.9.vcf.gz --out SNP_PASS_QC1_0.9_0.8_P --remove SNP_PASS_QC1_0.9_imiss_0.2.txt --maf 0.00001 --max-maf 0.99992 --recode --recode-INFO-all &&

# Zip vcf files
mv SNP_PASS_QC1_0.9_0.8_P.recode.vcf SNP_PASS_QC1_0.9_0.8_P.vcf &&
bgzip SNP_PASS_QC1_0.9_0.8_P.vcf &&

# Index vcf file
tabix -p vcf SNP_PASS_QC1_0.9_0.8_P.vcf.gz &&

# Get the AB annotation for your vcf
java -Xmx23g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T VariantAnnotator -R $ref -V SNP_PASS_QC1_0.9_0.8_P.vcf.gz -o SNP_PASS_QC1_0.9_0.8_P_ab.vcf.gz -A AlleleBalance -A AlleleBalanceBySample

# Let me know it is done
echo "Job is done"

date
