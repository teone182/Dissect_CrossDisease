#!/bin/bash -l
#SBATCH -A sens2017142
#SBATCH -p core -n 1
#SBATCH -J FourthFilter
#SBATCH -t 7-00:00:00
#SBATCH -o /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/output/PanDisease_FourthFilter-210202.output
#SBATCH -e /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/output/PanDisease_FourthFilter-210202.error
#SBATCH --mail-user matteo.bianchi@imbim.uu.se
#SBATCH --mail-type=ALL

# Load modules
module load bioinfo-tools
module load vcftools
module load bcftools
module load tabix
module load samtools
module load R
module load plink/1.90b4.9

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
vcftools --gzvcf SNP_PASS_QC1_0.9_0.8_P_ab.vcf.gz --keep PanDisease_AllControlsAFTERSampleFiltering.list --hardy --out SNP_PASS_QC1_0.9_0.8_P_ab_controls &&

# Get positions that have HWE p-value>0.000001 (include those)
awk '($6 + 0) >= 0.000001 {if (NR!=1) print $1 "\t" $2}' SNP_PASS_QC1_0.9_0.8_P_ab_controls.hwe > PanDisease_SNP_PASS_QC1_0.9_0.8_P_ab_controls_hwe_pass.txt &&

# Filter on HWE (include positions that pass HWE)
vcftools --gzvcf SNP_PASS_QC1_0.9_0.8_P_ab.vcf.gz --positions PanDisease_SNP_PASS_QC1_0.9_0.8_P_ab_controls_hwe_pass.txt  --recode --recode-INFO-all --out SNP_PASS_QC1_0.9_0.8_P_ab_HWE &&

# Zip vcf files
mv SNP_PASS_QC1_0.9_0.8_P_ab_HWE.recode.vcf SNP_PASS_QC1_0.9_0.8_P_ab_HWE.vcf &&
bgzip SNP_PASS_QC1_0.9_0.8_P_ab_HWE.vcf &&

# Index vcf file
tabix -p vcf SNP_PASS_QC1_0.9_0.8_P_ab_HWE.vcf.gz &&




#########################


# For final filtered vcf file


# Update AN, AF, AC field
zcat SNP_PASS_QC1_0.9_0.8_P_ab_HWE.vcf.gz | fill-an-ac | bgzip -c > SNP_PASS_QC_0.9_0.8_u.vcf.gz &&


# Index vcf file
tabix -p vcf SNP_PASS_QC_0.9_0.8_u.vcf.gz &&


# Annotate with RS number
java -Xmx23g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T VariantAnnotator -R $ref -V SNP_PASS_QC_0.9_0.8_u.vcf.gz --dbsnp $dbsnp -o SNP_PASS_QC_0.9_0.8_u_rs.vcf.gz &&

# Exchange missing RS number by CHR:POS - DOES NOT WORK, RUN SEPARATELY
bcftools annotate --set-id +'%CHROM\:%POS\' SNP_PASS_QC_0.9_0.8_u_rs.vcf.gz -Oz -o SNP_PASS_QC_0.9_0.8_u_rs2.vcf.gz &&


# Index vcf file
tabix -p vcf SNP_PASS_QC_0.9_0.8_u_rs2.vcf.gz &&


# Differential missingness filter

# Count number of variants in vcf file
n_var=$(zgrep -v "#" SNP_PASS_QC_0.9_0.8_u_rs2.vcf.gz | wc -l) &&

# Calculate p-value threshold (Bonferrroni alfa 0.05)
p_value=$(echo "scale=10; 0.05/$n_var" | bc) &&

# plink test-missing
plink --vcf SNP_PASS_QC_0.9_0.8_u_rs2.vcf.gz --test-missing --pheno PanDisease_samples.pheno --recode --double-id --allow-extra-chr --allow-no-sex --out SNP_PASS_QC_0.9_0.8 &&

# Find variants with significantly different missingness between cases and controls (Bonferrroni alfa 0.05)
awk -v var="$p_value" '($5 + 0) < var { if (NR!=1) print $2 }' SNP_PASS_QC_0.9_0.8.missing > SNP_PASS_QC_0.9_0.8_sign_dm.txt &&

# Exclude variants in list using SNP id (sign in DM case/control)
vcftools --gzvcf SNP_PASS_QC_0.9_0.8_u_rs2.vcf.gz --exclude SNP_PASS_QC_0.9_0.8_sign_dm.txt --recode --recode-INFO-all --out PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM &&

# Zip vcf file
mv PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM.recode.vcf PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM.vcf &&
bgzip PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM.vcf &&

# Index vcf file
tabix -p vcf PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM.vcf.gz


# Let me know it is done
echo "Job is done"

date
