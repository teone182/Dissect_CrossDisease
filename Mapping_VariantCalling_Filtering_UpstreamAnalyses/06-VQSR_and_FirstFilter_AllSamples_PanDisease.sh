#!/bin/bash -l
#SBATCH -A sens2017142
#SBATCH -p core -n 8
#SBATCH -J VQSR_FirstFilter_PanDisease
#SBATCH -t 7-00:00:00
#SBATCH -o /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/output/PanDisease_VQSR_and_FirstFilter-210114.output
#SBATCH -e /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/output/PanDisease_VQSR_and_FirstFilter-210114.error
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

###########STEP 1###############

### VQSR ###
java -Xmx23g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T VariantRecalibrator -R $ref \
 --input $InDir/PanDisease_AllSamples.vcf.gz \
 -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
 -resource:omni,known=false,training=true,truth=true,prior=12.0 $omni \
 -resource:1000G,known=false,training=true,truth=false,prior=10.0 $snphc \
 -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
 -an QD -an FS -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 \
 -tranche 99.0 -tranche 90.0 -recalFile recalibrate_SNP.recal \
 -tranchesFile recalibrate_SNP.tranches -rscriptFile recalibrate_SNP_plots.R --maxGaussians 4 

# Use recalibration model on snp call set
java -Xmx23g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T ApplyRecalibration -R $ref \
 --input $InDir/PanDisease_AllSamples.vcf.gz -mode SNP \
 --ts_filter_level 99.0 -recalFile recalibrate_SNP.recal -tranchesFile recalibrate_SNP.tranches -o PanDisease_SNP_raw_indels.vcf.gz 


# Get only positions that PASSED VQSR (only SNPs have been recalibrated)
vcftools --gzvcf PanDisease_SNP_raw_indels.vcf.gz --keep-filtered PASS --remove-filtered-geno-all --remove-filtered-all --out PanDisease_SNP_PASS --recode --recode-INFO-all &&

# Zip vcf
mv PanDisease_SNP_PASS.recode.vcf PanDisease_SNP_PASS.vcf &&
bgzip PanDisease_SNP_PASS.vcf &&

# Index vcf file
tabix -p vcf PanDisease_SNP_PASS.vcf.gz &&



# Min DP and Min GQ filter (does not filter out positions, but exchanges genotypes for ./. if criteria not fulfilled)
vcftools --gzvcf PanDisease_SNP_PASS.vcf.gz --minDP 8 --minGQ 20 --recode --recode-INFO-all --out PanDisease_SNP_PASS_DP_GQ &&

# Zip vcf
mv PanDisease_SNP_PASS_DP_GQ.recode.vcf PanDisease_SNP_PASS_DP_GQ.vcf &&
bgzip PanDisease_SNP_PASS_DP_GQ.vcf &&

# Index vcf file
tabix -p vcf PanDisease_SNP_PASS_DP_GQ.vcf.gz &&


# Get only polymorphic positions
vcftools --gzvcf PanDisease_SNP_PASS_DP_GQ.vcf.gz --maf 0.00001 --max-maf 0.99992 --out PanDisease_SNP_PASS_DP_GQ_P --recode --recode-INFO-all &&

# Zip vcf
mv PanDisease_SNP_PASS_DP_GQ_P.recode.vcf PanDisease_SNP_PASS_DP_GQ_P.vcf &&
bgzip PanDisease_SNP_PASS_DP_GQ_P.vcf &&

# Index vcf file
tabix -p vcf PanDisease_SNP_PASS_DP_GQ_P.vcf.gz &&



# Get the AB annotation for your vcf
java -Xmx23g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T VariantAnnotator -R $ref -V PanDisease_SNP_PASS_DP_GQ_P.vcf.gz -o PanDisease_SNP_PASS_DP_GQ_P_ab.vcf.gz -A AlleleBalance -A AlleleBalanceBySample &&

# Filter on AB (0.2 < AB < 0.8)
zgrep -v ^# PanDisease_SNP_PASS_DP_GQ_P_ab.vcf.gz | grep "ABHet=" | awk '{print $1, $2, $8}' | awk -F ';' '{print $1}' | awk -F '=' '{print $1,$2}' | awk '$4<0.2 {print $1,$2,$2,$3,$4}' > PanDisease_SNP_PASS_AB2.bed &&
zgrep -v ^# PanDisease_SNP_PASS_DP_GQ_P_ab.vcf.gz | grep "ABHet=" | awk '{print $1, $2, $8}' | awk -F ';' '{print $1}' | awk -F '=' '{print $1,$2}' | awk '$4>0.8 {print $1,$2,$2,$3,$4}' > PanDisease_SNP_PASS_AB8.bed &&

vcftools --gzvcf PanDisease_SNP_PASS_DP_GQ_P_ab.vcf.gz --exclude-bed PanDisease_SNP_PASS_AB2.bed --recode --recode-INFO-all --out PanDisease_SNP_PASS_DP_GQ_P_AB2 &&
vcftools --vcf PanDisease_SNP_PASS_DP_GQ_P_AB2.recode.vcf --exclude-bed PanDisease_SNP_PASS_AB8.bed  --recode --recode-INFO-all --out PanDisease_SNP_PASS_DP_GQ_P_AB &&

# Zip vcf
mv PanDisease_SNP_PASS_DP_GQ_P_AB.recode.vcf PanDisease_SNP_PASS_DP_GQ_P_AB.vcf &&
bgzip PanDisease_SNP_PASS_DP_GQ_P_AB.vcf &&

# Index vcf file
tabix -p vcf PanDisease_SNP_PASS_DP_GQ_P_AB.vcf.gz &&

# Let me know it is done
echo "Job is done"

date
