#!/bin/bash
#SBATCH -A sens2017142
#SBATCH -p core -n 8
#SBATCH -J PanDisease_AllSamples
#SBATCH -t 0-84:00:00
#SBATCH -o /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/output/PanDisease_AfterJointGeno-210109.output
#SBATCH -e /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/output/PanDisease_AfterJointGeno-210109.error
#SBATCH --mail-user matteo.bianchi@imbim.uu.se
#SBATCH --mail-type=ALL


# Load modules
module load bioinfo-tools
module load vcftools
module load tabix 

name='PanDisease_JointGenotyping'
ref='/proj/sens2017142/sens2017107/DISSECT/Reference/homo_sapiens_bwa.0.7.12.fa'

cd /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/004_JointGenotyping 

# Combine per chromosome VCF files 
java -Xmx23g -cp /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
  -R $ref  -V chr1-$name.vcf.gz -V chr2-$name.vcf.gz -V chr3-$name.vcf.gz -V chr4-$name.vcf.gz -V chr5-$name.vcf.gz \
  -V chr6-$name.vcf.gz -V chr7-$name.vcf.gz -V chr8-$name.vcf.gz -V chr9-$name.vcf.gz -V chr10-$name.vcf.gz \
  -V chr11-$name.vcf.gz -V chr12-$name.vcf.gz -V chr13-$name.vcf.gz -V chr14-$name.vcf.gz -V chr15-$name.vcf.gz \
  -V chr16-$name.vcf.gz -V chr17-$name.vcf.gz -V chr18-$name.vcf.gz -V chr19-$name.vcf.gz -V chr20-$name.vcf.gz \
  -V chr21-$name.vcf.gz -V chr22-$name.vcf.gz -V chrX-$name.vcf.gz -V chrY-$name.vcf.gz -V chrM-$name.vcf.gz \
  -out /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/005_AllSamples_JointGenotyping/PanDisease_AllSamples.vcf.gz -assumeSorted 

# Validate variants
java -Xmx23g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -R $ref \
  -T ValidateVariants --validationTypeToExclude ALL --variant /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/005_AllSamples_JointGenotyping/PanDisease_AllSamples.vcf.gz

# Let me know it is done
echo "Job is done"

date
