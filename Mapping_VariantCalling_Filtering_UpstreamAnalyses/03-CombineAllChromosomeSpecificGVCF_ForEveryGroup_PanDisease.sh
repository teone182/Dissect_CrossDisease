#!/bin/bash

while read line
do

#### Variables to change ####
name=$line	
ref=/proj/sens2017142/sens2017107/DISSECT/Reference/homo_sapiens_bwa.0.7.12.fa

echo "#!/bin/bash -l" > MergeChromosomes.sc;
echo "#SBATCH -A sens2017142" >> MergeChromosomes.sc;
echo "#SBATCH -p core -n 4" >> MergeChromosomes.sc;
echo "#SBATCH -J $name-MergeChromosomes" >> MergeChromosomes.sc;
echo "#SBATCH -t 6-00:00:00" >> MergeChromosomes.sc;
echo "#SBATCH -o /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/output/$name-MergeChromosomes-210101.output" >> MergeChromosomes.sc;
echo "#SBATCH -e /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/output/$name-MergeChromosomes-210101.error" >> MergeChromosomes.sc;
echo "#SBATCH --mail-user matteo.bianchi@imbim.uu.se" >> MergeChromosomes.sc;
echo "#SBATCH --mail-type=ALL" >> MergeChromosomes.sc;

# Load modules
echo "module load bioinfo-tools" >> MergeChromosomes.sc;
echo "module load vcftools" >> MergeChromosomes.sc;
echo "module load tabix" >> MergeChromosomes.sc; 

# Variables, same as above #
echo "name=$name"  >> MergeChromosomes.sc;
echo "ref=$ref"  >> MergeChromosomes.sc;


### Make one file to use in genotyping ###

#Which directory are you working in?
echo "cd /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/002_Grouped_gVCF_PerChromosome"  >> MergeChromosomes.sc;

# Combine chromosome files 
echo "java -Xmx23g -cp /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
  -R $ref  -V chr1-$name.vcf.gz -V chr2-$name.vcf.gz -V chr3-$name.vcf.gz -V chr4-$name.vcf.gz -V chr5-$name.vcf.gz \
  -V chr6-$name.vcf.gz -V chr7-$name.vcf.gz -V chr8-$name.vcf.gz -V chr9-$name.vcf.gz -V chr10-$name.vcf.gz \
  -V chr11-$name.vcf.gz -V chr12-$name.vcf.gz -V chr13-$name.vcf.gz -V chr14-$name.vcf.gz -V chr15-$name.vcf.gz \
  -V chr16-$name.vcf.gz -V chr17-$name.vcf.gz -V chr18-$name.vcf.gz -V chr19-$name.vcf.gz -V chr20-$name.vcf.gz \
  -V chr21-$name.vcf.gz -V chr22-$name.vcf.gz -V chrX-$name.vcf.gz -V chrY-$name.vcf.gz -V chrM-$name.vcf.gz \
  -out /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/003_Grouped_gVCF_ChromosomesCombined/$name.vcf.gz -assumeSorted"  >> MergeChromosomes.sc; 

  
echo "java -Xmx23g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -R $ref \
  -T ValidateVariants --validationTypeToExclude ALL --variant /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/003_Grouped_gVCF_ChromosomesCombined/$name.vcf.gz"  >> MergeChromosomes.sc;


# Let me know it is done
echo "echo "Job is done""  >> MergeChromosomes.sc;

echo "date"  >> MergeChromosomes.sc;

sbatch MergeChromosomes.sc;

done<gVCF_group.list
