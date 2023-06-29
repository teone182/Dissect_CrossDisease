#!/bin/bash

while read line
do

######## Variables to change #########
name=$line
ref=/proj/sens2017142/sens2017107/DISSECT/Reference/homo_sapiens_bwa.0.7.12.fa

echo "#!/bin/bash -l" > JointGenotyping.sc;
echo "#SBATCH -A sens2017142" >> JointGenotyping.sc;
echo "#SBATCH -p node -n 16" >> JointGenotyping.sc;
echo "#SBATCH -J $name-JointGenotyping-210103" >> JointGenotyping.sc;
echo "#SBATCH -t 5-00:00:00" >> JointGenotyping.sc; #Change the time
echo "#SBATCH -o /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/output/$name-JointGenotyping-210103.output" >> JointGenotyping.sc;
echo "#SBATCH -e /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/output/$name-JointGenotyping-210103.error" >> JointGenotyping.sc;

# Load modules
echo "module load bioinfo-tools" >> JointGenotyping.sc;
echo "module load tabix" >> JointGenotyping.sc;

# Variables, same as above #
echo "ref=$ref" >> JointGenotyping.sc;
echo "name=$name" >> JointGenotyping.sc;
echo "line=$line" >> JointGenotyping.sc;


# Where are the group.vcf files?
echo "cd /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/003_Grouped_gVCF_ChromosomesCombined" >> JointGenotyping.sc;

# Genotype (list the group.vcf files for all samples in cohort)
echo "java -Xmx128g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $ref -nt 16 \
-V group01.vcf.gz \
-V group02.vcf.gz \
-V group03.vcf.gz \
-V group04.vcf.gz \
-V group05.vcf.gz \
-V group06.vcf.gz \
-V group07.vcf.gz \
-V group08.vcf.gz \
-V group09.vcf.gz \
-V group10.vcf.gz \
-V group11.vcf.gz \
-V group12.vcf.gz \
-V group13.vcf.gz \
-V group14.vcf.gz \
-V group15.vcf.gz \
-V group16.vcf.gz \
-V group17.vcf.gz \
-V group18.vcf.gz \
-V group19.vcf.gz \
-V group20.vcf.gz \
-V group21.vcf.gz \
-V group22.vcf.gz \
-V group23.vcf.gz \
-V group24.vcf.gz \
-V group25.vcf.gz \
-V group26.vcf.gz \
-V group27.vcf.gz \
-V group28.vcf.gz \
-V group29.vcf.gz \
-V group30.vcf.gz \
-V group31.vcf.gz \
-V group32.vcf.gz \
-V group33.vcf.gz \
-V group34.vcf.gz \
-V group35.vcf.gz \
-V group36.vcf.gz \
-V group37.vcf.gz \
-V group38.vcf.gz \
-V group39.vcf.gz \
-V group40.vcf.gz \
-V group41.vcf.gz \
-o /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/004_JointGenotyping/$name-PanDisease_JointGenotyping.vcf.gz -L $name" >> JointGenotyping.sc;


# Let me know it is done
echo "echo "Job is done"" >> JointGenotyping.sc;

echo "date" >> JointGenotyping.sc;

# sbatch the source file
sbatch JointGenotyping.sc;

# Must be the same as the file input when calling the script
done<chromosome.list


