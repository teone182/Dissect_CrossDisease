#!/bin/bash

while read line
do

######## Variables to change #########
group=group01
name=$line
ref=/proj/sens2017142/sens2017107/DISSECT/Reference/homo_sapiens_bwa.0.7.12.fa

echo "#!/bin/bash -l" > combineGVCF_Group01.sc;
echo "#SBATCH -A sens2017142" >> combineGVCF_Group01.sc;
echo "#SBATCH -p core -n 1" >> combineGVCF_Group01.sc;
echo "#SBATCH -J $name-$group-combineGVCF-201230" >> combineGVCF_Group01.sc;
echo "#SBATCH -t 30:00:00" >> combineGVCF_Group01.sc;
echo "#SBATCH -o /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/output/$name-$group-combineGVCF-201230.output" >> combineGVCF_Group01.sc;
echo "#SBATCH -e /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/output/$name-$group-combineGVCF-201230.error" >> combineGVCF_Group01.sc;

# Load modules
echo "module load bioinfo-tools" >> combineGVCF_Group01.sc;
echo "module load tabix" >> combineGVCF_Group01.sc;

# Variables, same as above
echo "ref=$ref" >> combineGVCF_Group01.sc;
echo "group=$group" >> combineGVCF_Group01.sc;
echo "name=$name" >> combineGVCF_Group01.sc;
echo "line=$line" >> combineGVCF_Group01.sc;

# Where should the group.vcf files be created?
echo "cd /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/002_Grouped_gVCF_PerChromosome" >> combineGVCF_Group01.sc;

# Combine GVCFs (list all samples that should be combined in this group)

echo "java -Xmx7g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T CombineGVCFs -R $ref \
-V /proj/sens2017142/sens2017107/DISSECT/SLE/individual_vcf/sample_0001.g.vcf.gz \
-V /proj/sens2017142/sens2017107/DISSECT/Scand_Myositis/individual_vcf/sample_0002.g.vcf.gz \
-V /proj/sens2017142/sens2017107/DISSECT/Sjogren/individual_vcf/sample_nnnn.g.vcf.gz \
-V /proj/sens2017142/sens2017107/DISSECT/Controls/individual_vcf/sample_0100.g.vcf.gz \
-o $name-$group.vcf.gz -L $name" >> combineGVCF_Group01.sc;

# Validate variants
echo "java -Xmx7g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -R $ref \
   -T ValidateVariants  --validationTypeToExclude ALL --variant $name-$group.vcf.gz" >> combineGVCF_Group01.sc;


# Let me know it is done
echo "echo "Job is done"" >> combineGVCF_Group01.sc;

echo "date" >> combineGVCF_Group01.sc;

# sbatch the source file
sbatch combineGVCF_Group01.sc;

# Must be the same as the file input when calling the script
done<chromosome.list

#######This for all the 41 groups (~1400 samples in total)#####
