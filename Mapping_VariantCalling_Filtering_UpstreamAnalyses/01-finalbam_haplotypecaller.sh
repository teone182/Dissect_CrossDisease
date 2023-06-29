#!/bin/bash

# This script is used with a list "finalbam.list" 

while read line
do

######## Variables to change #########
pos=$(echo $line | awk 'END{print index($0,"pool")}')
pos=$(expr $pos - 1)
# At what position does the "pool" start?
name=${line:$pos}
# Extract the sample name from the line
location=/proj/b2015143/nobackup/samples/fastq
# Path to the samples in the finalbam.list	
ref=/proj/b2013060/ReferenceTools/homo_sapiens_bwa.0.7.12.fa
# The reference you want to use (indexed with script: /proj/b2013060/Script/Index_ref.sh)
known=/proj/b2013060/ReferenceTools/sorted_common_dbsnp147.bed
# Known sites for the recalibration with GATK (row 111)
bait=/proj/b2013060/ReferenceTools/bait_hg19.bed
# The bait region for HS metrics
target=/proj/b2013060/ReferenceTools/target_hg19.bed	
# The target region for HS metrics


echo "#!/bin/bash -l" > Finalbam_sam.sc;
echo "#SBATCH -A b2015143" >> Finalbam_sam.sc;
echo "#SBATCH -p core -n 2" >> Finalbam_sam.sc;
echo "#SBATCH -J $name-bam" >> Finalbam_sam.sc;
echo "#SBATCH -t 48:00:00" >> Finalbam_sam.sc;
echo "#SBATCH -o /proj/b2015143/nobackup/samples/output/$name-bam.output" >> Finalbam_sam.sc;
echo "#SBATCH -e /proj/b2015143/nobackup/samples/output/$name-bam.error" >> Finalbam_sam.sc;
echo "#SBATCH --mail-user lina.hultin-rosenberg@imbim.uu.se" >> Finalbam_sam.sc;
echo "#SBATCH --mail-type=ALL" >> Finalbam_sam.sc;


# Load modules
echo "module load bioinfo-tools" >> Finalbam_sam.sc;
echo "module load samtools" >> Finalbam_sam.sc;
echo "module load bwa/0.7.12" >> Finalbam_sam.sc;
echo "module load tabix" >> Finalbam_sam.sc;

# Variables, same as above #
echo "ref=$ref" >> Finalbam_sam.sc;
echo "location=$location" >> Finalbam_sam.sc;
echo "line=$line" >> Finalbam_sam.sc;
echo "name=$name" >> Finalbam_sam.sc;
echo "known=$known" >> Finalbam_sam.sc;
echo "bait=$bait" >> Finalbam_sam.sc;
echo "target=$target" >> Finalbam_sam.sc;


echo "cd /proj/b2015143/nobackup/samples/running/" >> Finalbam_sam.sc;


# Go from fastQ to BAM #

# Create SAM files in $SNIC_TMP folder to save space, automatically deleted after job finishes

# Align the read 1 data for barcode 1 and 2 [Note -t is the number of threads]
echo "bwa mem -t 6 -M -R '@RG\tID:$name\tSM:bar' $ref $location/$line/*_R1_* $location/$line/*_R2_* > \$SNIC_TMP/$name.sam" >> Finalbam_sam.sc;

  
# Use PICARD, start with the SAM file and ask it to simultaneously add read groups, sort the file, and spit it out as BAM. 
echo "java -Xmx15g -jar /sw/apps/bioinfo/picard/1.92/milou/AddOrReplaceReadGroups.jar \
 INPUT=\$SNIC_TMP/$name.sam OUTPUT=$name.bam SORT_ORDER=coordinate RGID=$name-id RGLB=$name-lib \
 RGPL=ILLUMINA RGPU=$name-01 RGSM=$name VALIDATION_STRINGENCY=LENIENT" >> Finalbam_sam.sc;


# Index the Picard derived BAM
echo "java -Xmx15g -jar /sw/apps/bioinfo/picard/1.92/milou/BuildBamIndex.jar INPUT=$name.bam \
 VALIDATION_STRINGENCY=LENIENT" >> Finalbam_sam.sc;


# Process reads with GATK #

# Step 1 - realign locally around potential indels. First, identify possible sites to realign.
echo "java -Xmx15g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -I $name.bam -R $ref \
 -T RealignerTargetCreator -o $name.intervals" >> Finalbam_sam.sc;

# Step 2 - Feed the intervals file back into GATK with a different argument to actually do the realignments
echo "java -Xmx15g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -I $name.bam -R $ref \
 -T IndelRealigner -o $name.realigned.bam -targetIntervals $name.intervals" >> Finalbam_sam.sc;

echo "rm $name.ba* $name.intervals" >> Finalbam_sam.sc;
 
# Use Picard to mark duplicate reads
echo "java -Xmx15g -jar /sw/apps/bioinfo/picard/1.92/milou/MarkDuplicates.jar INPUT=$name.realigned.bam \
 OUTPUT=$name.realigned.marked.bam METRICS_FILE=$name.realigned.marked.metrics VALIDATION_STRINGENCY=LENIENT" >> Finalbam_sam.sc;
 
echo "rm $name.realigned.ba*" >> Finalbam_sam.sc;

# Reindexing the picard file
echo "java -Xmx15g -jar /sw/apps/bioinfo/picard/1.92/milou/BuildBamIndex.jar INPUT=$name.realigned.marked.bam \
 VALIDATION_STRINGENCY=LENIENT" >> Finalbam_sam.sc;


# Perform quality recalibration with GATK #

# Step 1 - Use a list of known sites
echo "java -Xmx15g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T BaseRecalibrator -I $name.realigned.marked.bam \
 -R $ref -knownSites $known -o $name.recal_data.table -rf BadCigar" >> Finalbam_sam.sc; 

# Step 2 - Base recalibration
echo "java -Xmx15g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T BaseRecalibrator -I $name.realigned.marked.bam \
 -R $ref -knownSites $known -BQSR $name.recal_data.table -o $name.post_recal_data.table " >> Finalbam_sam.sc;

# Step 3 - Print reads
echo "java -Xmx15g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T PrintReads -I $name.realigned.marked.bam \
 -R $ref -BQSR $name.recal_data.table -o $name.realigned.marked.recalibrated.bam" >> Finalbam_sam.sc;


# Quality test #

# Flagstat
echo "samtools flagstat $name.realigned.marked.recalibrated.bam > ../quality_test/$name.flagstat" >> Finalbam_sam.sc;

# HSmetrics
echo "java -Xmx15g -jar /sw/apps/bioinfo/picard/1.92/milou/CalculateHsMetrics.jar R=$ref BAIT_INTERVALS=$bait \
 TARGET_INTERVALS=$target INPUT=$name.realigned.marked.recalibrated.bam OUTPUT=../quality_test/$name.hybrid_selection_metrics \
 VALIDATION_STRINGENCY=LENIENT" >> Finalbam_sam.sc;


# Delete files #
# Delete files to free up space
echo "rm $name.recal_data.table $name.post_recal_data.table $name.realigned.marked.ba* $name.realigned.marked.metrics" >> Finalbam_sam.sc;


# Move and rename BAM files #
echo "cd /proj/b2015143/nobackup/samples/finalbam/" >> Finalbam_sam.sc;
echo "mv /proj/b2015143/nobackup/samples/running/$name.realigned.marked.recalibrated.bam $name.final.bam" >> Finalbam_sam.sc;
echo "mv /proj/b2015143/nobackup/samples/running/$name.realigned.marked.recalibrated.bai $name.final.bai" >> Finalbam_sam.sc;


# Go to haplotypecaller folder
echo "cd /proj/b2015143/nobackup/samples/haplotypecaller/" >> Finalbam_sam.sc;


# Do genotype calling in each individual #

echo "java -Xmx15g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T HaplotypeCaller -R $ref \
 -I /proj/b2015143/nobackup/samples/finalbam/$name.final.bam --emitRefConfidence GVCF --variant_index_type LINEAR \
 --variant_index_parameter 128000 -o $name.g.vcf.gz -rf BadCigar" >> Finalbam_sam.sc;

echo "java -Xmx15g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -R $ref \
  -T ValidateVariants  --validationTypeToExclude ALL --variant $name.g.vcf.gz" >> Finalbam_sam.sc;

# Zip g.vcf file
#echo "bgzip $name.g.vcf" >> Finalbam_sam.sc;


# Let me know it is done
echo "echo "Job is done"" >> Finalbam_sam.sc;

echo "date" >> Finalbam_sam.sc;

sbatch Finalbam_sam.sc;

done <finalbam.list
