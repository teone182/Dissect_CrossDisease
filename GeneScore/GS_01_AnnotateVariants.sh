#!/bin/bash

while read line; do

#Set variables
chromosome=$line
InputFilesFolder='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/007_AllSamples_Cleaned_OverlappingWith2021DatabaseInfo'
WorkFolder='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/040_Cleaned_OverlappingWith2021DatabaseInfo_+-2Kb_GeneRiskScore'
Annovar='/proj/sens2017142/Annovar'
Scripts='/proj/sens2017142/nobackup/matteob/SCRIPTS/02_PAN_DISEASE_PROJECT/40_GeneRiskScore_Scripts'

echo "#!/bin/bash -l" > 01_PanDisease_GeneRiskScore.sc;
echo "#SBATCH -A sens2017142" >> 01_PanDisease_GeneRiskScore.sc;
echo "#SBATCH -p core -n 8" >> 01_PanDisease_GeneRiskScore.sc;  
echo "#SBATCH -J GeneRiskScore_PanDisease_${chromosome}" >> 01_PanDisease_GeneRiskScore.sc;
echo "#SBATCH -t 0-24:00:00" >> 01_PanDisease_GeneRiskScore.sc;
echo "#SBATCH -o /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/040_Cleaned_OverlappingWith2021DatabaseInfo_+-2Kb_GeneRiskScore/01_PanDisease_GeneRiskScore_2Kb_${chromosome}.output" >> 01_PanDisease_GeneRiskScore.sc;
echo "#SBATCH -e /proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/040_Cleaned_OverlappingWith2021DatabaseInfo_+-2Kb_GeneRiskScore/01_PanDisease_GeneRiskScore_2Kb_${chromosome}.error" >> 01_PanDisease_GeneRiskScore.sc;
echo "#SBATCH --mail-user matteo.bianchi@imbim.uu.se" >> 01_PanDisease_GeneRiskScore.sc;
echo "#SBATCH --mail-type=ALL" >> 01_PanDisease_GeneRiskScore.sc;

echo "module load bioinfo-tools" >> 01_PanDisease_GeneRiskScore.sc;
echo "module load vcftools" >> 01_PanDisease_GeneRiskScore.sc;
echo "module load bcftools" >> 01_PanDisease_GeneRiskScore.sc;
echo "module load snpEff" >> 01_PanDisease_GeneRiskScore.sc;
echo "module load python/2.7.9" >> 01_PanDisease_GeneRiskScore.sc;
echo "module load R" >> 01_PanDisease_GeneRiskScore.sc;
echo "module load tabix" >> 01_PanDisease_GeneRiskScore.sc;
echo "module load BEDTools" >> 01_PanDisease_GeneRiskScore.sc;

echo "chromosome=$chromosome" >> 01_PanDisease_GeneRiskScore.sc;
echo "InputFilesFolder=$InputFilesFolder" >> 01_PanDisease_GeneRiskScore.sc;
echo "WorkFolder=$WorkFolder" >> 01_PanDisease_GeneRiskScore.sc;
echo "Annovar=$Annovar" >> 01_PanDisease_GeneRiskScore.sc;
echo "Scripts=$Scripts" >> 01_PanDisease_GeneRiskScore.sc;

echo "cd $WorkFolder" >> 01_PanDisease_GeneRiskScore.sc;

#####split vcf files into chromosomes##########
echo "bcftools view -r ${chromosome} -Oz -o PanDisease_${chromosome}.vcf.gz $InputFilesFolder/PAN_DISEASE_SNP_PASS_QC_0.9_0.8_DM_Cleaned_OverlappingWithDB_NoFalsePositives.vcf.gz" >> 01_PanDisease_GeneRiskScore.sc;

####make annovar-ready file####
echo "$Annovar/convert2annovar.pl -format vcf4 $WorkFolder/PanDisease_${chromosome}.vcf.gz -outfile PanDisease_${chromosome}.input -allsample -withfreq -include 2> PanDisease_annovar_${chromosome}.log" >> 01_PanDisease_GeneRiskScore.sc;

####Annotate with Annovar different annotations######
echo "$Annovar/table_annovar.pl PanDisease_${chromosome}.input $Annovar/humandb/ --buildver hg19 -out PanDisease_${chromosome} -remove -protocol refGene,cadd13 -operation g,f --argument '--neargene 2000, --neargene 2000' --thread 40 --maxgenethread 40 -nastring . >> PanDisease_annovar_${chromosome}.log " >> 01_PanDisease_GeneRiskScore.sc;

echo "cut -f 18- PanDisease_${chromosome}.input > a1_${chromosome}" >> 01_PanDisease_GeneRiskScore.sc;
echo "zgrep '^#CHR' $WorkFolder/PanDisease_${chromosome}.vcf.gz | cut -f 10- > b1_${chromosome}" >> 01_PanDisease_GeneRiskScore.sc;
echo "cat b1_${chromosome} a1_${chromosome} > PanDisease_temp_${chromosome}" >> 01_PanDisease_GeneRiskScore.sc;
echo "paste PanDisease_${chromosome}.hg19_multianno.txt PanDisease_temp_${chromosome} > PanDisease_${chromosome}.meta" >> 01_PanDisease_GeneRiskScore.sc;
echo "rm a1_${chromosome} b1_${chromosome} PanDisease_temp_${chromosome}" >> 01_PanDisease_GeneRiskScore.sc;

#Let me know it is done
echo "echo "Job is done"" >> 01_PanDisease_GeneRiskScore.sc;

#sbatch the source file
sbatch 01_PanDisease_GeneRiskScore.sc;

#Must be the same as the file input when calling the script
done <List_chromosomes.list

