#!/bin/bash

while read line; do

##Set variables###
comparison=$line
in_Folder='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/060_SKATO/PrepareFiles'
WorkFolder='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/063_SKATO_CADDWeighted_CaddWithLinearConversion'
ScriptsFolder='/proj/sens2017142/nobackup/matteob/SCRIPTS/02_PAN_DISEASE_PROJECT/41_SKATO_Scripts'
Output_Error='/proj/sens2017142/nobackup/matteob/02_PAN_DISEASE_PROJECT/060_SKATO/Output_Error'

echo "#!/bin/bash -l" > PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "#SBATCH -A sens2017142" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "#SBATCH -p core -n 1" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "#SBATCH -J SKATO_PanDisease_${comparison}" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "#SBATCH -t 0-10:00:00" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "#SBATCH -o $Output_Error/SkatO_CaddWeighted_CaddWithLinearConversion_221208_PanDisease_${comparison}.output" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "#SBATCH -e $Output_Error/SkatO_CaddWeighted_CaddWithLinearConversion_221208_PanDisease_${comparison}.error" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "#SBATCH --mail-user matteo.bianchi@imbim.uu.se" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "#SBATCH --mail-type=ALL" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

echo "module load bioinfo-tools" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "module load bcftools" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "module load vcftools" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "module load plink/1.90b4.9" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "module load R_packages/4.0.0" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "module load BEDTools" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "module load snpEff/4.3t" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "module load tabix" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

echo "comparison=${comparison}" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "in_Folder=$in_Folder" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "WorkFolder=$WorkFolder" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "ScriptsFolder=$ScriptsFolder" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "Output_Error=$Output_Error" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

###THE FAM HAS TO HAVE THE PHENOTYPE#####PRODUCE IT BEFORE######
###Copy/Generate all the neededfiles to the Pop folder: bed, bim, fam_WithPheno, SetID, Cov####

echo "cd $WorkFolder" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "mkdir ${comparison}/" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "cd ${comparison}/" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

####Create SSD and Info File####
echo "touch PanDisease_${comparison}.SSD" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "touch PanDisease_${comparison}.Info" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

######ACTUAL SKAT-O#############

echo "R --vanilla <<EOF" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

echo "library(SKAT)" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "library(dplyr)" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

###THE FAM HAS TO HAVE THE PHENOTYPE#####PRODUCE IT######
#######IDs MUST HAVE THE SAME ORDER AS IN THE BED,BIM,FAM FILES PRODUCED!!!!!######

echo "data <- read.csv(\"$in_Folder/PanDisease_ALLSubPhenotypesMatrix_CorrectOrderBasedOnFam.txt\", head=T, sep=\"\t\")" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "data_${comparison} <- data %>% select(\"FID\",\"IID\",\"${comparison}\")" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "data_${comparison}\\\$V3 <- 0" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "data_${comparison}\\\$V4 <- 0" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "data_${comparison}\\\$V5 <- -9" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "data_${comparison}_ok <- data_${comparison} %>% select(\"FID\",\"IID\",\"V3\",\"V4\",\"V5\",\"${comparison}\")" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "write.table(data_${comparison}_ok, \"$in_Folder/PanDisease_${comparison}.fam\", row.names=F, quote=F, col.names=F)" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

echo "File.Bed <- \"$in_Folder/PanDisease.bed\"" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "File.Bim <- \"$in_Folder/PanDisease.bim\"" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "File.Fam <- \"$in_Folder/PanDisease_${comparison}.fam\"" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "File.SetID <- \"$in_Folder/PanDisease_SetID_4SKATO\"" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts; #######TO CHANGE####

echo "Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, \"PanDisease_${comparison}.SSD\", \"PanDisease_${comparison}.Info\")" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

echo "File.SSD <- \"PanDisease_${comparison}.SSD\"" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "File.Info <- \"PanDisease_${comparison}.Info\"" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

echo "File.Cov <- \"$in_Folder/PanDisease_COVARIATES_File\"" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

echo "FAM_COV=Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary=FALSE)" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

echo "y = FAM_COV\\\$Phenotype - 1" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "y[y < 0] <- NA" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

#Generate covariate vectors (SEX and First 4 PCAs)
echo "X1 = FAM_COV\\\$SEX" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "X2 = FAM_COV\\\$PC1" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "X3 = FAM_COV\\\$PC2" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "X4 = FAM_COV\\\$PC3" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "X5 = FAM_COV\\\$PC4" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

echo "SSD.INFO <- Open_SSD(File.SSD, File.Info)" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

#SSD.INFO\$nSample #Number of samples
#SSD.INFO\$nSets   #Number of sets, transcripts (or genes)
#SSD.INFO\$nSNPs   #Number of variants thet is in the set

##########
#FAM=Read_Plink_FAM(File.Fam, Is.binary=TRUE) #Loading the FAM file
#y=FAM$PHENO #Saving the AS phenotype in y
#########

echo "skat.obj = SKAT_Null_Model(y ~ 1 + X1 + X2 + X3 + X4 + X5, out_type=\"D\")" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts; #Making a null model for dichotomous/binary, small sample adjustment (inflate data???)

echo "cadd_weights <- Read_SNP_WeightFile(\"$in_Folder/PanDisease_CaddWeights_WithLinearConversion_4SKAT\")" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

echo "skat.out <- SKATBinary.SSD.All(SSD.INFO, skat.obj, method=\"SKATO\", obj.SNPWeight=cadd_weights)" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "save(skat.out,file=\"ResultsObject_SKATO_COV_CaddWeighted_CaddWithLinearConversion_PanDisease_${comparison}\")" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

echo "write.table(skat.out\\\$results,file=\"Results_SKATO_COV_CaddWeighted_CaddWithLinearConversion_PanDisease_${comparison}.txt\",row.names=F, sep=\"\t\", quote=F)" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

echo "Close_SSD()" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

echo "warnings()" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

echo "EOF" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

echo "R --vanilla <<EOF" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "data <- read.csv(\"Results_SKATO_COV_CaddWeighted_CaddWithLinearConversion_PanDisease_${comparison}.txt\", sep=\"\t\", head=T)" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "data_ok <- subset(data, data\\\$N.Marker.All > 1)" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "data_ok\\\$Bonferroni <- p.adjust(data_ok\\\$P.value, method = \"bonferroni\")" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "data_ok\\\$FDR <- p.adjust(data_ok\\\$P.value, method = \"fdr\")" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "write.table(data_ok,file=\"Results_SKATO_COV_CaddWeighted_CaddWithLinearConversion_PanDisease_${comparison}_Complete.txt\",row.names=F, sep=\"\t\", quote=F)" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;
echo "EOF" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

echo "sort -gk2,2 Results_SKATO_COV_CaddWeighted_CaddWithLinearConversion_PanDisease_${comparison}_Complete.txt > Results_SKATO_COV_CaddWeighted_CaddWithLinearConversion_PanDisease_${comparison}_Complete_sorted.txt" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

echo "rm Results_SKATO_COV_CaddWeighted_CaddWithLinearConversion_PanDisease_${comparison}_Complete" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

#Let me know it is done
echo "echo "Job is done"" >> PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

#sbatch the source file
sbatch PanDisease_SKATO_CaddWeighted_CaddWithLinearConversion.scripts;

done <List_Comparisons_4SKATO_test

