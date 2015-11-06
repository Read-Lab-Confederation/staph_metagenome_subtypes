##########################################################################################
####### Set the working directory and make the output directory for beta estimates #######
##########################################################################################

CWD="$(pwd)"
mkdir $CWD/Beta_Estimates/

##########################################################################################
############# While loop for analyzing each sample one after the other ###################
##########################################################################################


while read sample
do

#############################################################################################
########Directly downloading the HMP fastq files from the HMP ftp site ######################
#############################################################################################


wget -m ftp://public-ftp.hmpdacc.org/Illumina/PHASEII/tongue_dorsum/${sample}.tar.bz2
mv $CWD/public-ftp.hmpdacc.org/Illumina/PHASEII/tongue_dorsum/${sample}.tar.bz2 .
tar jxf ${sample}.tar.bz2
cd ${sample}

#################################################################################################
############# Concatenating all the fastq files into one single fastq file ######################
#################################################################################################



cat ${sample}.denovo_duplicates_marked.trimmed.1.fastq ${sample}.denovo_duplicates_marked.trimmed.2.fastq ${sample}.denovo_duplicates_marked.trimmed.singleton.fastq >> All_${sample}.fastq


#################################################################################################
############# removing unwanted fastq files #####################################################
#################################################################################################

rm ${sample}.denovo_duplicates_marked.trimmed.1.fastq
rm ${sample}.denovo_duplicates_marked.trimmed.2.fastq
rm ${sample}.denovo_duplicates_marked.trimmed.singleton.fastq


#################################################################################################
############# making folders for keeping the bwa alignment files ################################
#################################################################################################


mkdir -p ${sample}_BWA/tmp/bwa
mkdir -p ${sample}_BWA/tmp/bwa/index


#################################################################################################
############# bring up the reference genome for reference mapping using bwa #####################
#################################################################################################

cp $CWD/Staph_ASR.fasta ${sample}_BWA/tmp/bwa/index/ref.fa

#################################################################################################
############# performing the BWA alignment ######################################################
#################################################################################################


/usr/bin/bwa index ${sample}_BWA/tmp/bwa/index/ref.fa
java -jar /usr/local/bin/CreateSequenceDictionary.jar R=${sample}_BWA/tmp/bwa/index/ref.fa O=${sample}_BWA/tmp/bwa/index/ref.dict
/usr/bin/samtools faidx ${sample}_BWA/tmp/bwa/index/ref.fa
/usr/bin/bwa aln -e 10 -f ${sample}_BWA/tmp/bwa/bwa.sai ${sample}_BWA/tmp/bwa/index/ref.fa All_${sample}.fastq
/usr/bin/bwa samse -f ${sample}_BWA/tmp/bwa/bwa.sam ${sample}_BWA/tmp/bwa/index/ref.fa ${sample}_BWA/tmp/bwa/bwa.sai All_${sample}.fastq

#################################################################################################
############# calling Samstools to generate the mpileup file ####################################
#################################################################################################


cd $CWD/${sample}/${sample}_BWA/tmp/bwa/
samtools view -b -S -o ${sample}_BWA.bam bwa.sam
samtools sort ${sample}_BWA.bam ${sample}_BWA_sorted
samtools mpileup -f $CWD/Staph_ASR.fasta ${sample}_BWA_sorted.bam > ${sample}_BWA_sorted_mpileup.txt


#################################################################################################
############# Formating mpileup file to generate binstrain input ################################
#################################################################################################

awk 'BEGIN{print "Value\tRefNt\tCoverage\t.\t,\tA/a\tG/g\tC/c\tT/t"}    { s=tolower($5);l=split(s,a,"");

 for (i=1;i<=l;i++) b[a[i]]++;

 print $2,$3,$4,b["."]+0,b[","]+0,b["a"]+0,b["g"]+0,b["c"]+0,b["t"]+0;

 delete a;delete b}' OFS="\t" ${sample}_BWA_sorted_mpileup.txt > Counts_${sample}_BWA_sorted_mpileup.txt
cd ..
cd ..
cd ..
cd ..

#################################################################################################
############# Calling in and running the R script of binstrain ##################################
#################################################################################################

sed 's/\$sample/'"$sample"'/g' Scripts_R.R > Scripts_Staph_R_NEW.R
R-3.0.1 CMD BATCH Scripts_Staph_R_NEW.R
rm Scripts_Staph_R_NEW.R

#################################################################################################
############# List of Samples that need to be analyzed ##########################################
#################################################################################################

done < $CWD/list_samples.txt
