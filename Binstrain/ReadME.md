#Running BinStrain for SNP based genotyping of Staphylococcus aureus subtypes from metagenome shotgun data

Version: 1.2

Date: 2016-08-10

Author: Sandeep J. Joseph and Ben Li

Maintainer: Sandeep J. Joseph <sandeepjoseph@emory.edu> and Ben Li <ben.li@emory.edu>

binstrain uses a binomial mixture model to describe the observed alternative allele (SNP) derived from prior comparative genomic analysis to estimate the proportion of 40 S. aureus subtypes in metagenome samples. The 2 main input files are the SNP pattern file that contains the position information on subtype-specific SNPs as well as SNPs shared by other subtypes and coverage counts at these SNP sites, obtained by reference mapping of the reads to a reference genome. A well-established two step method is used to estimate the proportion of serovar strain-specific SNPs present in pure or mixed infection of genotypes. First step is a direct estimate by using the sparsity of the design matrix. Quadratic optimization method is involved in the second step.  Originally developed for direct sequencing of C. trachomatis from clinical samples.

License: GPL-2

###Dependencies

bwa
samtools
R-3-0.1
awk

###Usage

You need to have the following files in the same folder where you plan to run the analysis:

*[Functions_BinStrain.R](https://github.com/Read-Lab-Confederation/staph_metagenome_subtypes/blob/master/Binstrain/Functions_BinStrain.R)

*[ReferenceStrains.txt](https://github.com/Read-Lab-Confederation/staph_metagenome_subtypes/blob/master/Binstrain/ReferenceStrains.txt)

*[Scripts_R.R](https://github.com/Read-Lab-Confederation/staph_metagenome_subtypes/blob/master/Binstrain/Scripts_R.R)

*[Staph_ASR.fasta](https://github.com/Read-Lab-Confederation/staph_metagenome_subtypes/blob/master/Binstrain/Staph_ASR.fasta)

*[Staph_SNP_Pattern.txt](https://github.com/Read-Lab-Confederation/staph_metagenome_subtypes/blob/master/Binstrain/Staph_SNP_Pattern.txt_

*[list_samples.txt](https://github.com/Read-Lab-Confederation/staph_metagenome_subtypes/blob/master/Binstrain/list_samples.txt)

*[binStrain.sh](https://github.com/Read-Lab-Confederation/staph_metagenome_subtypes/blob/master/Binstrain/binStrain.sh)
        
Details regarding each of the above scripts are commented out in the [binStrain.sh shell script](https://github.com/Read-Lab-Confederation/staph_metagenome_subtypes/blob/master/Binstrain/binStrain.sh)

###Example Usage

If you want to download the raw sequence data generated from Human Microbiome Project, you need to get the SRS id from [HMP site](http://hmpdacc.org/HMASM/) and list the SRS ids in the file list_samples.txt.

Then run the following command in the terminal

    ./binStrain.sh

The example list\_samples.txt file provided contains 3 SRS ids - SRS016969, SRS020571, SRS020628. The shell script binStrain.sh will download the corresponding SRS raw reads from the ftp site of HMP and process the data and estimate the beta values. The final beta estimates will be generated into a folder named Beta_Estimates.

###How to generate the SNP pattern Matrix.

1. Generate a whole genome alignment from the selected reference genomes using [ProgreesiveMauve](http://darlinglab.org/mauve/user-guide/introduction.html)
2. Use the "export SNPs" function in MAUVE to export the SNPs from the alignment.
3. Trim the positions where missing or gaps (-) are present at any of the SNP sites of the reference genomes in order to make sure we are getting the core SNPs.
4. Use the awk script Format\_SNP\_Pattern.awk to format the output file from step 3 into the input SNPpattern file for binstrain. Usage: 

        awk Format_SNP_Pattern.awk inputfile > Final_SNP_patternFile.txt

