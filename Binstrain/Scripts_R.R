###########################################################
############# Set the work directory ######################
###########################################################


wd = getwd()
setwd(dir=wd)


###########################################################
################  Read Data ###############################
###########################################################
##### Be careful the format for the data!!! ###############

SNP.Pattern = 
  read.delim("Staph_SNP_Pattern.txt",header = F)

Sample = 
  read.delim("$sample/$sample_BWA/tmp/bwa/Counts_$sample_BWA_sorted_mpileup.txt",header = T)

RefStrains = 
    read.delim("ReferenceStrains.txt",header = T)





##########################################################
################ Load Functions ##########################
##########################################################

source('Functions_BinStrain.R')
ix0 = which(Sample$Coverage==0)
if(length(ix0) > 0)
{
    Sample = Sample[-ix0,]
}


Sample2 = Preprocessing(Sample, SNP.Pattern)
Result2 = Single_est(Sample, SNP.Pattern)
betas = Whole_est(Result2,Sample2)
betas[betas<1e-10] = 0
betas = betas/sum(betas,na.rm= T)
res = cbind(Reference_Strain = RefStrains, Beta_Estimate = betas)
write.table(res, file="Beta_Estimates/$sample_beta.txt",col.names= T,row.names = F,quote= F)



