Preprocessing = 
  function (Counts_Sample, SNP_Pattern, N = 102057)
  {
    IDs = SNP_Pattern$V1
    
    Sample = Counts_Sample[which(Counts_Sample$Value %in% IDs), ]
    
    Sample_num_reads = Sample$Coverage
    
    
    ### if have no "coverage"
    #Sample_num_reads = rowSums(Sample[,4:9])
    
    names(Sample)[6:9]=c("A","G","C","T")
    
    
    
    
    
    Mat_Zs = matrix(0, nrow = dim(Sample)[1], ncol = 40)
    
    if (dim(Sample)[1] == N) 
    {   Missing = NA
        Zs = SNP_Pattern[,4:dim(SNP_Pattern)[2]]
        for (i in 1:N)
        { index = as.numeric(Zs[i,])
          index = index[!is.na(index)]
          Mat_Zs[i,index-1] = 1
        }
    }  else
      
    {
      Missing = IDs[which(!(IDs %in% Sample$Value ) == T)]
      Zs = SNP_Pattern[-which(SNP_Pattern$V1 %in% Missing),4:dim(SNP_Pattern)[2]]
      for (i in 1:dim(Sample)[1])
      { index = as.numeric(Zs[i,])
        index = index[!is.na(index)]
        Mat_Zs[i,index-1] = 1
      }                          
    }
    
    Sample_num_allele = numeric(dim(Sample)[1])
    if (all(is.na(Missing)) == FALSE)
    {SNP_Pattern_del_missing = SNP_Pattern[-which(SNP_Pattern$V1 %in% Missing),]
     Sample_num_allele = sapply(1:dim(Sample)[1], function(x) 
     {Sample[[as.character(SNP_Pattern_del_missing$V3)[x]]][x]})
    } else
      Sample_num_allele = sapply(1:dim(Sample)[1], function(x) 
      {Sample[[as.character(SNP_Pattern$V3)[x]]][x]})
    
    
    
    
    list(Sample = Sample, Sample_num_reads = Sample_num_reads,
         Sample_num_allele = Sample_num_allele, Mat_Zs = Mat_Zs,
         Missing = Missing   
    )
    
  }



Single_est = 
  function (Counts_Sample, SNP_Pattern, N = 102057,filename =NA)
  {
    Counts_sample_dealed = Preprocessing(Counts_Sample, SNP_Pattern, N)
    single_row = which(rowSums(Counts_sample_dealed$Mat_Zs) == 1)
    single_all_info = Counts_sample_dealed$Mat_Zs[single_row,]
    single_num_read = as.numeric(as.character(Counts_sample_dealed$Sample_num_reads[single_row] ))
    single_num_allele = as.numeric(Counts_sample_dealed$Sample_num_allele[single_row] )
    single_vari_name = sapply(1: nrow(single_all_info), function(x){which(single_all_info[x,]==1)})
    pi_est = single_num_allele/single_num_read
    single_est = data.frame(vari_name = single_vari_name, p_est = pi_est)
    single_names = single_vari_name[!duplicated(single_vari_name)]
    single_names = single_names[order(single_names)]
    
    para_est_mean = sapply(single_names,function(x)
    {mean(single_est[single_est$vari_name == x,]$p_est)})
    para_est_sd = sapply(single_names,function(x)
    {sd(single_est[single_est$vari_name == x,]$p_est)})
    
    
    para_details = lapply(single_names,function(x)
    {Counts_sample_dealed$Sample[as.numeric(row.names(single_est[single_est$vari_name == x,])),]})
    
    para_means = lapply(single_names,function(x)
    {single_est[as.numeric(row.names(single_est[single_est$vari_name == x,])),]})
    
    if (is.na(filename) == 0)
    { filename = as.character(paste(filename,".pdf",sep=""))
      pdf(filename)
      plot(single_est$vari_name,single_est$p_est,main = filename,
           xlab = "beta_i",ylab="beta_i_estimate")
      dev.off()
    }

    results=list(para_means = para_means ,
                 para_details = para_details,
                 para_est_mean = data.frame(vari_name = single_names,para_est_mean = para_est_mean) ,
                 para_est_sd =  data.frame(vari_name = single_names,para_est_sd = para_est_sd) 
                 
    )
  }






Whole_est = function(results,data)
{
  require(quadprog)
  #results = Mix1_results
  #data = Mix1_dealed
  beta = rep(-1,dim(data$Mat_Zs)[2])
  for (i in 1:length(results$para_means))
  { 
    if (max(results$para_means[[i]]$p_est)<0.05)
    {beta[unique(results$para_means[[i]]$vari_name)]=0}
    
  }
  
  #plot(results$para_means[[1]]$vari_name,
  #     results$para_means[[1]]$p_est,xlim=c(1,32),ylim=c(0,1))

  
  Mat_Zs = as.data.frame(data$Mat_Zs)
  Z_new0 = Mat_Zs[,which(beta==-1),drop = F]
  Z_new = Z_new0[which(rowSums(Z_new0) ==1),,drop=F]
  
  index = as.numeric(row.names(Z_new))
  n = as.numeric(as.character((data$Sample_num_reads)))
  x = as.numeric(data$Sample_num_allele)
  pi = x/n
  pi = pi[index]
  
  all.names = names(Z_new)
  v.names = sapply(1:dim(Z_new)[1], function(i) all.names[which(Z_new[i,] ==1)])
  aa = strsplit(v.names,split="V")
  v.names = sapply(1:length(aa),function(i) aa[[i]][2]) 
  
  est_means = data.frame(var.name = v.names, pi.est = pi)
  
  
  para_means = lapply(unique(est_means$var.name),function(x)
    est_means[which(est_means$var.name == x),])
  
  for (i in 1:length(para_means))
  { 
    if (max(para_means[[i]]$pi.est)<0.05)
    { index0 = as.numeric(as.character(unique(para_means[[i]]$var.name)))
      beta[index0]=0}
    
  }
  
  
  X = Mat_Zs[,which(beta==-1),drop=F]
  index11 = !rowSums(X)==0
  index22 = !colSums(X)==0
  ixNA = which(colSums(X)==0)
  if (length(ixNA) > 0)
  { ixNANA=sapply(strsplit(names(ixNA),"V"),function(x){x[2]})
    beta[as.numeric(ixNANA)]=NA }
  X = X[index11,index22,drop=F]
  pi = x/n
  pi = pi[index11]
  aa = lm.restricted(pi,as.matrix(X))
  
  beta[which(beta== -1)] =aa
  
  beta
}



lm.restricted=function(y,X) {
  ## X must have full column rank. deal with duplicated columns in X                                                                  
  ix=!duplicated(t(X))
  X2=X[,ix,drop=F]
  
  
  n=nrow(X2); p=ncol(X2)
  Dmat=t(X2) %*% X2
  dvec=t(X2) %*% y
  bvec=rep(0, p)
  Amat=diag(p)
  
  
  a=solve.QP(Dmat, dvec, Amat, bvec, meq=0)
  res=rep(NA, ncol(X))
  res[ix]= a[[1]]
  res
}

Uniquecounts = function(Mat_Zs,Sample_num_allele)
{
  ix1 = which(rowSums(Mat_Zs)  ==1 )
  ix2 = which(Sample_num_allele>0)
  ix = intersect(ix1,ix2)
  submat = Mat_Zs[ix,]
  submat1 = Mat_Zs[ix1,]
  
  
  data.frame(counts = colSums(submat), percent = paste(round(colSums(submat)/colSums(submat1)*100,3),sep="\t"))
}
