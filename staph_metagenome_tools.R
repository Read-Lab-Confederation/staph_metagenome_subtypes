#functions for the metagenome analysis
bintr <- function(mat,cutoff){
  mat[which(mat > cutoff)] <- 1
  mat[which(!(mat > cutoff))] <- 0
  return(mat)
}

calc_FTS <- function(pop,mini) {  
  #pop and mini are the outputs of the calc_hits funciton - 2 element vectors
  fishmat <- matrix(c(mini[1],mini[2],pop[1]-mini[1],pop[2]-mini[2]),ncol=2,nrow=2) 
  print(fishmat)
  print(c("Sum fishmat: ",sum(fishmat)))
  return(fisher.test(fishmat, alternative = 'g'))
}

calc_hits <-function(nameset,mat) {
  #takes a binary matrix, returns number of hits in lower triangle of matrix.  Can subset the larger matrix by a smaller set
  minimat <- select(mat,one_of(nameset))[rownames(mat) %in% nameset,] 
  minimat.size <- (length(minimat)^2 - length(minimat) )/2
  minimat.hits <-sum(minimat)/2
  minimat.nohits <- minimat.size - minimat.hits
  return(c(minimat.hits, minimat.nohits))
}

#this is for hamming distances - returns total hamming distances
calc_hits_slice <-function(nameset,mat) {
  #this is for hamming distances - returns total hamming distances
  minimat <- slice(mat,nameset)[,nameset] 
  minimat.hits <-sum(minimat)/2
  return(minimat.hits)
}

create_cooccur_mat <- function(mat){
  #method below from Tyler Rinker: Stack overflow
  library(reshape2)
  dat2 <- melt(mat)
  w <- dcast(dat2, V2~V1)
  x <- as.matrix(w[,-1])
  x[is.na(x)] <- 0
  x <- apply(x, 2,  function(x) as.numeric(x > 0))
  v <- x %*% t(x) 
  diag(v) <- 0                                      #repalce diagonal
  dimnames(v) <- list(w[, 1], w[,1])  
  u <- data.frame((v > 0)*1)  #convert to binary
  return(u)
}

genotypes_plot <- function(mat,tit) {
  top_genos <- c("CC_30","CC_8","CC_45","CC_398","CC_133","CC_59","CC_15","CC_97","CC_5","CC_9", "CC_22", "CC_239", "CC_1", "CC_121")
  cS <- colSums(mat)
  cStop <- cS[top_genos]
  cStop[["others"]] <- sum(cS) - sum(cStop)  
  barplot(cStop, main = tit, las = 3, cex.names = 0.8, col = (brewer.pal(12,name = 'Set3')))
}

all_genotypes_plot <- function(mat,tit) {
  cS <- colSums(mat)
  barplot(cS, main = tit, las = 3, cex.names = 0.8, col = "gray")
}


run_bs_subj_adonis <-function(df,bs_vec,subj_vec) {
  #runs a eset of adonis comparisons
  #df is a dataframe of distances
  #bs_vec and sub_vec are vectors of factors
  library(e1071)
  library(vegan)
  # look at relationships by body site
  body_site_adonis <- adonis(df ~ bs_vec)
  print(body_site_adonis)
  #test for differences in dsipersion between groups
  body_site_betadisper <- betadisper(dist(df),bs_vec)
  print(body_site_betadisper)
  p=permutest(body_site_betadisper)
  print(p$tab)
  # now by subject_id
  subject_adonis <- adonis(df ~ subj_vec)
  print(subject_adonis)
  subject_betadisper <- betadisper(dist(df),subj_vec)
  print(subject_betadisper)
  p=permutest(subject_betadisper)
  print(p$tab)
  # now interaction of body_site and subject
  body_site_subj_adonis <- adonis(df ~ bs_vec + subj_vec)
  print(body_site_subj_adonis)
}

make_subtype_matrix <-function(df) {
  library(dplyr)
  mat <- select(df,matches("CC")) %>% as.matrix
  assert_that(dim(mat)[2] == 33)
  return(mat)
}

plot_coverages <-function(combined.df,titl) {
  check_staph_df(combined.df)
  par(mar = c(12,4,4,2),cex = 0.8)
  with(combined.df,boxplot(Staph_cov ~ Body.site, border = "red", outline = FALSE, boxlty = 0, whisklty = 0, staplelty = 0, ylim = c(0,10.5), las = 2, main= titl, ylab = "S. aureus coverage"))
  with(combined.df,stripchart(Staph_cov ~ Body.site, add = TRUE,vertical = TRUE, pch= 16, method = "stack" ))
}

plot_adjusted_coverages <- function(combined.df,titl){
  check_staph_df(combined.df)
  stcols <- grep("CC|MLST",colnames(combined.df))
  adjusted_cov.df <- combined.df[,stcols]*combined.df$Staph_cov
  adjusted_cov.df <- adjusted_cov.df[,order(names(adjusted_cov.df))]
  #plot 1
  boxplot(adjusted_cov.df,border = "red", outline = FALSE, boxlty = 0, whisklty = 0, staplelty = 0, ylim = c(0,10), las = 2, main= titl,ylab = "Adjusted S. aureus coverage")
  stripchart(adjusted_cov.df, add = TRUE, vertical = TRUE, pch= 16)
}

plot_mecA <- function(combined.df,titl) {
  check_staph_df(combined.df)
  with(combined.df, plot(Staph_cov,mecA_cov, col = Body.site, pch = 16, main = titl, xlim = c(0,10)))
}

plot_diversity_vers_cov <- function (combined.df,titl){
  library(vegan)
  check_staph_df(combined.df)
  stcols <- grep("CC|MLST",colnames(combined.df))
  Shannon <- diversity(combined.df[,stcols])
  lm_Shannon <- with(combined.df,lm(Staph_cov ~ Shannon))
  print(summary(lm_Shannon))
  with(combined.df, plot(Staph_cov,Shannon, pch = 16, main = titl, xlim = c(0,10)))
  abline(lm_Shannon, col = "red")
}
# plot_CC_cov <- function(CC) {
#  BROKEn   
# assert_that(not_empty(combined))
#   CCeval <- eval(parse(text=CC))
#   t1 <- which(select(combined,CCeval) > 0.9) #need to do this to evaluated argument (http://stackoverflow.com/questions/11880906/pass-subset-argument-through-a-function-to-subset)
#   with(combined,plot(Staph_cov[t1]*CCeval[t1], mecA_cov[t1]), xlim = c(0,12), ylim = c(0,12))
# }
check_staph_df <- function(df){
  library(assertthat)
  assert_that(length(grep("Body.site",colnames(df))) == 1)
  assert_that(length(grep("Subject.Id",colnames(df))) == 1)
  assert_that(length(grep("CC",colnames(df))) == 33)
}

subject_perm <- function(df,multiSubjects,hamming_mat){
  library(gdata)
  check_staph_df(df)
  #mult
  sub1.hits = 0
  sub1.cells = 0  # total number of cells queried
  for (i in multiSubjects$Subject.Id) {
    sub1 <- which(df$Subject.Id == as.character(i))    
    sub1.hits <- calc_hits_slice(sub1,dat6) + sub1.hits
    sub1.cells <- sub1.cells + ((length(sub1)^2 - length(sub1))/2)
  }
  dat7 <- lowerTriangle(hamming_mat)
  sub1.perms <- replicate(10000, sum(dat7[sample.int(length(dat7),sub1.cells)]),simplify = 'array')
  tle <- paste("10,000 samples of",sub1.cells,"random pairs")
  hist(sub1.perms, xlab = "Hamming scores", main = tle, xlim = c(1400,1700))
  cat("Score for intraperson hits = ",sub1.hits, "\n")
  print("Quartlies for random distribution")
  print(quantile(sub1.perms))
  cat("Empirical p value = ",ecdf(sub1.perms)(sub1.hits))
}

by_factor_perm <- function(bs,df,hamming_mat){
  check_staph_df(df)
  for (i in bs) {
      bss_rows <- which(df$Body.site == i)
      if (length(bss_rows) > 3) {
      perms <- replicate(10000, calc_hits_slice(sample(x = 1:nrow(hamming_mat), size = length(bss_rows)),dat6),simplify = 'array')
      print(i)
      print(c("Number of samples ",length(bss_rows)))
      print("Distribution of random hits")
      print(quantile(perms))
      num_hits <- calc_hits_slice(bss_rows,hamming_mat)
      print(num_hits)
      cat("Empirical p value ")
      print(ecdf(perms)(num_hits))
      cat("\n")
    }
    else (cat("Zero samples in ", i))
  }
}

intra_body_FTS <- function(body1,body2,df,multiSubjects,u) {
  library(dplyr)
  check_staph_df(df)
  temp.an <- filter(df,Body.site == body1) %>% filter(Subject.Id %in% multiSubjects$Subject.Id) %>% subset(select = c(Subject.Id,Body.site,Sample.Id))
  temp.td <- filter(df,Body.site == body2) %>% filter(Subject.Id %in% multiSubjects$Subject.Id) %>% subset(select = c(Subject.Id,Body.site,Sample.Id))
  temp_an_td_join <- inner_join(temp.td,temp.an, by = "Subject.Id")
  
  #sum of all the interactions between the body sites
  if(nrow(temp_an_td_join) > 0) {
    sum_an_td_hits <- sum(u[unique(temp_an_td_join[,3]), unique(temp_an_td_join[,5])])
    an_td <- length(unique(temp_an_td_join[,3])) * length(unique(temp_an_td_join[,5]))
    an_td_no_hits <- an_td - sum_an_td_hits
    intra_sample <- length(temp_an_td_join$Sample.Id.x)
    intra_sample.hits <- sum(sapply(1:intra_sample, function(x) u[temp_an_td_join$Sample.Id.x[x],temp_an_td_join$Sample.Id.y[x]]))
    intra_sample.nohits <- intra_sample -intra_sample.hits
    print(c(body1,body2))
    fishmat <- matrix(c(intra_sample.hits,intra_sample.nohits,sum_an_td_hits-intra_sample.hits,an_td_no_hits -intra_sample.nohits),ncol = 2,nrow = 2)
    print(fishmat)
    print(c("Sum fishmat: ",sum(fishmat)))
    if(sum(fishmat) > 24)  {
      return(fisher.test(fishmat, alternative = 'g'))
    }
    else {cat("Too low overlap between sites\n\n")}
  }
  else {
    cat("No shared rare SNPs between ",body1,body2)
    cat("\n")
    }
}

merge_CCs <- function(in_data,CC) {
  new_col <- select(in_data,matches(CC)) %>% rowSums()
  in_data <- select(in_data,-(matches(CC)))
  in_data <- mutate(in_data,temp = new_col)
  names(in_data)[names(in_data) == 'temp'] <- CC
  return(in_data)
}

plot_CC_types <-function(CC, CCcol = "red", mat, SRA_file, map11, map10, plotdir, cutoff = 0.2) {
  library(RgoogleMaps)
  crows <- which(SRA_file[[CC]] > cutoff)
  #cruns <-select(mat, Run) %>% slice(crows) 
  if (length(crows) > 0) {
    CC_df <- slice(SRA_file, crows) %>% select(Latitude,Logitude)
    plotname <- paste(plotdir,CC,"_z11.png",sep= "", collapse = "")
    png(plotname,width=640, height =640, res = 75)
    PlotOnStaticMap(map11, lat = as.numeric(CC_df$Latitude) , lon = as.numeric(CC_df$Logitude), cex=1.5,pch=20, col = CCcol)
    dev.off()
    
    plotname <- paste(plotdir,CC,"_z10.png",sep= "", collapse = "")
    png(plotname,width=640, height =640, res = 75)
    PlotOnStaticMap(map10, lat = as.numeric(CC_df$Latitude) , lon = as.numeric(CC_df$Logitude), cex=1.5,pch=20, col = CCcol)
    dev.off() 
  }
  
}

avg_geog_dist <- function(p){
  mat <- as.data.frame(combinations(nrow(p),2))
  dist_vec <- sapply(1:nrow(mat), function(x) distance.chord(as.numeric(p[mat[x,1],]),as.numeric(p[mat[x,2],])))
  return(mean(dist_vec))
}

rand_distances <- function(n,latlon,perms = 1000) {
  res_vec <- replicate(perms,avg_geog_dist(sample_n(latlon,n)),simplify = "vector")
  return(res_vec)
  
}

#below from http://www.biostat.umn.edu/~sudiptob/Software/distonearth.R
distance.chord <- function(point1,point2){
  
  R <- 6371
  
  p1rad <- point1 * pi/180
  p2rad <- point2 * pi/180
  
  lat <- p1rad[2]
  lon <- p1rad[1]
  
  u1 <- c(cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat))
  
  lat <- p2rad[2]
  lon <- p2rad[1]
  
  u2 <- c(cos(lat)*cos(lon),cos(lat)*sin(lon),sin(lat))
  
  R*sqrt(sum((u1-u2)^2))
  
}

CC_geog_perm_test <- function(SRA_file, CC, cutoff, s = 234523, reps= 1000) {
  crows <- which(SRA_file[[CC]] > cutoff)
  CC_df <- slice(SRA_file, crows) %>% select(Latitude,Logitude)
  av_geog <- avg_geog_dist(CC_df)
  cat("Av geog. distance of ",CC," = ",av_geog,"\n")
  set.seed(s)
  rand_dists <- rand_distances(length(crows),select(SRA_file,Latitude,Logitude),reps)
  cat("Quartile of random distribution of dists ",ecdf(rand_dists)(av_geog),"\n")
}

decorate_staph_tree <- function(CC, tree, strains, cutoff = 0.65, deco = "red"){
  tag_list <- filter(strains,grepl(CC,Reference.CC)) %>% filter(Beta > cutoff) %>% select(Sample.Id.of.0.75X)
  tps <- which(tree$tip.label %in% tag_list$Sample.Id.of.0.75X)
  plot(tree, "unrooted", show.tip.label = FALSE, main = CC)
  tiplabels(tip = tps, pch= 20, col = deco)
}

dist_between_stations <- function(pairs,geog.mat){
  #pairs are piars of points
  #geor mat contains cols with point,longditude, latitue
  p1 <- filter(geog.mat,Run == pairs[1]) %>% select(Logitude,Latitude) %>% t() %>% as.vector() %>% as.numeric()
  p2 <- filter(geog.mat,Run == pairs[2]) %>% select(Logitude,Latitude) %>% t() %>% as.vector() %>% as.numeric()
  return(distance.chord(p1,p2))
}
