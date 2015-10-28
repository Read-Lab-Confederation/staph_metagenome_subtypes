---
title: "Distribution analysis of 0.5X filtered data"
author: "Tim Read"
date: "August 13, 2015"
output: html_document
---
Starting form the filtered table from 'HMP_coverage.Rmd'.  Run a series of analysis to look at relationships between body site and subjects.



```r
print(date())
```

```
## [1] "Wed Oct 28 16:17:35 2015"
```

```r
library(reshape2)
#library(igraph)
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
## 
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## 
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
#library(biomod2)
library(e1071)
library(RColorBrewer)
library(gdata)
```

```
## gdata: read.xls support for 'XLS' (Excel 97-2004) files ENABLED.
## 
## gdata: read.xls support for 'XLSX' (Excel 2007+) files ENABLED.
## 
## Attaching package: 'gdata'
## 
## The following objects are masked from 'package:dplyr':
## 
##     combine, first, last
## 
## The following object is masked from 'package:stats':
## 
##     nobs
## 
## The following object is masked from 'package:utils':
## 
##     object.size
```

```r
library(vegan)
```

```
## Loading required package: permute
## Loading required package: lattice
## This is vegan 2.3-0
```

```r
library(assertthat)
source('./staph_metagenome_tools.R', echo=TRUE)
```

```
## 
## > bintr <- function(mat, cutoff) {
## +     mat[which(mat > cutoff)] <- 1
## +     mat[which(!(mat > cutoff))] <- 0
## +     return(mat)
## + }
## 
## > calc_FTS <- function(pop, mini) {
## +     fishmat <- matrix(c(mini[1], mini[2], pop[1] - mini[1], pop[2] - 
## +         mini[2]), ncol = 2, nrow = 2)
## +  .... [TRUNCATED] 
## 
## > calc_hits <- function(nameset, mat) {
## +     minimat <- select(mat, one_of(nameset))[rownames(mat) %in% 
## +         nameset, ]
## +     minimat.size <- ( .... [TRUNCATED] 
## 
## > calc_hits_slice <- function(nameset, mat) {
## +     minimat <- slice(mat, nameset)[, nameset]
## +     minimat.hits <- sum(minimat)/2
## +     return(minima .... [TRUNCATED] 
## 
## > create_cooccur_mat <- function(mat) {
## +     library(reshape2)
## +     dat2 <- melt(mat)
## +     w <- dcast(dat2, V2 ~ V1)
## +     x <- as.matrix(w[, -1])
##  .... [TRUNCATED] 
## 
## > genotypes_plot <- function(mat, tit) {
## +     top_genos <- c("CC_30", "CC_8", "CC_45", "CC_398", "CC_133", 
## +         "CC_59", "CC_15", "CC_97", "CC_ ..." ... [TRUNCATED] 
## 
## > all_genotypes_plot <- function(mat, tit) {
## +     cS <- colSums(mat)
## +     barplot(cS, main = tit, las = 3, cex.names = 0.8, col = "gray")
## + }
## 
## > run_bs_subj_adonis <- function(df, bs_vec, subj_vec) {
## +     library(e1071)
## +     library(vegan)
## +     body_site_adonis <- adonis(df ~ bs_vec)
## +     .... [TRUNCATED] 
## 
## > make_subtype_matrix <- function(df) {
## +     library(dplyr)
## +     mat <- select(df, matches("CC")) %>% as.matrix
## +     assert_that(dim(mat)[2] == 33) .... [TRUNCATED] 
## 
## > plot_coverages <- function(combined.df, titl) {
## +     check_staph_df(combined.df)
## +     par(mar = c(12, 4, 4, 2), cex = 0.8)
## +     with(combined.df, .... [TRUNCATED] 
## 
## > plot_adjusted_coverages <- function(combined.df, titl) {
## +     check_staph_df(combined.df)
## +     stcols <- grep("CC|MLST", colnames(combined.df))
## +  .... [TRUNCATED] 
## 
## > plot_mecA <- function(combined.df, titl) {
## +     check_staph_df(combined.df)
## +     with(combined.df, plot(Staph_cov, mecA_cov, col = Body.site, 
## +   .... [TRUNCATED] 
## 
## > plot_diversity_vers_cov <- function(combined.df, titl) {
## +     library(vegan)
## +     check_staph_df(combined.df)
## +     stcols <- grep("CC|MLST", coln .... [TRUNCATED] 
## 
## > check_staph_df <- function(df) {
## +     library(assertthat)
## +     assert_that(length(grep("Body.site", colnames(df))) == 1)
## +     assert_that(length( .... [TRUNCATED] 
## 
## > subject_perm <- function(df, multiSubjects, hamming_mat) {
## +     library(gdata)
## +     check_staph_df(df)
## +     sub1.hits = 0
## +     sub1.cells = 0
## +  .... [TRUNCATED] 
## 
## > by_factor_perm <- function(bs, df, hamming_mat) {
## +     check_staph_df(df)
## +     for (i in bs) {
## +         bss_rows <- which(df$Body.site == i)
## +    .... [TRUNCATED] 
## 
## > intra_body_FTS <- function(body1, body2, df, multiSubjects, 
## +     u) {
## +     library(dplyr)
## +     check_staph_df(df)
## +     temp.an <- filter(df, Bo .... [TRUNCATED] 
## 
## > merge_CCs <- function(in_data, CC) {
## +     new_col <- select(in_data, matches(CC)) %>% rowSums()
## +     in_data <- select(in_data, -(matches(CC)))
## +  .... [TRUNCATED] 
## 
## > plot_CC_types <- function(CC, CCcol = "red", mat, 
## +     SRA_file, map11, map10, plotdir, cutoff = 0.2) {
## +     library(RgoogleMaps)
## +     crows <-  .... [TRUNCATED] 
## 
## > avg_geog_dist <- function(p) {
## +     mat <- as.data.frame(combinations(nrow(p), 2))
## +     dist_vec <- sapply(1:nrow(mat), function(x) distance.chord .... [TRUNCATED] 
## 
## > rand_distances <- function(n, latlon, perms = 1000) {
## +     res_vec <- replicate(perms, avg_geog_dist(sample_n(latlon, 
## +         n)), simplify = "v ..." ... [TRUNCATED] 
## 
## > distance.chord <- function(point1, point2) {
## +     R <- 6371
## +     p1rad <- point1 * pi/180
## +     p2rad <- point2 * pi/180
## +     lat <- p1rad[2]
## +   .... [TRUNCATED] 
## 
## > CC_geog_perm_test <- function(SRA_file, CC, cutoff, 
## +     s = 234523, reps = 1000) {
## +     crows <- which(SRA_file[[CC]] > cutoff)
## +     CC_df <- s .... [TRUNCATED] 
## 
## > decorate_staph_tree <- function(CC, tree, strains, 
## +     cutoff = 0.65, deco = "red") {
## +     tag_list <- filter(strains, grepl(CC, Reference.CC))  .... [TRUNCATED] 
## 
## > dist_between_stations <- function(pairs, geog.mat) {
## +     p1 <- filter(geog.mat, Run == pairs[1]) %>% select(Logitude, 
## +         Latitude) %>% t() .... [TRUNCATED] 
## 
## > H_distance_between_stations <- function(pairs, mat) {
## +     h1 <- filter(mat, Run == pairs[1])[, 2:ncol(mat)] %>% t() %>% 
## +         as.vector()
## +   .... [TRUNCATED]
```

###Read in data file created in earlier pipeline


```r
dat4 <- read.table("./Data/cov0.5")
```

###Create data files

```r
#list of all subjects with more than one sample
multiSubjects <- count(dat4,Subject.Id) %>% filter(n > 1) %>% select(Subject.Id ) 
dat5 <- make_subtype_matrix(dat4)
#create Hamming dist matrices with and without cutof  min value of 0.2
dat4$Subject.Id <- as.factor(dat4$Subject.Id)
dat6 <- make_subtype_matrix(dat4) %>% bintr(0.2) %>% hamming.distance %>% data.frame 
dat8 <- make_subtype_matrix(dat4) %>% hamming.distance %>% data.frame 
```
### PERMANOVA

test for significant associations of subtype with with bodysite and subject.  us e Hamming dist. matrix. Two levels, one with a beta cutoff for all samples > 0.2 and one without

```r
set.seed(344098)
run_bs_subj_adonis(dat6,dat4$Body.site,dat4$Subject.Id)
```

```
## 
## Call:
## adonis(formula = df ~ bs_vec) 
## 
## Permutation: free
## Number of permutations: 999
## 
## Terms added sequentially (first to last)
## 
##            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
## bs_vec      8    0.5059 0.063235  2.3899 0.12016  0.002 **
## Residuals 140    3.7043 0.026459         0.87984          
## Total     148    4.2102                  1.00000          
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## 	Homogeneity of multivariate dispersions
## 
## Call: betadisper(d = dist(df), group = bs_vec)
## 
## No. of Positive Eigenvalues: 26
## No. of Negative Eigenvalues: 0
## 
## Average distance to median:
##               anterior nares attached keratinized gingiva 
##                       11.043                        5.292 
##                buccal mucosa                  hard palate 
##                        8.523                        0.000 
##   left retroauricular crease  right retroauricular crease 
##                       10.414                       11.862 
##                        stool         supragingival plaque 
##                        0.000                        8.000 
##                tongue dorsum 
##                       11.617 
## 
## Eigenvalues for PCoA axes:
##      PCoA1      PCoA2      PCoA3      PCoA4      PCoA5      PCoA6 
## 11128.1137  5169.6086  1500.9591  1299.3574   859.2078   652.1567 
##      PCoA7      PCoA8 
##   402.3531   316.9454 
##            Df    Sum Sq  Mean Sq        F N.Perm Pr(>F)
## Groups      8  462.3259 57.79074 4.212607    999  0.014
## Residuals 140 1920.5931 13.71852       NA     NA     NA
## 
## Call:
## adonis(formula = df ~ subj_vec) 
## 
## Permutation: free
## Number of permutations: 999
## 
## Terms added sequentially (first to last)
## 
##            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
## subj_vec   81    2.5697 0.031725  1.2957 0.61036  0.025 *
## Residuals  67    1.6405 0.024485         0.38964         
## Total     148    4.2102                  1.00000         
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## 	Homogeneity of multivariate dispersions
## 
## Call: betadisper(d = dist(df), group = subj_vec)
## 
## No. of Positive Eigenvalues: 26
## No. of Negative Eigenvalues: 0
## 
## Average distance to median:
##         1         2         5         6         8         9        10 
## 0.000e+00 0.000e+00 0.000e+00 0.000e+00 6.055e-14 0.000e+00 0.000e+00 
##        11        12        13        15        18        19        20 
## 0.000e+00 8.413e+00 7.722e+00 1.230e+01 8.537e+00 0.000e+00 0.000e+00 
##        21        22        24        25        26        27        28 
## 0.000e+00 0.000e+00 0.000e+00 0.000e+00 6.110e+00 0.000e+00 8.322e+00 
##        29        30        31        32        33        34        35 
## 9.657e+00 6.874e+00 0.000e+00 1.192e+01 1.092e+01 8.944e+00 7.410e+00 
##        37        38        39        40        41        42        43 
## 9.631e+00 3.873e+00 1.034e+01 1.034e+01 7.009e+00 0.000e+00 6.103e+00 
##        44        45        46        47        49        51        52 
## 5.292e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 
##        54        58        60        62        63        64        65 
## 0.000e+00 0.000e+00 5.935e+00 0.000e+00 1.285e+01 0.000e+00 0.000e+00 
##        66        67        68        69        70        72        73 
## 0.000e+00 0.000e+00 0.000e+00 1.257e+01 0.000e+00 1.027e+01 1.319e+01 
##        74        75        76        79        80        81        82 
## 0.000e+00 0.000e+00 6.095e+00 3.389e-14 0.000e+00 3.873e+00 0.000e+00 
##        83        84        85        86        88        89        90 
## 3.255e-14 6.252e-14 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 
##        93        94        95        96        97        98       100 
## 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 
##       102       103       104       107       110 
## 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 
## 
## Eigenvalues for PCoA axes:
##      PCoA1      PCoA2      PCoA3      PCoA4      PCoA5      PCoA6 
## 11128.1137  5169.6086  1500.9591  1299.3574   859.2078   652.1567 
##      PCoA7      PCoA8 
##   402.3531   316.9454 
##           Df   Sum Sq  Mean Sq        F N.Perm Pr(>F)
## Groups    81 3294.148 40.66850 1.113509    999  0.398
## Residuals 67 2447.029 36.52283       NA     NA     NA
## 
## Call:
## adonis(formula = df ~ bs_vec + subj_vec) 
## 
## Permutation: free
## Number of permutations: 999
## 
## Terms added sequentially (first to last)
## 
##            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
## bs_vec      8    0.5059 0.063235  2.6677 0.12016  0.002 **
## subj_vec   80    2.2821 0.028526  1.2034 0.54203  0.075 . 
## Residuals  60    1.4222 0.023704         0.33781          
## Total     148    4.2102                  1.00000          
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
run_bs_subj_adonis(dat8,dat4$Body.site,dat4$Subject.Id)
```

```
## 
## Call:
## adonis(formula = df ~ bs_vec) 
## 
## Permutation: free
## Number of permutations: 999
## 
## Terms added sequentially (first to last)
## 
##            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
## bs_vec      8    1.2090 0.151120  12.344 0.41363  0.001 ***
## Residuals 140    1.7139 0.012242         0.58637           
## Total     148    2.9228                  1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## 	Homogeneity of multivariate dispersions
## 
## Call: betadisper(d = dist(df), group = bs_vec)
## 
## No. of Positive Eigenvalues: 145
## No. of Negative Eigenvalues: 0
## 
## Average distance to median:
##               anterior nares attached keratinized gingiva 
##                       57.438                        6.218 
##                buccal mucosa                  hard palate 
##                       10.928                        0.000 
##   left retroauricular crease  right retroauricular crease 
##                       25.589                       20.978 
##                        stool         supragingival plaque 
##                        0.000                        5.099 
##                tongue dorsum 
##                       17.546 
## 
## Eigenvalues for PCoA axes:
##       PCoA1       PCoA2       PCoA3       PCoA4       PCoA5       PCoA6 
## 250759.9737  11740.0204   3185.8316   2129.2399   1629.7697   1494.0862 
##       PCoA7       PCoA8 
##   1192.9320    938.1961 
##            Df   Sum Sq   Mean Sq        F N.Perm Pr(>F)
## Groups      8 32230.11 4028.7641 8.703185    999  0.002
## Residuals 140 64806.96  462.9069       NA     NA     NA
## 
## Call:
## adonis(formula = df ~ subj_vec) 
## 
## Permutation: free
## Number of permutations: 999
## 
## Terms added sequentially (first to last)
## 
##            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
## subj_vec   81    1.6469 0.020332  1.0676 0.56345  0.398
## Residuals  67    1.2760 0.019044         0.43655       
## Total     148    2.9228                  1.00000       
## 
## 	Homogeneity of multivariate dispersions
## 
## Call: betadisper(d = dist(df), group = subj_vec)
## 
## No. of Positive Eigenvalues: 145
## No. of Negative Eigenvalues: 0
## 
## Average distance to median:
##         1         2         5         6         8         9        10 
## 0.000e+00 0.000e+00 0.000e+00 0.000e+00 4.187e+01 0.000e+00 0.000e+00 
##        11        12        13        15        18        19        20 
## 0.000e+00 2.825e+01 5.394e+01 2.108e+01 2.994e+01 0.000e+00 0.000e+00 
##        21        22        24        25        26        27        28 
## 0.000e+00 0.000e+00 0.000e+00 0.000e+00 9.317e+00 0.000e+00 2.247e+01 
##        29        30        31        32        33        34        35 
## 7.616e+00 1.704e+01 0.000e+00 2.051e+01 9.811e+00 6.312e+01 2.465e+01 
##        37        38        39        40        41        42        43 
## 4.441e+01 1.859e+01 2.010e+01 1.920e+01 1.413e+01 0.000e+00 1.230e+02 
##        44        45        46        47        49        51        52 
## 5.305e+01 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 
##        54        58        60        62        63        64        65 
## 0.000e+00 0.000e+00 9.454e+00 0.000e+00 2.662e+01 0.000e+00 0.000e+00 
##        66        67        68        69        70        72        73 
## 0.000e+00 0.000e+00 0.000e+00 8.201e+00 0.000e+00 1.496e+01 6.964e+00 
##        74        75        76        79        80        81        82 
## 0.000e+00 0.000e+00 1.132e+01 7.036e+00 0.000e+00 7.018e+00 0.000e+00 
##        83        84        85        86        88        89        90 
## 2.888e-13 1.819e+01 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 
##        93        94        95        96        97        98       100 
## 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 
##       102       103       104       107       110 
## 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 
## 
## Eigenvalues for PCoA axes:
##       PCoA1       PCoA2       PCoA3       PCoA4       PCoA5       PCoA6 
## 250759.9737  11740.0204   3185.8316   2129.2399   1629.7697   1494.0862 
##       PCoA7       PCoA8 
##   1192.9320    938.1961 
##           Df   Sum Sq  Mean Sq       F N.Perm Pr(>F)
## Groups    81 60889.57 751.7230 1.46006    999  0.251
## Residuals 67 34495.47 514.8578      NA     NA     NA
## 
## Call:
## adonis(formula = df ~ bs_vec + subj_vec) 
## 
## Permutation: free
## Number of permutations: 999
## 
## Terms added sequentially (first to last)
## 
##            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
## bs_vec      8   1.20896 0.151120 11.9318 0.41363  0.001 ***
## subj_vec   80   0.95396 0.011924  0.9415 0.32638  0.602    
## Residuals  60   0.75992 0.012665         0.25999           
## Total     148   2.92283                  1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
### Permutation tests

```r
#test ffor whether intra-subject distance greater than intersubject
subject_perm(dat4,multiSubjects,dat6)
```

![plot of chunk perm_tests](figure/perm_tests-1.png) 

```
## Score for intraperson hits =  328 
## [1] "Quartlies for random distribution"
##   0%  25%  50%  75% 100% 
##  319  365  375  384  421 
## Empirical p value =  5e-04
```

```r
#now look at the same test between body sites
bs <- levels(dat4$Body.site)
by_factor_perm(bs,dat4,dat6)
```

```
## [1] "anterior nares"
## [1] "Number of samples " "22"                
## [1] "Distribution of random hits"
##   0%  25%  50%  75% 100% 
##  357  537  574  610  767 
## [1] 447
## Empirical p value [1] 0.0116
## 
## Zero samples in  attached keratinized gingiva[1] "buccal mucosa"
## [1] "Number of samples " "9"                 
## [1] "Distribution of random hits"
##   0%  25%  50%  75% 100% 
##   16   80   90  100  142 
## [1] 66
## Empirical p value [1] 0.0635
## 
## Zero samples in  hard palate[1] "left retroauricular crease"
## [1] "Number of samples " "15"                
## [1] "Distribution of random hits"
##   0%  25%  50%  75% 100% 
##  110  242  262  282  368 
## [1] 266
## Empirical p value [1] 0.5723
## 
## [1] "right retroauricular crease"
## [1] "Number of samples " "23"                
## [1] "Distribution of random hits"
##   0%  25%  50%  75% 100% 
##  382  590  630  668  826 
## [1] 690
## Empirical p value [1] 0.8647
## 
## Zero samples in  stoolZero samples in  supragingival plaque[1] "tongue dorsum"
## [1] "Number of samples " "73"                
## [1] "Distribution of random hits"
##   0%  25%  50%  75% 100% 
## 5512 6350 6524 6696 7472 
## [1] 6482
## Empirical p value [1] 0.4379
```
###Plots of subtype distribution

```r
presence_mat <- as.data.frame(bintr(dat5,0.2))
top_score_mat <- as.data.frame(bintr(dat5,0.5))
# png("~/Dropbox/ARTICLES_BY_TDR/2015-staph-metagenome/HMP_barchart.png",width=640, height =640, res = 75)
# dev.off()
genotypes_plot(presence_mat,"All samples, subtypes present > 0.2")
```

![plot of chunk 0.5Xsubtype_plots](figure/0.5Xsubtype_plots-1.png) 

```r
genotypes_plot(top_score_mat,"All samples, subtypes present > 0.5")
```

![plot of chunk 0.5Xsubtype_plots](figure/0.5Xsubtype_plots-2.png) 

```r
all_genotypes_plot(presence_mat,"All CCs, 0.5X cutoff, subtypes present > 0.2")
```

![plot of chunk 0.5Xsubtype_plots](figure/0.5Xsubtype_plots-3.png) 

```r
all_genotypes_plot(top_score_mat,"All CCs, 0.5X cutoff, subtypes present > 0.5")
```

![plot of chunk 0.5Xsubtype_plots](figure/0.5Xsubtype_plots-4.png) 

```r
for (i in bs) {
  bss_rows <- which(dat4$Body.site == i)
  if(length(bss_rows) > 0) {
    bs_df <- slice(presence_mat,bss_rows)
    genotypes_plot(bs_df,paste(">0.2 beta: ", i))
  }
}
```

![plot of chunk 0.5Xsubtype_plots](figure/0.5Xsubtype_plots-5.png) ![plot of chunk 0.5Xsubtype_plots](figure/0.5Xsubtype_plots-6.png) ![plot of chunk 0.5Xsubtype_plots](figure/0.5Xsubtype_plots-7.png) ![plot of chunk 0.5Xsubtype_plots](figure/0.5Xsubtype_plots-8.png) ![plot of chunk 0.5Xsubtype_plots](figure/0.5Xsubtype_plots-9.png) ![plot of chunk 0.5Xsubtype_plots](figure/0.5Xsubtype_plots-10.png) ![plot of chunk 0.5Xsubtype_plots](figure/0.5Xsubtype_plots-11.png) ![plot of chunk 0.5Xsubtype_plots](figure/0.5Xsubtype_plots-12.png) ![plot of chunk 0.5Xsubtype_plots](figure/0.5Xsubtype_plots-13.png) 

```r
for (i in bs) {
  bss_rows <- which(dat4$Body.site == i)
  if(length(bss_rows) > 0) {
    bs_df <- slice(top_score_mat,bss_rows)
    genotypes_plot(bs_df,paste(">0.5 beta: ", i))
  }
}
```

![plot of chunk 0.5Xsubtype_plots](figure/0.5Xsubtype_plots-14.png) ![plot of chunk 0.5Xsubtype_plots](figure/0.5Xsubtype_plots-15.png) ![plot of chunk 0.5Xsubtype_plots](figure/0.5Xsubtype_plots-16.png) ![plot of chunk 0.5Xsubtype_plots](figure/0.5Xsubtype_plots-17.png) ![plot of chunk 0.5Xsubtype_plots](figure/0.5Xsubtype_plots-18.png) ![plot of chunk 0.5Xsubtype_plots](figure/0.5Xsubtype_plots-19.png) ![plot of chunk 0.5Xsubtype_plots](figure/0.5Xsubtype_plots-20.png) ![plot of chunk 0.5Xsubtype_plots](figure/0.5Xsubtype_plots-21.png) ![plot of chunk 0.5Xsubtype_plots](figure/0.5Xsubtype_plots-22.png) 
###PCA

```r
par(mfrow=c(2,2))
pcobj <- prcomp(dat6)
tr_gray <- rgb(0.5,.5,.5,.15)

for (i in bs) {
  prcols <- rep(tr_gray,nrow(dat6))
  prcols[which(dat4$Body.site == i)] <- "red"
  plot(pcobj$x,col = prcols, pch = 16, main = i)
}
```

![plot of chunk 0.5XPCA](figure/0.5XPCA-1.png) ![plot of chunk 0.5XPCA](figure/0.5XPCA-2.png) 

```r
for (i in multiSubjects$Subject.Id) {
  sub_rows = which(dat4$Subject.Id == as.character(i))
  if (length(sub_rows) > 3){
    prcols <- rep(tr_gray,nrow(dat6))
    prcols[sub_rows] <- "blue"
    plot(pcobj$x,col = prcols, pch = 16, main = c("Subject",i))
  }
}
```

![plot of chunk 0.5XPCA](figure/0.5XPCA-3.png) ![plot of chunk 0.5XPCA](figure/0.5XPCA-4.png) ![plot of chunk 0.5XPCA](figure/0.5XPCA-5.png) 
###Session Info

```r
sessionInfo()
```

```
## R version 3.2.1 (2015-06-18)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: OS X 10.10.5 (Yosemite)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  base     
## 
## other attached packages:
##  [1] assertthat_0.1     vegan_2.3-0        lattice_0.20-33   
##  [4] permute_0.8-4      gdata_2.17.0       RColorBrewer_1.1-2
##  [7] e1071_1.6-7        dplyr_0.4.2        reshape2_1.4.1    
## [10] knitr_1.11        
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.0     cluster_2.0.3   magrittr_1.5    MASS_7.3-44    
##  [5] R6_2.1.1        stringr_1.0.0   plyr_1.8.3      tools_3.2.1    
##  [9] parallel_3.2.1  grid_3.2.1      nlme_3.1-122    mgcv_1.8-7     
## [13] DBI_0.3.1       class_7.3-14    gtools_3.5.0    lazyeval_0.1.10
## [17] Matrix_1.2-2    formatR_1.2     evaluate_0.7.2  stringi_0.5-5  
## [21] methods_3.2.1
```

