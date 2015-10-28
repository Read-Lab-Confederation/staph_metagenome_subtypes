---
title: "NYC Data analysis"
author: "Tim Read"
date: "August 31, 2015"
output: html_document
---

Analysis of NYC subway data.  Based on a older script called 'NYsubwaySa.R'.


```r
print(date())
```

```
## [1] "Wed Oct 28 16:51:36 2015"
```

```r
library(RgoogleMaps)
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
library(vegan)
```

```
## Loading required package: permute
## Loading required package: lattice
## This is vegan 2.3-0
```

```r
library(ade4)
```

```
## 
## Attaching package: 'ade4'
## 
## The following object is masked from 'package:vegan':
## 
##     cca
```

```r
library(e1071)
library(gtools)
```

```
## 
## Attaching package: 'gtools'
## 
## The following object is masked from 'package:e1071':
## 
##     permutations
## 
## The following object is masked from 'package:permute':
## 
##     permute
```

```r
#library(biomod2)
library(RColorBrewer)
library(assertthat)
source('./staph_metagenome_tools.R')
```

### Load files


```r
#Public data from the original NYC subway publication: Afshinnekoo E, Meydan C, Chowdhury S, Jaroudi D, Boyer C, Bernstein N, Maritz JM, Reeves D, Gandara J, Chhangawala S, Ahsanuddin S, Simmons A, Nessel T, Sundaresh B, Pereira E, Jorgensen E, Kolokotronis S-O, Kirchberger N, Garcia I, Gandara D, Dhanraj S, Nawrin T, Saletore Y, Alexander N, Vijay P, Hénaff EM, Zumbo P, Walsh M, O’Mullan GD, Tighe S, Dudley JT, Dunaif A, Ennis S, O’Halloran E, Magalhaes TR, Boone B, Jones AL, Muth TR, Paolantonio KS, Alter E, Schadt EE, Garbarino J, Prill RJ, Carlton JM, Levy S, Mason CE. Geospatial Resolution of Human and Bacterial Diversity with City-Scale Metagenomics. Cell Systems [Internet]. Elsevier; 2015 Jul 29;1(1):72–87.
NYCdata <- (read.csv("./Data/DataTable5-metaphlan-metadata_v19 .csv",stringsAsFactors = FALSE, header = TRUE))[1:4]
```

```
## Warning in file(file, "rt"): cannot open file './Data/DataTable5-metaphlan-
## metadata_v19 .csv': No such file or directory
```

```
## Error in file(file, "rt"): cannot open the connection
```

```r
# we made this table from parsing the SRA
strain_SRA <- read.table("./Data/runs-to-samples.txt", header = TRUE)
colnames(strain_SRA) <- c("Run", "Sample.ID") #to make join easier
NYCdata_SRA <- left_join(NYCdata, strain_SRA, by = "Sample.ID")
```

```
## Error in left_join(NYCdata, strain_SRA, by = "Sample.ID"): object 'NYCdata' not found
```


```r
#Our table of binstrain results
# Staph_betas <- read.csv("./Data/Final_Beta_NYC_STAPH.csv",stringsAsFactors = FALSE, header = TRUE)
# #remove the final row, which is actually a bead wash control (suggests some carrover contamination)
# Staph_betas <- Staph_betas[-130,]
```

### import CC, merge CC, rename and tidy, create various forms of the data to be fed into funciotns 
(reflects the piecemeal development of the fucnitons over several months with changing specs)

(Note used some similar commands in HMP_coverage.Rmd).  Need to merge becasue we originally separated CC30, CC8 and CC5 into multiple groups but this turned out not to be specific.


```r
Staph_betas <- read.csv("./Data/FInal_NYC_Staph_MecA_Table.csv", stringsAsFactors = FALSE, header = TRUE) %>% filter(Staph_Coverage > 0.025) %>% filter(!(Sample.Id == "SRR1750088")) 
colnames(Staph_betas)[colnames(Staph_betas) == "CC_8_72"] <- "CC_72"
Staph_betas <- merge_CCs(Staph_betas,"CC_8_")
Staph_betas <- merge_CCs(Staph_betas,"CC_30_")
Staph_betas <- merge_CCs(Staph_betas,"CC_5_")

names(Staph_betas)[names(Staph_betas) == "MLST_93"] <- "CC_93"
Staph_betas <- Staph_betas[,mixedsort(colnames(Staph_betas))]
colnames(Staph_betas) <- gsub("(CC_.{0,3})_.{0,4}$","\\1",colnames(Staph_betas))
colnames(Staph_betas)[colnames(Staph_betas) == "Sample.Id"] <- "Run"
#remove SRR1748847, which has an incorrect coord
Staph_betas <- filter(Staph_betas, Run != "SRR1748847")
staph_mat <- make_subtype_matrix(Staph_betas)

staph_df <- as.data.frame(staph_mat)
staph_df <- cbind(staph_df,Staph_betas$Run)
colnames(staph_df)[colnames(staph_df) == "Staph_betas$Run"] <- "Run"
staph_df$Run <- as.factor(staph_df$Run)
staph_df_coords <- inner_join(NYCdata_SRA, staph_df, by = "Run") #one of the samples does not have a lat , lon posiiotn
```

```
## Error in inner_join(NYCdata_SRA, staph_df, by = "Run"): object 'NYCdata_SRA' not found
```

```r
coords_df_subtype_mat <- make_subtype_matrix(staph_df_coords)
```

```
## Error in select_(.data, .dots = lazyeval::lazy_dots(...)): object 'staph_df_coords' not found
```

### genotype plots


```r
presence_mat <- as.data.frame(bintr(staph_mat,0.2))
top_score_mat <- as.data.frame(bintr(staph_mat,0.5))
# png("~/Dropbox/ARTICLES_BY_TDR/2015-staph-metagenome/HMP_barchart.png",width=640, height =640, res = 75)
# dev.off()
all_genotypes_plot(presence_mat,"All NYC samples, subtypes present > 0.2")
```

![plot of chunk NYC_genotypes](figure/NYC_genotypes-1.png) 

```r
all_genotypes_plot(top_score_mat,"All NYC samples, subtypes present > 0.5")
```

![plot of chunk NYC_genotypes](figure/NYC_genotypes-2.png) 

### Color stations reporting S. aureus

```r
locations_with_runs <- filter(NYCdata_SRA, !is.na(Run) ) %>%
  select(Latitude, Logitude)
```

```
## Error in filter_(.data, .dots = lazyeval::lazy_dots(...)): object 'NYCdata_SRA' not found
```

```r
lats <- as.numeric(as.character(locations_with_runs$Latitude))
```

```
## Error in eval(expr, envir, enclos): object 'locations_with_runs' not found
```

```r
lons <- as.numeric(as.character(locations_with_runs$Logitude))
```

```
## Error in eval(expr, envir, enclos): object 'locations_with_runs' not found
```

```r
cols_staph <- rep(NULL,length(lats))
```

```
## Error in eval(expr, envir, enclos): object 'lats' not found
```

```r
#all the statiitons which had S. aureus
cols_staph[which(NYCdata_SRA$Run %in% Staph_betas$Run)] <- "red"
```

```
## Error in cols_staph[which(NYCdata_SRA$Run %in% Staph_betas$Run)] <- "red": object 'cols_staph' not found
```

```r
alllats <- as.numeric(as.character(NYCdata$Latitude)) #every non-numeric value is converted to NA
```

```
## Error in eval(expr, envir, enclos): object 'NYCdata' not found
```

```r
alllons <- as.numeric(as.character(NYCdata$Logitude))
```

```
## Error in eval(expr, envir, enclos): object 'NYCdata' not found
```
### Generate google maps

```r
gmap1 <- GetMap(center = c(lat = 40.7127, lon = -74.0059), size = c(640, 640), zoom = 11, GRAYSCALE = TRUE)
#centered on Queens
gmap2 <- GetMap(center = c(lat = 40.7500, lon = -73.8667), size = c(640, 640), zoom = 11, GRAYSCALE = TRUE)
gmap3 <- GetMap(center = c(lat = 40.7500, lon = -73.8667), size = c(640, 640), zoom = 10, GRAYSCALE = TRUE)
```
### Make plots of overall coverage
 

```r
png("./NYC_subway_plots/all_stations_z11.png",width=640, height =640, res = 75)
PlotOnStaticMap(gmap2, lat = lats , lon = lons, cex=1.5,pch=20, col = "blue")
```

```
## Error in LatLon2XY(lat, lon, zoom): object 'lats' not found
```

```r
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
png("./NYC_subway_plots/all_stations_z10.png",width=640, height =640, res = 75)
PlotOnStaticMap(gmap3, lat = lats , lon = lons, cex=1.5,pch=20, col = "blue")
```

```
## Error in LatLon2XY(lat, lon, zoom): object 'lats' not found
```

```r
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
png("./NYC_subway_plots/all_staph_z11.png",width=640, height =640, res = 75)
PlotOnStaticMap(gmap2, lat = lats , lon = lons, cex=1.5,pch=20, col = cols_staph)
```

```
## Error in LatLon2XY(lat, lon, zoom): object 'lats' not found
```

```r
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
png("./NYC_subway_plots/all_staph_z10.png",width=640, height =640, res = 75)
PlotOnStaticMap(gmap3, lat = lats , lon = lons, cex=1.5,pch=20, col = cols_staph)
```

```
## Error in LatLon2XY(lat, lon, zoom): object 'lats' not found
```

```r
dev.off()
```

```
## quartz_off_screen 
##                 2
```

### PLot each CC type


```r
for (i in colnames(staph_df)[1:33]){
    plot_CC_types(CC = i, mat = staph_df, map10 = gmap3, map11 = gmap2, plotdir = "./NYC_subway_plots/", SRA_file = staph_df_coords)

}
```

```
## Error in which(SRA_file[[CC]] > cutoff): object 'staph_df_coords' not found
```
### Mantel test for spatial autocorrelation of Staph aureus genotypes



```r
##First get Bray curtis matrix of dissimilarities
# braymat <- vegdist(staph_mat)
# jacmat <- vegdist(staph_mat, method = "jaccard")
# ##now get geographical distance
# 
# points_staph <- inner_join(NYCdata_SRA, Staph_betas, by = "Run") %>% select(Logitude,Latitude) 
# points_staph <- points_staph[-32,]#row 32 is mislabeled
# geogdist_staph <- dist(points_staph)
# set.seed(45678)
# mantel.rtest(geogdist_staph,braymat,nrepet = 9999)
# 
# mantel.rtest(geogdist_staph,jacmat,nrepet = 9999)
# 
```

### calculate hamming distances as an alternative using e1071 function


```r
binmat <- bintr(coords_df_subtype_mat,0.2)
```

```
## Error in mat[which(mat > cutoff)] <- 1: object 'coords_df_subtype_mat' not found
```

```r
#hammmingmat <- hamming.distance(binmat)
bin_df <- as.data.frame(cbind(staph_df_coords$Run,binmat))
```

```
## Error in cbind(staph_df_coords$Run, binmat): object 'staph_df_coords' not found
```

```r
colnames(bin_df)[1] <- "Run"
```

```
## Error in colnames(bin_df)[1] <- "Run": object 'bin_df' not found
```

```r
# mantel.rtest(geogdist_staph,hammmingmat,nrepet = 9999)
run_geog <- select(staph_df_coords,Run, Latitude, Logitude)
```

```
## Error in select_(.data, .dots = lazyeval::lazy_dots(...)): object 'staph_df_coords' not found
```

```r
combs <- combinations(r=2,v=run_geog$Run, n = length(run_geog$Run))
```

```
## Error in mode(n): object 'run_geog' not found
```

```r
geog_distance_vector <- sapply(1:nrow(combs), function(x) dist_between_stations(combs[x,],run_geog))
```

```
## Error in nrow(combs): object 'combs' not found
```

```r
hamm_distance_vector <- sapply(1:nrow(combs), function(x) H_distance_between_stations(combs[x,],bin_df))
```

```
## Error in nrow(combs): object 'combs' not found
```

```r
#some were taken from the same location - need to filter these out
zero_stations <- which(geog_distance_vector == 0)
```

```
## Error in which(geog_distance_vector == 0): object 'geog_distance_vector' not found
```

```r
geog_distance_vector <- geog_distance_vector[-(zero_stations)]
```

```
## Error in eval(expr, envir, enclos): object 'geog_distance_vector' not found
```

```r
hamm_distance_vector <- hamm_distance_vector[-(zero_stations)]
```

```
## Error in eval(expr, envir, enclos): object 'hamm_distance_vector' not found
```

```r
hist(geog_distance_vector, breaks = 50)
```

```
## Error in hist(geog_distance_vector, breaks = 50): object 'geog_distance_vector' not found
```

```r
hist(hamm_distance_vector)
```

```
## Error in hist(hamm_distance_vector): object 'hamm_distance_vector' not found
```

```r
boxplot(geog_distance_vector ~ hamm_distance_vector, xlab = "Hamming distance", ylab = "km")
```

```
## Error in eval(expr, envir, enclos): object 'geog_distance_vector' not found
```

```r
#try regression
reg <- lm(geog_distance_vector ~ as.numeric(hamm_distance_vector))
```

```
## Error in eval(expr, envir, enclos): object 'geog_distance_vector' not found
```

```r
segments(x0=1,x1=7,y0=reg$coefficients[1],y1=reg$coefficients[1]+(reg$coefficients[2]*6), col = "red")
```

```
## Error in segments(x0 = 1, x1 = 7, y0 = reg$coefficients[1], y1 = reg$coefficients[1] + : object 'reg' not found
```

```r
summary(reg)
```

```
## Error in summary(reg): object 'reg' not found
```

#functions for looking at geog distance of individual subtypes by permutation

```r
CC_geog_perm_test(SRA_file = staph_df_coords, CC = "CC_8", cutoff = 0.2, s = 345, reps= 1000)
```

```
## Error in which(SRA_file[[CC]] > cutoff): object 'staph_df_coords' not found
```

```r
CC_geog_perm_test(SRA_file = staph_df_coords, CC = "CC_30", cutoff = 0.2, s = 4564, reps= 1000)
```

```
## Error in which(SRA_file[[CC]] > cutoff): object 'staph_df_coords' not found
```

```r
CC_geog_perm_test(SRA_file = staph_df_coords, CC = "CC_5", cutoff = 0.2, s = 23, reps= 1000)
```

```
## Error in which(SRA_file[[CC]] > cutoff): object 'staph_df_coords' not found
```

```r
CC_geog_perm_test(SRA_file = staph_df_coords, CC = "CC_45", cutoff = 0.2, s = 455, reps= 1000)
```

```
## Error in which(SRA_file[[CC]] > cutoff): object 'staph_df_coords' not found
```

```r
CC_geog_perm_test(SRA_file = staph_df_coords, CC = "CC_15", cutoff = 0.2, s = 989, reps= 1000)
```

```
## Error in which(SRA_file[[CC]] > cutoff): object 'staph_df_coords' not found
```

```r
CC_geog_perm_test(SRA_file = staph_df_coords, CC = "CC_1", cutoff = 0.2, s = 467856765, reps= 1000)
```

```
## Error in which(SRA_file[[CC]] > cutoff): object 'staph_df_coords' not found
```



### session info


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
##  [1] assertthat_0.1      RColorBrewer_1.1-2  gtools_3.5.0       
##  [4] e1071_1.6-7         ade4_1.7-2          vegan_2.3-0        
##  [7] lattice_0.20-33     permute_0.8-4       dplyr_0.4.2        
## [10] RgoogleMaps_1.2.0.7 knitr_1.11         
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.0     cluster_2.0.3   magrittr_1.5    MASS_7.3-44    
##  [5] R6_2.1.1        stringr_1.0.0   tools_3.2.1     parallel_3.2.1 
##  [9] grid_3.2.1      nlme_3.1-122    mgcv_1.8-7      png_0.1-7      
## [13] DBI_0.3.1       class_7.3-14    lazyeval_0.1.10 RJSONIO_1.3-0  
## [17] Matrix_1.2-2    formatR_1.2     evaluate_0.7.2  stringi_0.5-5  
## [21] methods_3.2.1
```


