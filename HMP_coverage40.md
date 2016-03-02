---
title: "Creating Combined starting file using all 40 original CCs"
author: "Tim Read"
date: "March 1, 2016"
output: html_document
---

The original "combined", "cov0.025" and "cov0.5", were made with a step that merged all the CC5, CC8 and CC30 subgroups together.  Prompted by a reviewer's comment, this script remakes the files, which are used in downstream processing.


```r
devtools::install_github("jennybc/googlesheets")
```

```
## Skipping install for github remote, the SHA1 (debf4a81) has not changed since last install.
##   Use `force = TRUE` to force installation
```


```r
library(devtools)
library(googlesheets)
library()
library("dplyr")
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(assertthat)
library(gtools)
library(readr)
source('./staph_metagenome_tools.R')
```

Downlad data from google sheet

```r
covs <- gs_title("Coverage_Staph_MeCA") # note: might need to register app with google here .
```

```
## Auto-refreshing stale OAuth token.
```

```
## Sheet successfully identified: "Coverage_Staph_MeCA"
```

```r
print(covs$updated)
```

```
## [1] "2015-08-05 17:41:50 GMT"
```

```r
mapping <-gs_read(covs)
```

```
## Accessing worksheet titled 'anterior_nares_MecA_Staph_Posit'.
```

```
## No encoding supplied: defaulting to UTF-8.
```

Get the binstrain output


```r
dat3 <- read.csv("./Data/Final_HMP_Matrix.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(dat3)[colnames(dat3) == "CC_8_72"] <- "CC_72"
dat3 <- dat3[,-(2:5)] #drop redundant cols
names(dat3)[names(dat3) == "MLST_93"] <- "CC_93"
dat3 <- dat3[,mixedsort(colnames(dat3))]
#colnames(dat3) <- gsub("(CC_.{0,3})_.{0,4}$","\\1",colnames(dat3))
#remove missing samples
missing_samples <- c("SRS021960","SRS098620") #data consistency issue with these samples
dat4 <- filter(dat3,!(Sample.Id %in% missing_samples))
combined40 <- inner_join(mapping,dat4,by = c("Sample_id" = "Sample.Id"))

#remove duplicated records tidy up col names
combined40 <- combined40[-c(46,95,96),]
colnames(combined40)[which(colnames(combined40) == "Body_Site")] <-  "Body.site"
colnames(combined40)[which(colnames(combined40) == "Sample_id")]  <- "Sample.Id"
```
Now filter by cverage and save

```r
cov40_0.025 <- filter(combined40, Staph_cov > 0.025)
cov40_0.5 <- filter(combined40, Staph_cov > 0.5)
write.table(cov40_0.025,"./Data/cov40_0.025")
write.table(cov40_0.5,"./Data/cov40_0.5")
write.table(combined40,"./Data/combined40")
```


