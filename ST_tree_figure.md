---
title: "STs plotted on subtyping tree"
author: "Tim Read"
date: "March 2016"
output: html_document
---



```r
library(phangorn)
```

```
## Loading required package: ape
```

```r
library(dplyr)
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
source('./staph_metagenome_tools.R')
```




```r
STstrains <- read.csv("./Data/S2_Data.csv", header = TRUE, stringsAsFactors = FALSE)
```

Load tree


```r
load("~/dm")
NJ <- nj(dm)
```
 Find and label ST5 tips
 


```r
plotST <- function(ST){
  strains <- filter(STstrains, MLST.type == ST) %>% select(SRA.id)
  tips2label <- which(NJ$tip.label %in% strains[,"SRA.id"])
  tit <- paste("ST",ST," strains on tree",sep = "")
  plot(NJ, "unrooted", show.tip.label = FALSE, main = tit)
  tiplabels(tip = tips2label, pch= 20, col = "red")
}
```


```r
MLST_vec <- select(STstrains,MLST.type) %>% unique() 
```



```r
for(i in MLST_vec[,"MLST.type"]){
  plotST(i)
}
```

![plot of chunk ST_plots](figure/ST_plots-1.png)![plot of chunk ST_plots](figure/ST_plots-2.png)![plot of chunk ST_plots](figure/ST_plots-3.png)![plot of chunk ST_plots](figure/ST_plots-4.png)![plot of chunk ST_plots](figure/ST_plots-5.png)![plot of chunk ST_plots](figure/ST_plots-6.png)![plot of chunk ST_plots](figure/ST_plots-7.png)![plot of chunk ST_plots](figure/ST_plots-8.png)![plot of chunk ST_plots](figure/ST_plots-9.png)![plot of chunk ST_plots](figure/ST_plots-10.png)![plot of chunk ST_plots](figure/ST_plots-11.png)![plot of chunk ST_plots](figure/ST_plots-12.png)![plot of chunk ST_plots](figure/ST_plots-13.png)![plot of chunk ST_plots](figure/ST_plots-14.png)![plot of chunk ST_plots](figure/ST_plots-15.png)![plot of chunk ST_plots](figure/ST_plots-16.png)![plot of chunk ST_plots](figure/ST_plots-17.png)![plot of chunk ST_plots](figure/ST_plots-18.png)![plot of chunk ST_plots](figure/ST_plots-19.png)![plot of chunk ST_plots](figure/ST_plots-20.png)![plot of chunk ST_plots](figure/ST_plots-21.png)![plot of chunk ST_plots](figure/ST_plots-22.png)![plot of chunk ST_plots](figure/ST_plots-23.png)![plot of chunk ST_plots](figure/ST_plots-24.png)![plot of chunk ST_plots](figure/ST_plots-25.png)![plot of chunk ST_plots](figure/ST_plots-26.png)![plot of chunk ST_plots](figure/ST_plots-27.png)![plot of chunk ST_plots](figure/ST_plots-28.png)![plot of chunk ST_plots](figure/ST_plots-29.png)![plot of chunk ST_plots](figure/ST_plots-30.png)![plot of chunk ST_plots](figure/ST_plots-31.png)![plot of chunk ST_plots](figure/ST_plots-32.png)![plot of chunk ST_plots](figure/ST_plots-33.png)![plot of chunk ST_plots](figure/ST_plots-34.png)![plot of chunk ST_plots](figure/ST_plots-35.png)![plot of chunk ST_plots](figure/ST_plots-36.png)![plot of chunk ST_plots](figure/ST_plots-37.png)![plot of chunk ST_plots](figure/ST_plots-38.png)![plot of chunk ST_plots](figure/ST_plots-39.png)![plot of chunk ST_plots](figure/ST_plots-40.png)![plot of chunk ST_plots](figure/ST_plots-41.png)![plot of chunk ST_plots](figure/ST_plots-42.png)![plot of chunk ST_plots](figure/ST_plots-43.png)![plot of chunk ST_plots](figure/ST_plots-44.png)![plot of chunk ST_plots](figure/ST_plots-45.png)![plot of chunk ST_plots](figure/ST_plots-46.png)![plot of chunk ST_plots](figure/ST_plots-47.png)![plot of chunk ST_plots](figure/ST_plots-48.png)![plot of chunk ST_plots](figure/ST_plots-49.png)![plot of chunk ST_plots](figure/ST_plots-50.png)![plot of chunk ST_plots](figure/ST_plots-51.png)![plot of chunk ST_plots](figure/ST_plots-52.png)![plot of chunk ST_plots](figure/ST_plots-53.png)![plot of chunk ST_plots](figure/ST_plots-54.png)![plot of chunk ST_plots](figure/ST_plots-55.png)![plot of chunk ST_plots](figure/ST_plots-56.png)

Load in the binstrain assignments for the ST5 and dive into those that were assigned to ST5, those that had a strong evidecne for assingment to another strain (binstrain > 0.8), and thse that had little evidence fr either


```r
St5_assign <- read.table("./Data/Final_highest_Beta_real_ST_5.txt",header = FALSE, stringsAsFactors = FALSE)
CC5_St5 <- filter(St5_assign, grepl("CC_5",V2)) %>% select(V1)
ST5_not_CC5 <- filter(St5_assign, !grepl("CC_5",V2)) %>% filter(V3 > 0.8) %>% select(V1)
ST5_not_anything <- filter(St5_assign, !grepl("CC_5",V2)) %>% filter(V3 <= 0.8) %>% select(V1)
```
### Now print the different assignments

```r
tips2label <- which(NJ$tip.label %in% CC5_St5[,"V1"])
tit <- "ST5 strains binstrain type CC5"
plot(NJ, "unrooted", show.tip.label = FALSE, main = tit)
tiplabels(tip = tips2label, pch= 20, col = "blue")
```

![plot of chunk ST5_classed_as_CC5](figure/ST5_classed_as_CC5-1.png)


```r
tips2label <- which(NJ$tip.label %in% ST5_not_CC5[,"V1"])
tit <- "ST5 strains binstrain type another CC low score"
plot(NJ, "unrooted", show.tip.label = FALSE, main = tit)
tiplabels(tip = tips2label, pch= 20, col = "blue")
```

![plot of chunk ST5_classed_as_not_CC5](figure/ST5_classed_as_not_CC5-1.png)


```r
tips2label <- which(NJ$tip.label %in% ST5_not_anything[,"V1"])
tit <- "ST5 strains binstrain type another CC"
plot(NJ, "unrooted", show.tip.label = FALSE, main = tit)
tiplabels(tip = tips2label, pch= 20, col = "blue")
```

![plot of chunk ST5](figure/ST5-1.png)


