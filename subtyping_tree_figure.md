---
title: "Subtyping tree figure"
author: "Tim Read"
date: "October 13, 2015"
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
strains <- read.csv("./Data/2114_strain_subtypes.csv", header = TRUE, stringsAsFactors = FALSE)
CCs <- select(strains, Reference.CC) %>% unique() %>% arrange()
```

Load tree


```r
load("~/dm")
NJ <- nj(dm)
plot(NJ, "unrooted", show.tip.label = FALSE)
```

![plot of chunk load_tree](figure/load_tree-1.png)
 Find and label CC_30 lables
 
 Run major groups with beta cutoff of .65

```r
decorate_staph_tree("CC_30",NJ,strains)
```

![plot of chunk plots_trees0.65](figure/plots_trees0.65-1.png)

```r
decorate_staph_tree("CC_5_5",NJ,strains)
```

![plot of chunk plots_trees0.65](figure/plots_trees0.65-2.png)

```r
decorate_staph_tree("CC_8_",NJ,strains)
```

![plot of chunk plots_trees0.65](figure/plots_trees0.65-3.png)

```r
for (i in CCs$Reference.CC){
  decorate_staph_tree(i,NJ,strains)
}
```

![plot of chunk plots_trees0.65](figure/plots_trees0.65-4.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-5.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-6.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-7.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-8.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-9.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-10.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-11.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-12.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-13.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-14.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-15.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-16.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-17.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-18.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-19.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-20.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-21.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-22.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-23.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-24.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-25.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-26.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-27.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-28.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-29.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-30.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-31.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-32.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-33.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-34.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-35.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-36.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-37.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-38.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-39.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-40.png)![plot of chunk plots_trees0.65](figure/plots_trees0.65-41.png)
Run major groups with beta cutoff of .80

```r
decorate_staph_tree("CC_30",NJ,strains, cutoff = 0.8, deco = "blue")
```

![plot of chunk plots_trees0.85](figure/plots_trees0.85-1.png)

```r
decorate_staph_tree("CC_5_5",NJ,strains)
```

![plot of chunk plots_trees0.85](figure/plots_trees0.85-2.png)

```r
decorate_staph_tree("CC_8_",NJ,strains)
```

![plot of chunk plots_trees0.85](figure/plots_trees0.85-3.png)

```r
for (i in CCs$Reference.CC){
  decorate_staph_tree(i,NJ,strains,cutoff = 0.8, deco = "blue")
}
```

![plot of chunk plots_trees0.85](figure/plots_trees0.85-4.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-5.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-6.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-7.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-8.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-9.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-10.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-11.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-12.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-13.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-14.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-15.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-16.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-17.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-18.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-19.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-20.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-21.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-22.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-23.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-24.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-25.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-26.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-27.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-28.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-29.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-30.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-31.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-32.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-33.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-34.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-35.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-36.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-37.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-38.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-39.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-40.png)![plot of chunk plots_trees0.85](figure/plots_trees0.85-41.png)
Decorated tree figure for grant


```r
plot(NJ, "unrooted", show.tip.label = FALSE)
#cc30
tl <- filter(strains,grepl("CC_30",Reference.CC)) %>% filter(Beta > 0.80) %>% select(Sample.Id.of.0.75X)
tps <- which(NJ$tip.label %in% tl$Sample.Id.of.0.75X)
tiplabels(tip = tps, pch= 20, col = "red")
##CC_5_5
tl <- filter(strains,grepl("CC_5_5",Reference.CC)) %>% filter(Beta > 0.80) %>% select(Sample.Id.of.0.75X)
tps <- which(NJ$tip.label %in% tl$Sample.Id.of.0.75X)
tiplabels(tip = tps, pch= 20, col = "green")
#st8
tl <- filter(strains,grepl("CC_8_8_2",Reference.CC)) %>% filter(Beta > 0.80) %>% select(Sample.Id.of.0.75X)
tps <- which(NJ$tip.label %in% tl$Sample.Id.of.0.75X)
tiplabels(tip = tps, pch= 20, col = "blue")
#st75 argentius
tl <- filter(strains,grepl("CC_75",Reference.CC)) %>% filter(Beta > 0.80) %>% select(Sample.Id.of.0.75X)
tps <- which(NJ$tip.label %in% tl$Sample.Id.of.0.75X)
tiplabels(tip = tps, pch= 20, col = "gray")
#cc133
tl <- filter(strains,grepl("CC_133",Reference.CC)) %>% filter(Beta > 0.80) %>% select(Sample.Id.of.0.75X)
tps <- which(NJ$tip.label %in% tl$Sample.Id.of.0.75X)
tiplabels(tip = tps, pch= 20, col = "orange")
add.scale.bar()
```

![plot of chunk pretty_tree](figure/pretty_tree-1.png)
Check for senstivity


```r
tl <- filter(strains,Beta > 0.80) %>% select(Sample.Id.of.0.75X)
plot(NJ, "unrooted", show.tip.label = FALSE, main = "Sensitivity: beta > 0.65")
tps <- which(NJ$tip.label %in% tl$Sample.Id.of.0.75X)
tiplabels(tip = tps, pch= 20, col = "red")
```

![plot of chunk sensitivity_tree](figure/sensitivity_tree-1.png)

Look at SNP#1752540, which is common in ST398

```r
SNP <- read.table("./Data/SNP1752540_sample_tags.txt", header = FALSE, stringsAsFactors = FALSE)
SNPtps <- which(NJ$tip.label %in% SNP$V1)
plot(NJ, "unrooted", show.tip.label = FALSE, main = "SNP#1752540")
tiplabels(tip = SNPtps, pch= 20, col = "red")
```

![plot of chunk SNP#1752540](figure/SNP#1752540-1.png)
