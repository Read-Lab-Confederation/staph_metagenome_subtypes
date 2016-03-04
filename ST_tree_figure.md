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
source('./staph_metagenome_tools.R')
```




```r
STstrains <- read.csv("./Data/S2_Data.csv", header = TRUE, stringsAsFactors = FALSE)
```

Load tree


```r
load("~/dm")
```

```
## Warning in readChar(con, 5L, useBytes = TRUE): cannot open compressed file
## '/Users/timothyread/dm', probable reason 'No such file or directory'
```

```
## Error in readChar(con, 5L, useBytes = TRUE): cannot open the connection
```

```r
NJ <- nj(dm)
```

```
## Error in nj(dm): object 'dm' not found
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

```
## Error in NJ$tip.label: object of type 'closure' is not subsettable
```



