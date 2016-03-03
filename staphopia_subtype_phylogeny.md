---
title: "Phylogenetic tree of subtype tested strains"
author: "Tim Read"
date: "October 6, 2015"
output: html_document
---


```r
source("~/.staphopia_logon.R")
db <- staphopia_logon()
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

```
## Error: RPostgreSQL package required to connect to postgres db
```

```r
library(staphopiaRtools) 
```

```
## Error in library(staphopiaRtools): there is no package called 'staphopiaRtools'
```

```r
library(assertthat)
library(IRanges)
```

```
## Loading required package: methods
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:dplyr':
## 
##     combine, intersect, setdiff, union
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
##     grep, grepl, intersect, is.unsorted, lapply, lengths, Map,
##     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
##     pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
##     setdiff, sort, table, tapply, union, unique, unlist, unsplit
```

```
## Loading required package: S4Vectors
```

```
## Loading required package: stats4
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:dplyr':
## 
##     rename
```

```
## 
## Attaching package: 'IRanges'
```

```
## The following objects are masked from 'package:dplyr':
## 
##     collapse, desc, slice
```

```r
library(Biostrings)
```

```
## Loading required package: XVector
```

```r
library(phangorn)
```

```
## Loading required package: ape
```
test

```r
# P <- pull_ids(db,project = "PRJNA239001") %>% filter(st_stripped == 105)
# staphopia_alignment_pipeline(db,P,fasta_file = "./aligned_105_CDC.fasta", ref = "N315")
```


pull list of strains used

```r
strain_tag_df <- read.csv("./Data/strains_used_for_subtype_tests.csv",stringsAsFactors = FALSE, header = TRUE, col.names = c("tags","st"))
```
download data


```r
sample_tab <- tbl(db,"sample_metadata") %>% select(id,sample_tag)
selected_samples_tab <- filter(sample_tab,sample_tag %in% strain_tag_df$tags)
staphopia_alignment_pipeline(db,selected_samples_tab,fasta_file = "~/aligned_subtyping_strains.fasta", ref = "N315") #store locally b/c of file size
```

```
## Error in postgresqlExecStatement(conn, statement, ...): RS-DBI driver: (could not Retrieve the result : FATAL:  terminating connection due to administrator command
## SSL connection has been closed unexpectedly
## )
```

make phyDat file and becasue it is enormous cache it


```r
staph_align <- read.phyDat("~/aligned_subtyping_strains.fasta", format = "fasta", type = "DNA")
```

```
## Warning in file(con, "rb"): cannot open file '/Users/timread/
## aligned_subtyping_strains.fasta': No such file or directory
```

```
## Error in file(con, "rb"): cannot open the connection
```

```r
save(staph_align,file = "~/staph_align")
```

```
## Error in save(staph_align, file = "~/staph_align"): object 'staph_align' not found
```

distance matrix - this takes about 4 hours

```r
dm = dist.dna(as.DNAbin(staph_align))
```

```
## Error in as.DNAbin(staph_align): object 'staph_align' not found
```

```r
save(dm,file="~/dm")
```

```
## Error in save(dm, file = "~/dm"): object 'dm' not found
```
and NJ tree


```r
treeNJ = NJ(dm)
```

```
## Error in nj(x): object 'dm' not found
```

```r
layout(matrix(c(1,1), 2, 1), height=c(1,2))
par(mar = c(.1,.1,.1,.1))
plot(treeNJ, "unrooted", show.tip.label = FALSE)
```

```
## Error in plot(treeNJ, "unrooted", show.tip.label = FALSE): error in evaluating the argument 'x' in selecting a method for function 'plot': Error: object 'treeNJ' not found
```

```r
sessionInfo()
```

```
## R version 3.2.3 (2015-12-10)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: OS X 10.11.3 (El Capitan)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats4    parallel  methods   stats     graphics  grDevices utils    
## [8] datasets  base     
## 
## other attached packages:
##  [1] phangorn_2.0.2      ape_3.4             Biostrings_2.38.3  
##  [4] XVector_0.10.0      IRanges_2.4.6       S4Vectors_0.8.5    
##  [7] BiocGenerics_0.16.1 assertthat_0.1      dplyr_0.4.3        
## [10] knitr_1.12.3       
## 
## loaded via a namespace (and not attached):
##  [1] igraph_1.0.1    Rcpp_0.12.3     magrittr_1.5    zlibbioc_1.16.0
##  [5] nnls_1.4        lattice_0.20-33 R6_2.1.2        quadprog_1.5-5 
##  [9] stringr_1.0.0   tools_3.2.3     grid_3.2.3      nlme_3.1-125   
## [13] DBI_0.3.1       digest_0.6.9    Matrix_1.2-3    formatR_1.2.1  
## [17] evaluate_0.8    stringi_1.0-1
```

