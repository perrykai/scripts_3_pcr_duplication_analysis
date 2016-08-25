**Script:** `3_24illumina_create_dge_normalize_mirna_exp.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/2_24illumina_mirna_expression_characterization/scripts`

**Date:**  8/18/16

**Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/2_24illumina_mirna_expression_characterization/`

**Input File(s):** 
1. `1_24illumina_filtered_rounded_mean_mature_mirna_exp.Rdata`
2. `2_24illumina_mature_mirna_annotation.Rdata`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/2_24illumina_mirna_expression_characterization/`

**Output File(s):** `4_24illumina_dge_normalized_mirna_expression.Rdata`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Save data](#save-data)

## Objectives

The objective of this script is to create a dge object of the subset of 24 Illumina libraries. 
Additionally, the read counts of the miRNAs will be normalized using the cpm function of edgeR, filtered for expression (rough filter: < 1cpm in > 6 libraries removed) 
and the calcNormFactors and estimateCommonDisp functions of edgeR will be applied to the non-cpm read counts prior to the final cpm normalization of the read counts. 

THIS ANALYSIS COMPLETED USING R/3.2.0
## Install libraries


```r
library(limma)
```

```
## Loading required package: methods
```

```r
library(edgeR)

rm(list=ls())

setwd("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/2_24illumina_mirna_expression_characterization/scripts/")
```

## Load data
Load the miRDeep2 read count data for the 24 Illumina libraries


```r
load("../1_24illumina_filtered_rounded_mean_mature_mirna_exp.Rdata")
```

Load the annotation file from the mature miRNAs


```r
load("../2_24illumina_mature_mirna_annotation.Rdata")
```


## Analysis

### Create the dge object and transform the read counts: log-cpm, then filter genes by expression (1 cpm in 44 or more samples retained):

Create the dge object:


```r
dge<-DGEList(counts=illumina24.no.zero.dfmeanrcround,genes=illumina24.total.mature.annot2)
dim(dge)
```

```
## [1] 325  24
```

```r
dge[1:5,1:5]
```

```
## An object of class "DGEList"
## $counts
##                1034  1058 1080  1096  1116
## ssc-let-7a    24174 53646 6014 57852 56632
## ssc-let-7c     3001  5920 1765  8426  7084
## ssc-let-7d-3p   211   232 2040   643   661
## ssc-let-7d-5p  1246  3067  364  3617  3382
## ssc-let-7e     1026  2068  270  2582  2601
## 
## $samples
##      group lib.size norm.factors
## 1034     1  1837116            1
## 1058     1  3712023            1
## 1080     1  1838011            1
## 1096     1  3533423            1
## 1116     1  3390992            1
## 
## $genes
##                        Name  chr0     start       end width strand  type
## ssc-let-7a       ssc-let-7a  chr3  44864443  44864464    22      + miRNA
## ssc-let-7c       ssc-let-7c chr13 191559351 191559372    22      + miRNA
## ssc-let-7d-3p ssc-let-7d-3p  chr3  44867331  44867352    22      + miRNA
## ssc-let-7d-5p ssc-let-7d-5p  chr3  44867277  44867298    22      + miRNA
## ssc-let-7e       ssc-let-7e  chr6  51858385  51858406    22      + miRNA
##                      Alias          Precursors
## ssc-let-7a    MIMAT0013865 MI0017984,MI0013085
## ssc-let-7c    MIMAT0002151           MI0002445
## ssc-let-7d-3p MIMAT0025357           MI0022120
## ssc-let-7d-5p MIMAT0025356           MI0022120
## ssc-let-7e    MIMAT0013866           MI0013086
```

Calculate the read counts per million in order to filter miRNAs by normalized expression:


```r
cpm.dge<-cpm(dge)
dim(cpm.dge)
```

```
## [1] 325  24
```

```r
cpm.dge[1:5,1:5]
```

```
##                     1034        1058      1080       1096       1116
## ssc-let-7a    13158.6683 14451.95787 3272.0152 16372.7920 16700.7177
## ssc-let-7c     1633.5387  1594.81770  960.2772  2384.6565  2089.0642
## ssc-let-7d-3p   114.8539    62.49961 1109.8954   181.9765   194.9282
## ssc-let-7d-5p   678.2370   826.23410  198.0402  1023.6533   997.3483
## ssc-let-7e      558.4841   557.10862  146.8979   730.7362   767.0322
```

```r
if (sum(rownames(illumina24.no.zero.dfmeanrcround)!=rownames(cpm.dge))!=0) stop ("miRNAs not the same between read counts and cpm")
if (sum(colnames(illumina24.no.zero.dfmeanrcround)!=colnames(cpm.dge))!=0) stop ("animal ids not the same between read counts and cpm")
if (sum(rownames(illumina24.no.zero.dfmeanrcround)!=rownames(illumina24.total.mature.annot2))!=0) stop ("miRNAs not the same between read counts and annotation")
```

Filter miRNAs with at least 1 cpm in at least 1/4 of the samples (24/4=6)


```r
filtercpm<-rowSums(cpm.dge>=1)>=6
sum(filtercpm)
```

```
## [1] 266
```

```r
nrow(cpm.dge) - sum(filtercpm)
```

```
## [1] 59
```

We are removing 59 miRNA profiles from the analysis

So, keep the miRNA profiles in dge based on those retained in the cpm-filtering step:

This retains the rounded, filtered mean read counts, not the cpm (this will be done later):


```r
dge<-dge[filtercpm,]
names(dge)
```

```
## [1] "counts"  "samples" "genes"
```

```r
dge[1:5,1:5]
```

```
## An object of class "DGEList"
## $counts
##                1034  1058 1080  1096  1116
## ssc-let-7a    24174 53646 6014 57852 56632
## ssc-let-7c     3001  5920 1765  8426  7084
## ssc-let-7d-3p   211   232 2040   643   661
## ssc-let-7d-5p  1246  3067  364  3617  3382
## ssc-let-7e     1026  2068  270  2582  2601
## 
## $samples
##      group lib.size norm.factors
## 1034     1  1837116            1
## 1058     1  3712023            1
## 1080     1  1838011            1
## 1096     1  3533423            1
## 1116     1  3390992            1
## 
## $genes
##                        Name  chr0     start       end width strand  type
## ssc-let-7a       ssc-let-7a  chr3  44864443  44864464    22      + miRNA
## ssc-let-7c       ssc-let-7c chr13 191559351 191559372    22      + miRNA
## ssc-let-7d-3p ssc-let-7d-3p  chr3  44867331  44867352    22      + miRNA
## ssc-let-7d-5p ssc-let-7d-5p  chr3  44867277  44867298    22      + miRNA
## ssc-let-7e       ssc-let-7e  chr6  51858385  51858406    22      + miRNA
##                      Alias          Precursors
## ssc-let-7a    MIMAT0013865 MI0017984,MI0013085
## ssc-let-7c    MIMAT0002151           MI0002445
## ssc-let-7d-3p MIMAT0025357           MI0022120
## ssc-let-7d-5p MIMAT0025356           MI0022120
## ssc-let-7e    MIMAT0013866           MI0013086
```

```r
dim(dge$counts)
```

```
## [1] 266  24
```

```r
if (sum(colnames(dge)!=colnames(cpm.dge))!=0) stop ("colnames not the same between dge and cpm.dge")
```

Apply the TMM normalization:


```r
dge<-calcNormFactors(dge)
head(dge$samples)
```

```
##      group lib.size norm.factors
## 1034     1  1837116    1.0997476
## 1058     1  3712023    1.0954787
## 1080     1  1838011    0.5287203
## 1096     1  3533423    1.2071338
## 1116     1  3390992    1.0525756
## 1134     1  4031993    1.0476428
```

```r
hist(dge$samples$norm.factors)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png) 

This function (estimateCommonDisp) applies normalization factors, caluclates normalized expression based on robust count of normalized reads.


```r
dge<-estimateCommonDisp(dge,verbose=TRUE)
```

```
## Disp = 0.2614 , BCV = 0.5113
```

```r
dge$common.dispersion
```

```
## [1] 0.2613973
```

The cpm conversion is used to normalize the read counts between libraries (accounts for differences in sequencing depth)


```r
illumina24.dge.filter.cpm<-cpm(dge)
head(illumina24.dge.filter.cpm)
```

```
##                     1034        1058      1080       1096       1116
## ssc-let-7a    11965.1710 13192.36804 6188.5556 13563.3610 15866.5259
## ssc-let-7c     1485.3759  1455.81812 1816.2289  1975.4698  1984.7166
## ssc-let-7d-3p   104.4366    57.05233 2099.2107   150.7509   185.1917
## ssc-let-7d-5p   616.7206   754.22199  374.5651   848.0031   947.5313
## ssc-let-7e      507.8293   508.55268  277.8367   605.3481   728.7193
## ssc-let-7f     6875.9889  9344.28664 2190.7940  4599.8953  7031.6794
##                     1134       1154       1170        1194      1240
## ssc-let-7a    14521.0026 14468.0583 13908.5145 12906.41460 8236.7019
## ssc-let-7c     1653.6112  2010.1989  1535.3617  1449.76760 1534.9533
## ssc-let-7d-3p   116.7116   171.4653    81.1370    79.81704  356.8581
## ssc-let-7d-5p   764.1885   801.1485   757.3587   677.71919  526.4815
## ssc-let-7e      562.7250   597.4418   574.9205   399.81079  325.3434
## ssc-let-7f    11254.2623  4551.4029  9861.7465  9818.94656 5940.0652
##                     1278       1300        1426        1434        1458
## ssc-let-7a    11767.0266 12706.5273 15757.24079 11191.30202 13796.55136
## ssc-let-7c     1819.8801  1578.0898  1717.41269  1492.87561  1565.46844
## ssc-let-7d-3p   154.5153   139.9004    69.46328    83.57906    58.10252
## ssc-let-7d-5p   630.8878   756.9056   945.27779   595.43084   744.85998
## ssc-let-7e      468.0203   449.9498   620.63485   434.63145   513.31067
## ssc-let-7f     7860.8916  7168.9480 12038.37861 10355.25709  9659.36490
##                     1484       1502       1512        1534       1580
## ssc-let-7a    16654.5466 19606.3663 11612.4553 17922.39657 18699.2042
## ssc-let-7c     2161.7949  1915.0313  1999.3701  1864.80290  2291.5400
## ssc-let-7d-3p   173.8961    93.7773   457.9916    80.74084   124.4318
## ssc-let-7d-5p   865.4073  1228.9762   813.3299  1040.31469  1188.0388
## ssc-let-7e      672.7116   772.9224   554.3277   715.28105   744.9286
## ssc-let-7f     5935.3403 13359.8108  5854.3957 13539.10049  9171.8593
##                      1594        1640       1644        1662
## ssc-let-7a    15838.16227 18479.50948 20034.3404 14957.98372
## ssc-let-7c     1728.98945  1684.28365  2122.0257  1591.11549
## ssc-let-7d-3p    79.14388    93.03745   114.1596    76.46593
## ssc-let-7d-5p   940.84828  1234.82103  1205.3912   804.60881
## ssc-let-7e      638.22439   734.69012  1009.3053   540.25518
## ssc-let-7f    13187.85685 15681.83396 13773.0214 11365.02147
```

```r
dim(illumina24.dge.filter.cpm)
```

```
## [1] 266  24
```

## Save data


```r
save(illumina24.dge.filter.cpm, file="../4_24illumina_dge_normalized_mirna_expression.Rdata")
```

