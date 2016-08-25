**Script:** `1_create_24illumina_mirna_expression_df.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/2_24illumina_mirna_expression_characterization/scripts`

**Date:**  8/17/16

**Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/9_mirdeep2_core_quantify_predict_24illumina_output`

**Input File(s):** `miRNAs_expressed_all_samples_17_08_2016_t_11_33_33.csv`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/2_24illumina_mirna_expression_characterization`

**Output File(s):** `1_24illumina_filtered_rounded_mean_mature_mirna_exp.Rdata`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
6. [Save data](#save-data)

## Objectives

The objective of this script is to create a matrix of miRNA expression data (read counts) for incorporation into the PCR duplicate analysis, beginning with the known miRNA expression data output from the miRDeep2 core module.
To acheive one read count per miRNA per animal, I will need to take the average of the mature read counts in instances of multiple precursor sequences.

The result of this script will be a data frame containing one average read count per miRNA per animal.

So, what I need to do:

1. Extract the columns of the mature miRNA read counts for each animal
2. Use the 'by' function to apply the function colmeans to the entire data frame of read counts (What this will do is go down the columns looking at the index of grouped miRNA names and take the average of the read counts for that group of miRNA. The result of this will be a list containing the average read counts for each miRNA for each animal)
3. Transform the list output back into a data.frame using the plyr package to prepare for gblup function
4. Filter the data for expression threshold: The total read count for the miRNA needs to be greater than 0
5. Restore the pig IDs in place of the 3-digit codes as the column names of the data frame, for use with the gblup function of gwaR
6. Round the mean read counts to the nearest integer for use with the edgeR package.
## Install libraries


```r
library(plyr)


setwd("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/2_24illumina_mirna_expression_characterization/scripts")
rm(list=ls())
```

## Load data
### 1. Read in the raw read counts output from miRDeep2 quantifier module


```r
rawcounts<-read.csv("../../1_preprocess_24_illumina_samples/9_mirdeep2_core_quantify_predict_24illumina_output/miRNAs_expressed_all_samples_17_08_2016_t_11_33_33.csv", sep = "\t", header = TRUE, row.names=NULL)
head(rawcounts)
```

```
##         X.miRNA read_count    precursor   total  X001  X002 X003  X004
## 1    ssc-let-7a    1416064 ssc-let-7a-1 1416064 24170 53629 6007 57842
## 2    ssc-let-7a    1416269 ssc-let-7a-2 1416269 24179 53662 6022 57862
## 3    ssc-let-7c     174980   ssc-let-7c  174980  3001  5920 1765  8426
## 4 ssc-let-7d-5p      80709   ssc-let-7d   80709  1246  3067  364  3617
## 5 ssc-let-7d-3p      14107   ssc-let-7d   14107   211   232 2040   643
## 6    ssc-let-7e      56528   ssc-let-7e   56528  1026  2068  270  2582
##    X005  X006  X007  X008  X009  X010  X011  X012  X013   X014  X015
## 1 56616 61336 59220 57935 35572 35556 39448 36971 76453 220007 96167
## 2 56649 61339 59249 57945 35576 35534 39449 36962 76438 219991 96170
## 3  7084  6985  8230  6396  3996  6624  6101  4591  8332  29347 10912
## 4  3382  3228  3280  3155  1868  2272  2115  2202  4586  11705  5192
## 5   661   493   702   338   220  1540   518   407   337   1643   405
## 6  2601  2377  2446  2395  1102  1404  1569  1309  3011   8544  3578
##     X016  X017  X018  X019  X020   X021  X022  X023  X024 X001.norm.
## 1 106326 19857 14697 34623 78725 124880 42282 29813 47932   10087.68
## 2 106289 19866 14716 34633 78765 124868 42332 29854 47919   10091.43
## 3  13799  1940  2532  3603  9650  13632  3856  3160  5098    1252.51
## 4   5524  1245  1030  2010  5003   7418  2827  1795  2578     520.03
## 5   1110    95   580   156   524    624   213   170   245      88.06
## 6   4294   783   702  1382  3137   5032  1682  1503  1731     428.21
##   X002.norm. X003.norm. X004.norm. X005.norm. X006.norm. X007.norm.
## 1   11594.14    1810.71   11514.40   11947.50   11309.73   11797.42
## 2   11601.28    1815.23   11518.38   11954.47   11310.28   11803.19
## 3    1279.85     532.03    1677.33    1494.92    1287.96    1639.53
## 4     663.06     109.72     720.02     713.69     595.21     653.42
## 5      50.16     614.92     128.00     139.49      90.90     139.85
## 6     447.08      81.39     513.99     548.88     438.29     487.28
##   X008.norm. X009.norm. X010.norm. X011.norm. X012.norm. X013.norm.
## 1   12814.92    8919.38    5620.79    8990.69   10559.25   12958.81
## 2   12817.14    8920.39    5617.31    8990.92   10556.68   12956.27
## 3    1414.76    1001.96    1047.14    1390.49    1311.23    1412.28
## 4     697.87     468.39     359.16     482.04     628.91     777.33
## 5      74.76      55.16     243.45     118.06     116.24      57.12
## 6     529.76     276.32     221.95     357.59     373.86     510.37
##   X014.norm. X015.norm. X016.norm. X017.norm. X018.norm. X019.norm.
## 1   13372.06   10673.77   12318.71   13944.40    7432.16   12734.51
## 2   13371.09   10674.10   12314.42   13950.72    7441.77   12738.19
## 3    1783.72    1211.15    1598.72    1362.35    1280.41    1325.20
## 4     711.43     576.27     640.00     874.29     520.86     739.29
## 5      99.86      44.95     128.60      66.71     293.30      57.38
## 6     519.31     397.13     497.49     549.85     355.00     508.31
##   X020.norm. X021.norm. X022.norm. X023.norm. X024.norm.
## 1   12150.72   11132.34   12252.43   10955.78   10831.56
## 2   12156.89   11131.27   12266.92   10970.84   10828.63
## 3    1489.42    1215.22    1117.39    1161.25    1152.03
## 4     772.18     661.27     819.21     659.63     582.57
## 5      80.88      55.63      61.72      62.47      55.36
## 6     484.18     448.57     487.41     552.33     391.17
```

```r
dim(rawcounts)
```

```
## [1] 492  52
```

Set the name of the first column to 'miRNA'


```r
colnames(rawcounts)[[1]] <- "miRNA"
```

Remove the "X" character from the beginning of each column:


```r
colnames(rawcounts)<-gsub("X", "", colnames(rawcounts))
colnames(rawcounts)
```

```
##  [1] "miRNA"      "read_count" "precursor"  "total"      "001"       
##  [6] "002"        "003"        "004"        "005"        "006"       
## [11] "007"        "008"        "009"        "010"        "011"       
## [16] "012"        "013"        "014"        "015"        "016"       
## [21] "017"        "018"        "019"        "020"        "021"       
## [26] "022"        "023"        "024"        "001.norm."  "002.norm." 
## [31] "003.norm."  "004.norm."  "005.norm."  "006.norm."  "007.norm." 
## [36] "008.norm."  "009.norm."  "010.norm."  "011.norm."  "012.norm." 
## [41] "013.norm."  "014.norm."  "015.norm."  "016.norm."  "017.norm." 
## [46] "018.norm."  "019.norm."  "020.norm."  "021.norm."  "022.norm." 
## [51] "023.norm."  "024.norm."
```

View the data set:


```r
head(rawcounts[1:8])
```

```
##           miRNA read_count    precursor   total   001   002  003   004
## 1    ssc-let-7a    1416064 ssc-let-7a-1 1416064 24170 53629 6007 57842
## 2    ssc-let-7a    1416269 ssc-let-7a-2 1416269 24179 53662 6022 57862
## 3    ssc-let-7c     174980   ssc-let-7c  174980  3001  5920 1765  8426
## 4 ssc-let-7d-5p      80709   ssc-let-7d   80709  1246  3067  364  3617
## 5 ssc-let-7d-3p      14107   ssc-let-7d   14107   211   232 2040   643
## 6    ssc-let-7e      56528   ssc-let-7e   56528  1026  2068  270  2582
```

```r
colclass<-NULL

for (i in colnames(rawcounts)) {
   colclass<-c(colclass,class(rawcounts[,i]))
}

table(colclass)
```

```
## colclass
##  factor numeric 
##       2      50
```

```r
head(colclass)
```

```
## [1] "factor"  "numeric" "factor"  "numeric" "numeric" "numeric"
```

Notice columns 1 and 3 are factors - the mature miRNA names and the miRNA precursor names.
### 2. Read in the config file, maintaining the characters in the 3-digit code names


```r
configfile<-read.table("../../1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/6_prepped_fasta_24illumina_output/1_config_24illumina_mapper.txt", header = FALSE, sep = " ", row.names=NULL, colClasses = c('character','character'))
head(configfile)
```

```
##                     V1  V2
## 1 1034_mapper_input.fa 001
## 2 1058_mapper_input.fa 002
## 3 1080_mapper_input.fa 003
## 4 1096_mapper_input.fa 004
## 5 1116_mapper_input.fa 005
## 6 1134_mapper_input.fa 006
```

Remove the "_mapper_input.fa" from each file name to leave the pig id:


```r
configfile$V1<-gsub("_mapper_input.fa", "", configfile$V1)
head(configfile)
```

```
##     V1  V2
## 1 1034 001
## 2 1058 002
## 3 1080 003
## 4 1096 004
## 5 1116 005
## 6 1134 006
```

```r
colnames(configfile)<-c("pigid", "code")
head(configfile)
```

```
##   pigid code
## 1  1034  001
## 2  1058  002
## 3  1080  003
## 4  1096  004
## 5  1116  005
## 6  1134  006
```

## Analysis
### 1. Extract the columns of mature read counts for each miRNA for each animal


```r
mirquant<-rawcounts[,c(1,5:28)]
colnames(mirquant)
```

```
##  [1] "miRNA" "001"   "002"   "003"   "004"   "005"   "006"   "007"  
##  [9] "008"   "009"   "010"   "011"   "012"   "013"   "014"   "015"  
## [17] "016"   "017"   "018"   "019"   "020"   "021"   "022"   "023"  
## [25] "024"
```

```r
dim(mirquant)
```

```
## [1] 492  25
```

```r
head(mirquant[1:8])
```

```
##           miRNA   001   002  003   004   005   006   007
## 1    ssc-let-7a 24170 53629 6007 57842 56616 61336 59220
## 2    ssc-let-7a 24179 53662 6022 57862 56649 61339 59249
## 3    ssc-let-7c  3001  5920 1765  8426  7084  6985  8230
## 4 ssc-let-7d-5p  1246  3067  364  3617  3382  3228  3280
## 5 ssc-let-7d-3p   211   232 2040   643   661   493   702
## 6    ssc-let-7e  1026  2068  270  2582  2601  2377  2446
```

Take a subset of this data.frame for testing:


```r
test <- mirquant[1:20,1:8]

test[1:10,]
```

```
##            miRNA   001   002  003   004   005   006   007
## 1     ssc-let-7a 24170 53629 6007 57842 56616 61336 59220
## 2     ssc-let-7a 24179 53662 6022 57862 56649 61339 59249
## 3     ssc-let-7c  3001  5920 1765  8426  7084  6985  8230
## 4  ssc-let-7d-5p  1246  3067  364  3617  3382  3228  3280
## 5  ssc-let-7d-3p   211   232 2040   643   661   493   702
## 6     ssc-let-7e  1026  2068  270  2582  2601  2377  2446
## 7     ssc-let-7f 14082 38536 2150 19808 25335 48079 18801
## 8     ssc-let-7f 13703 37460 2108 19432 24862 46999 18466
## 9     ssc-let-7g  5005 12008  945  9915 14049 12378 10013
## 10    ssc-let-7i  4073  9769  567  7695  6382  9685  6804
```

### 2. Use the 'by' function to apply the function colMeans to the entire data frame of read counts:
(What this will do is go down the columns looking at the index of grouped miRNA names and take the average of the read counts for that miRNA)
The result of this will be a list containing the average read counts for each miRNA for each animal.

Example: by(data, index, function)


```r
bytst<-by(test[,2:ncol(test)], test[,1], colMeans)
bytst[1:25]
```

```
## $`ssc-let-7a`
##     001     002     003     004     005     006     007 
## 24174.5 53645.5  6014.5 57852.0 56632.5 61337.5 59234.5 
## 
## $`ssc-let-7c`
##  001  002  003  004  005  006  007 
## 3001 5920 1765 8426 7084 6985 8230 
## 
## $`ssc-let-7d-3p`
##  001  002  003  004  005  006  007 
##  211  232 2040  643  661  493  702 
## 
## $`ssc-let-7d-5p`
##  001  002  003  004  005  006  007 
## 1246 3067  364 3617 3382 3228 3280 
## 
## $`ssc-let-7e`
##  001  002  003  004  005  006  007 
## 1026 2068  270 2582 2601 2377 2446 
## 
## $`ssc-let-7f`
##     001     002     003     004     005     006     007 
## 13892.5 37998.0  2129.0 19620.0 25098.5 47539.0 18633.5 
## 
## $`ssc-let-7g`
##   001   002   003   004   005   006   007 
##  5005 12008   945  9915 14049 12378 10013 
## 
## $`ssc-let-7i`
##  001  002  003  004  005  006  007 
## 4073 9769  567 7695 6382 9685 6804 
## 
## $`ssc-miR-1`
##   001   002   003   004   005   006   007 
##  2372 10823   821  2637  3576  8515  2385 
## 
## $`ssc-miR-100`
##   001   002   003   004   005   006   007 
## 21034 35674  3049 38984 50318 38032 31860 
## 
## $`ssc-miR-101`
##    001    002    003    004    005    006    007 
## 1253.5 6102.0  322.5 1671.5 2064.0 5642.0 1887.0 
## 
## $`ssc-miR-103`
##    001    002    003    004    005    006    007 
##  865.5 1286.0  218.5 1881.0 1242.5 1057.5 1522.0 
## 
## $`ssc-miR-105-1`
## 001 002 003 004 005 006 007 
##   0   0   0   0   0   0   2 
## 
## $`ssc-miR-105-2`
## 001 002 003 004 005 006 007 
##   0   0   0   0   1   0   0 
## 
## $`ssc-miR-106a`
## 001 002 003 004 005 006 007 
##   1   7   0   4   3   4   5 
## 
## $`ssc-miR-107`
## 001 002 003 004 005 006 007 
## 165 312  49 448 317 293 403 
## 
## $`ssc-miR-10a-3p`
## NULL
## 
## $`ssc-miR-10a-5p`
## NULL
## 
## $`ssc-miR-10b`
## NULL
## 
## $`ssc-miR-122`
## NULL
## 
## $`ssc-miR-1224`
## NULL
## 
## $`ssc-miR-1249`
## NULL
## 
## $`ssc-miR-124a`
## NULL
## 
## $`ssc-miR-125a`
## NULL
## 
## $`ssc-miR-125b`
## NULL
```

Notice here that the rest of the miRNA names remain since the miRNA name is a factor, but since there is no data for them they are filled with NULL.

Apply the by function to the full dataframe:


```r
meanrc<-by(mirquant[,2:ncol(mirquant)], mirquant[,1], colMeans)
```

This should be 411 (the number of mature pig miRNAs in miRBase), meaning we have one expression profile for each mature miRNA:


```r
length(meanrc)
```

```
## [1] 411
```

```r
head(meanrc)
```

```
## $`ssc-let-7a`
##      001      002      003      004      005      006      007      008 
##  24174.5  53645.5   6014.5  57852.0  56632.5  61337.5  59234.5  57940.0 
##      009      010      011      012      013      014      015      016 
##  35574.0  35545.0  39448.5  36966.5  76445.5 219999.0  96168.5 106307.5 
##      017      018      019      020      021      022      023      024 
##  19861.5  14706.5  34628.0  78745.0 124874.0  42307.0  29833.5  47925.5 
## 
## $`ssc-let-7c`
##   001   002   003   004   005   006   007   008   009   010   011   012 
##  3001  5920  1765  8426  7084  6985  8230  6396  3996  6624  6101  4591 
##   013   014   015   016   017   018   019   020   021   022   023   024 
##  8332 29347 10912 13799  1940  2532  3603  9650 13632  3856  3160  5098 
## 
## $`ssc-let-7d-3p`
##  001  002  003  004  005  006  007  008  009  010  011  012  013  014  015 
##  211  232 2040  643  661  493  702  338  220 1540  518  407  337 1643  405 
##  016  017  018  019  020  021  022  023  024 
## 1110   95  580  156  524  624  213  170  245 
## 
## $`ssc-let-7d-5p`
##   001   002   003   004   005   006   007   008   009   010   011   012 
##  1246  3067   364  3617  3382  3228  3280  3155  1868  2272  2115  2202 
##   013   014   015   016   017   018   019   020   021   022   023   024 
##  4586 11705  5192  5524  1245  1030  2010  5003  7418  2827  1795  2578 
## 
## $`ssc-let-7e`
##  001  002  003  004  005  006  007  008  009  010  011  012  013  014  015 
## 1026 2068  270 2582 2601 2377 2446 2395 1102 1404 1569 1309 3011 8544 3578 
##  016  017  018  019  020  021  022  023  024 
## 4294  783  702 1382 3137 5032 1682 1503 1731 
## 
## $`ssc-let-7f`
##      001      002      003      004      005      006      007      008 
##  13892.5  37998.0   2129.0  19620.0  25098.5  47539.0  18633.5  41081.5 
##      009      010      011      012      013      014      015      016 
##  27063.5  25633.5  26353.0  20855.5  58404.5 203564.5  67330.0  37886.0 
##      017      018      019      020      021      022      023      024 
##  13534.0   7414.0  26159.0  38624.5 103977.5  35902.5  20510.0  36413.5
```


### 3. Transform the list output back into a data.frame using the plyr package to prepare for gblup function:
Example: ldply(.data, .fun, .id)

id = name of the index column (used if data is a named list). Pass NULL to avoid creation
     of the index column. For compatibility, omit this argument or pass "NA" to avoid converting the index column
     to a factor; in this case, ".id" is used as column name.


```r
dfmeanrc<-ldply(meanrc, fun=NULL, id=names(meanrc))

head(dfmeanrc[1:8])
```

```
##             .id     001     002    003   004     005     006     007
## 1    ssc-let-7a 24174.5 53645.5 6014.5 57852 56632.5 61337.5 59234.5
## 2    ssc-let-7c  3001.0  5920.0 1765.0  8426  7084.0  6985.0  8230.0
## 3 ssc-let-7d-3p   211.0   232.0 2040.0   643   661.0   493.0   702.0
## 4 ssc-let-7d-5p  1246.0  3067.0  364.0  3617  3382.0  3228.0  3280.0
## 5    ssc-let-7e  1026.0  2068.0  270.0  2582  2601.0  2377.0  2446.0
## 6    ssc-let-7f 13892.5 37998.0 2129.0 19620 25098.5 47539.0 18633.5
```

```r
dim(dfmeanrc)
```

```
## [1] 411  25
```

These dimensions are what would be expected, because there are 411 mature sus scrofa miRNA sequences in miRBase,
and there are 24 animals in the analysis, plus the miRNA column.

Check that the correct miRNA name went with the correct data:


```r
if (sum(names(meanrc)!=dfmeanrc[,1]) != 0) stop ("miRNA names are not the same")

colnames(dfmeanrc)[[1]]<-"miRNA"

if (sum(colnames(dfmeanrc)!=colnames(mirquant)) != 0) stop ("animal order not the same")

head(dfmeanrc[,1:10])
```

```
##           miRNA     001     002    003   004     005     006     007
## 1    ssc-let-7a 24174.5 53645.5 6014.5 57852 56632.5 61337.5 59234.5
## 2    ssc-let-7c  3001.0  5920.0 1765.0  8426  7084.0  6985.0  8230.0
## 3 ssc-let-7d-3p   211.0   232.0 2040.0   643   661.0   493.0   702.0
## 4 ssc-let-7d-5p  1246.0  3067.0  364.0  3617  3382.0  3228.0  3280.0
## 5    ssc-let-7e  1026.0  2068.0  270.0  2582  2601.0  2377.0  2446.0
## 6    ssc-let-7f 13892.5 37998.0 2129.0 19620 25098.5 47539.0 18633.5
##       008     009
## 1 57940.0 35574.0
## 2  6396.0  3996.0
## 3   338.0   220.0
## 4  3155.0  1868.0
## 5  2395.0  1102.0
## 6 41081.5 27063.5
```

### 4. Filter the data for expression threshold: The total read count for the miRNA needs to be greater than 0

Set first column of dfmeanrc (miRNA ids) as the row.names:


```r
rownames(dfmeanrc)<-dfmeanrc$miRNA
```

Eliminate column of row names:


```r
dfmeanrc<-dfmeanrc[,-c(1)]
head(dfmeanrc[,1:10])
```

```
##                   001     002    003   004     005     006     007     008
## ssc-let-7a    24174.5 53645.5 6014.5 57852 56632.5 61337.5 59234.5 57940.0
## ssc-let-7c     3001.0  5920.0 1765.0  8426  7084.0  6985.0  8230.0  6396.0
## ssc-let-7d-3p   211.0   232.0 2040.0   643   661.0   493.0   702.0   338.0
## ssc-let-7d-5p  1246.0  3067.0  364.0  3617  3382.0  3228.0  3280.0  3155.0
## ssc-let-7e     1026.0  2068.0  270.0  2582  2601.0  2377.0  2446.0  2395.0
## ssc-let-7f    13892.5 37998.0 2129.0 19620 25098.5 47539.0 18633.5 41081.5
##                   009     010
## ssc-let-7a    35574.0 35545.0
## ssc-let-7c     3996.0  6624.0
## ssc-let-7d-3p   220.0  1540.0
## ssc-let-7d-5p  1868.0  2272.0
## ssc-let-7e     1102.0  1404.0
## ssc-let-7f    27063.5 25633.5
```

```r
dim(dfmeanrc)
```

```
## [1] 411  24
```


How many miRNAs are not expressed (total read count across all libraries = 0)?


```r
head(rowSums(dfmeanrc))
```

```
##    ssc-let-7a    ssc-let-7c ssc-let-7d-3p ssc-let-7d-5p    ssc-let-7e 
##     1416166.5      174980.0       14107.0       80709.0       56528.0 
##    ssc-let-7f 
##      955617.5
```

```r
tail(rowSums(dfmeanrc))
```

```
## ssc-miR-9859-3p ssc-miR-9860-5p ssc-miR-9861-5p ssc-miR-9862-3p 
##               0             528               0               0 
##     ssc-miR-99a     ssc-miR-99b 
##          182058          325183
```

```r
table(rowSums(dfmeanrc)==0)
```

```
## 
## FALSE  TRUE 
##   325    86
```

So, 86 miRNA profiles contain 0 read counts total, meaning 0 expression

Filter the matrix to keep only those miRNAs whose total expression is greater than 0.


```r
illumina24.no.zero.dfmeanrc<-dfmeanrc[rowSums(dfmeanrc)>0,]
dim(illumina24.no.zero.dfmeanrc)
```

```
## [1] 325  24
```

```r
head(illumina24.no.zero.dfmeanrc[,1:10])
```

```
##                   001     002    003   004     005     006     007     008
## ssc-let-7a    24174.5 53645.5 6014.5 57852 56632.5 61337.5 59234.5 57940.0
## ssc-let-7c     3001.0  5920.0 1765.0  8426  7084.0  6985.0  8230.0  6396.0
## ssc-let-7d-3p   211.0   232.0 2040.0   643   661.0   493.0   702.0   338.0
## ssc-let-7d-5p  1246.0  3067.0  364.0  3617  3382.0  3228.0  3280.0  3155.0
## ssc-let-7e     1026.0  2068.0  270.0  2582  2601.0  2377.0  2446.0  2395.0
## ssc-let-7f    13892.5 37998.0 2129.0 19620 25098.5 47539.0 18633.5 41081.5
##                   009     010
## ssc-let-7a    35574.0 35545.0
## ssc-let-7c     3996.0  6624.0
## ssc-let-7d-3p   220.0  1540.0
## ssc-let-7d-5p  1868.0  2272.0
## ssc-let-7e     1102.0  1404.0
## ssc-let-7f    27063.5 25633.5
```

```r
if (sum(rowSums(illumina24.no.zero.dfmeanrc)==0)!= 0) stop ("expression filtering did not work correctly")

if (sum(colnames(illumina24.no.zero.dfmeanrc)!=colnames(dfmeanrc)) != 0) stop ("animal order not the same")
```

### 5. Restore the pig IDs in place of the 3-digit codes as the column names of the data frame, for use with the gblup function of gwaR


```r
head(configfile)
```

```
##   pigid code
## 1  1034  001
## 2  1058  002
## 3  1080  003
## 4  1096  004
## 5  1116  005
## 6  1134  006
```

Now I need to substitute the 3-digit code with the pig IDs, ensuring the names stay in the correct order:

Use match function to find positional index and match column names:

The object dfmeanrc has column names that need to be re-named. I have the config file which contains
the current column names and the desired column names. What I am doing in this code is re-ordering the config file 
based on where the config file "code" column matches the position of the dfmeanrc object's column names, then having it return the corresponding value in column "pigid". 

So, when using match, need to have the first argument be the matrix/dataframe you want to change or match, and the second argument be what you want to index it by or match it against. 

"Where does [vector] match in [matrix]?" or "Match the column names of illumina24.no.zero.dfmeanrc to the configfile "code" column, then return the corresponding pigid."


```r
configfile[match(colnames(illumina24.no.zero.dfmeanrc),configfile$code),"pigid"]
```

```
##  [1] "1034" "1058" "1080" "1096" "1116" "1134" "1154" "1170" "1194" "1240"
## [11] "1278" "1300" "1426" "1434" "1458" "1484" "1502" "1512" "1534" "1580"
## [21] "1594" "1640" "1644" "1662"
```

Assign the column names using match:


```r
colnames(illumina24.no.zero.dfmeanrc)<- configfile[match(colnames(illumina24.no.zero.dfmeanrc),configfile$code),"pigid"]
head(illumina24.no.zero.dfmeanrc[1:10])
```

```
##                  1034    1058   1080  1096    1116    1134    1154    1170
## ssc-let-7a    24174.5 53645.5 6014.5 57852 56632.5 61337.5 59234.5 57940.0
## ssc-let-7c     3001.0  5920.0 1765.0  8426  7084.0  6985.0  8230.0  6396.0
## ssc-let-7d-3p   211.0   232.0 2040.0   643   661.0   493.0   702.0   338.0
## ssc-let-7d-5p  1246.0  3067.0  364.0  3617  3382.0  3228.0  3280.0  3155.0
## ssc-let-7e     1026.0  2068.0  270.0  2582  2601.0  2377.0  2446.0  2395.0
## ssc-let-7f    13892.5 37998.0 2129.0 19620 25098.5 47539.0 18633.5 41081.5
##                  1194    1240
## ssc-let-7a    35574.0 35545.0
## ssc-let-7c     3996.0  6624.0
## ssc-let-7d-3p   220.0  1540.0
## ssc-let-7d-5p  1868.0  2272.0
## ssc-let-7e     1102.0  1404.0
## ssc-let-7f    27063.5 25633.5
```

```r
dim(illumina24.no.zero.dfmeanrc)
```

```
## [1] 325  24
```

```r
if (sum(colnames(illumina24.no.zero.dfmeanrc)!=(configfile$pigid))!=0) stop ("match function did not work correctly")
```

### 6. Round the mean read counts to the nearest integer for use with the voom function


```r
illumina24.no.zero.dfmeanrc[1:10,1:6]
```

```
##                  1034    1058   1080  1096    1116    1134
## ssc-let-7a    24174.5 53645.5 6014.5 57852 56632.5 61337.5
## ssc-let-7c     3001.0  5920.0 1765.0  8426  7084.0  6985.0
## ssc-let-7d-3p   211.0   232.0 2040.0   643   661.0   493.0
## ssc-let-7d-5p  1246.0  3067.0  364.0  3617  3382.0  3228.0
## ssc-let-7e     1026.0  2068.0  270.0  2582  2601.0  2377.0
## ssc-let-7f    13892.5 37998.0 2129.0 19620 25098.5 47539.0
## ssc-let-7g     5005.0 12008.0  945.0  9915 14049.0 12378.0
## ssc-let-7i     4073.0  9769.0  567.0  7695  6382.0  9685.0
## ssc-miR-1      2372.0 10823.0  821.0  2637  3576.0  8515.0
## ssc-miR-100   21034.0 35674.0 3049.0 38984 50318.0 38032.0
```

```r
illumina24.no.zero.dfmeanrcround<-round(illumina24.no.zero.dfmeanrc)

illumina24.no.zero.dfmeanrcround[1:10,1:6]
```

```
##                1034  1058 1080  1096  1116  1134
## ssc-let-7a    24174 53646 6014 57852 56632 61338
## ssc-let-7c     3001  5920 1765  8426  7084  6985
## ssc-let-7d-3p   211   232 2040   643   661   493
## ssc-let-7d-5p  1246  3067  364  3617  3382  3228
## ssc-let-7e     1026  2068  270  2582  2601  2377
## ssc-let-7f    13892 37998 2129 19620 25098 47539
## ssc-let-7g     5005 12008  945  9915 14049 12378
## ssc-let-7i     4073  9769  567  7695  6382  9685
## ssc-miR-1      2372 10823  821  2637  3576  8515
## ssc-miR-100   21034 35674 3049 38984 50318 38032
```

The final matrix of filtered, rounded, mean read counts needs to have miRNA rownames and Animal ID colnames


```r
head(rownames(illumina24.no.zero.dfmeanrcround))
```

```
## [1] "ssc-let-7a"    "ssc-let-7c"    "ssc-let-7d-3p" "ssc-let-7d-5p"
## [5] "ssc-let-7e"    "ssc-let-7f"
```

```r
head(colnames(illumina24.no.zero.dfmeanrcround))
```

```
## [1] "1034" "1058" "1080" "1096" "1116" "1134"
```

```r
if (sum(colnames(illumina24.no.zero.dfmeanrcround)!= configfile$pigid)!= 0) stop ("rownames do not match pigid")
```

## Save data
What I am saving here is the filtered, rounded, average read counts in an .Rdata object


```r
save(illumina24.no.zero.dfmeanrcround, file = "../1_24illumina_filtered_rounded_mean_mature_mirna_exp.Rdata")
```

