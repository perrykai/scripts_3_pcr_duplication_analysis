**Script:** `6_extract_24bioosci_PCRDUP_mirna_expression_df.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/scripts`

**Date:**  9/22/16

**Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/3_mirdeep2_core_output_24bioosci_PCRDUP_samples`

**Input File(s):** `miRNAs_expressed_all_samples_22_09_2016_t_14_48_56.csv`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/4_24bioosci_PCRDUP_mirna_expression_matrix/`

**Output File(s):** `1_24bioosci_PCRDUP_rounded_mean_mature_mirna_expression.Rdata`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Save data](#save-data)

## Objectives
The objective of this script is to create a matrix of miRNA expression data (read counts) for incorporation into the miRNA DE analysis, beginning with the known miRNA expression data output from the miRDeep2 core module.
Later, expression data of the putative novel miRNA candidates can also be included, also output from the miRDeep2 core module. 
To acheive one read count per miRNA per animal, I will need to take the average of the mature read counts in instances of multiple precursor sequences.

The result of this script will be a data frame containing one average read count per miRNA per animal.

So, what I need to do:

1. Extract the columns of the mature miRNA read counts for each animal
2. Use the 'by' function to apply the function colmeans to the entire data frame of read counts (What this will do is go down the columns looking at the index of grouped miRNA names and take the average of the read counts for that group of miRNA. The result of this will be a list containing the average read counts for each miRNA for each animal)
3. Transform the list output back into a data.frame using the plyr package to prepare for gblup function
4. Restore the pig IDs in place of the 3-digit codes as the column names of the data frame, for use with the gblup function of gwaR
5. Round the mean read counts to the nearest integer for use with the voom function
6. Save the final data.frame.
## Install libraries


```r
rm(list=ls())
library(plyr)
```

## Load data


```r
setwd("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/scripts")
```

### 1. Read in the csv of the miRNA expression profiles, using check.names default to correct column names with non-ASCII characters


```r
rc<-read.csv("../3_mirdeep2_core_output_24bioosci_PCRDUP_samples/miRNAs_expressed_all_samples_22_09_2016_t_14_48_56.csv", sep = "\t", header = TRUE, row.names=NULL)
```

Set the name of the first column to "miRNA":


```r
colnames(rc)[[1]]<-"miRNA"
```

Remove the "X" character from the beginning of each column:


```r
colnames(rc)<-gsub("X", "", colnames(rc))
colnames(rc)
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
head(rc[1:8])
```

```
##           miRNA read_count    precursor   total    001    002   003   004
## 1    ssc-let-7a    1910808 ssc-let-7a-1 1910808 121924 120597 99533 70953
## 2    ssc-let-7a    1908440 ssc-let-7a-2 1908440 121818 120381 99530 70862
## 3    ssc-let-7c     952357   ssc-let-7c  952357  52830  77943 51631 28774
## 4 ssc-let-7d-5p     117102   ssc-let-7d  117102   7371   6641  6233  4304
## 5 ssc-let-7d-3p      18254   ssc-let-7d   18254    489   1396  1213   535
## 6    ssc-let-7e      62426   ssc-let-7e   62426   3636   4248  3334  2285
```

Do a for loop to check the class of each column in the data.frame


```r
colclass<-NULL

for (i in colnames(rc)) {
   colclass<-c(colclass,class(rc[,i]))
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

Notice here that the mature miRNA names and the miRNA precursor names are considered factors (columns 1 and 3).
### 2. Read in the config file, maintaining the characters in the 3-digit code names


```r
configfile<-read.table("../1_mirdeep2_mapper_input_24bioosci_PCRDUP_samples/1_config_24bioosci_PCRDUP_mapper.txt", header = FALSE, sep = " ", row.names=NULL, colClasses = c('character','character'))
head(configfile)
```

```
##                                      V1  V2
## 1 1034_24bioosci_PCRDUP_mapper_input.fa 001
## 2 1058_24bioosci_PCRDUP_mapper_input.fa 002
## 3 1080_24bioosci_PCRDUP_mapper_input.fa 003
## 4 1096_24bioosci_PCRDUP_mapper_input.fa 004
## 5 1116_24bioosci_PCRDUP_mapper_input.fa 005
## 6 1134_24bioosci_PCRDUP_mapper_input.fa 006
```

Remove the "_express.fa" from each file name to leave the pig id:


```r
configfile$V1<-gsub("_24bioosci_PCRDUP_mapper_input.fa", "", configfile$V1)
```

Make filenames more informative:


```r
colnames(configfile)<-c("pigid","code")
colnames(configfile)
```

```
## [1] "pigid" "code"
```

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

## Analysis
### 1. Extract the columns of mature read counts for each miRNA for each animal


```r
mirquant<-rc[,c(1,5:28)]
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
##           miRNA    001    002   003   004   005    006    007
## 1    ssc-let-7a 121924 120597 99533 70953 88472 185087 105395
## 2    ssc-let-7a 121818 120381 99530 70862 88396 184823 105280
## 3    ssc-let-7c  52830  77943 51631 28774 52407  90522  69993
## 4 ssc-let-7d-5p   7371   6641  6233  4304  4989  11692   8054
## 5 ssc-let-7d-3p    489   1396  1213   535  1042   1594   1123
## 6    ssc-let-7e   3636   4248  3334  2285  3023   5894   3186
```

Take a subset of this data.frame for testing:


```r
test<-mirquant[1:20,1:8]

test[1:10,]
```

```
##            miRNA    001    002   003   004   005    006    007
## 1     ssc-let-7a 121924 120597 99533 70953 88472 185087 105395
## 2     ssc-let-7a 121818 120381 99530 70862 88396 184823 105280
## 3     ssc-let-7c  52830  77943 51631 28774 52407  90522  69993
## 4  ssc-let-7d-5p   7371   6641  6233  4304  4989  11692   8054
## 5  ssc-let-7d-3p    489   1396  1213   535  1042   1594   1123
## 6     ssc-let-7e   3636   4248  3334  2285  3023   5894   3186
## 7     ssc-let-7f  80027  35826 38816 56235 24187  82762  53927
## 8     ssc-let-7f  79541  36369 38915 56012 24477  82992  53867
## 9     ssc-let-7g  33758  25786 26002 24050 16175  45276  36200
## 10    ssc-let-7i  16206  10927  9434  8780  6770  16013  17468
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
##      001      002      003      004      005      006      007 
## 121871.0 120489.0  99531.5  70907.5  88434.0 184955.0 105337.5 
## 
## $`ssc-let-7c`
##   001   002   003   004   005   006   007 
## 52830 77943 51631 28774 52407 90522 69993 
## 
## $`ssc-let-7d-3p`
##  001  002  003  004  005  006  007 
##  489 1396 1213  535 1042 1594 1123 
## 
## $`ssc-let-7d-5p`
##   001   002   003   004   005   006   007 
##  7371  6641  6233  4304  4989 11692  8054 
## 
## $`ssc-let-7e`
##  001  002  003  004  005  006  007 
## 3636 4248 3334 2285 3023 5894 3186 
## 
## $`ssc-let-7f`
##     001     002     003     004     005     006     007 
## 79784.0 36097.5 38865.5 56123.5 24332.0 82877.0 53897.0 
## 
## $`ssc-let-7g`
##   001   002   003   004   005   006   007 
## 33758 25786 26002 24050 16175 45276 36200 
## 
## $`ssc-let-7i`
##   001   002   003   004   005   006   007 
## 16206 10927  9434  8780  6770 16013 17468 
## 
## $`ssc-miR-1`
##     001     002     003     004     005     006     007 
## 3000633  560629  635504 3194495  350213 1521754 1351120 
## 
## $`ssc-miR-100`
##   001   002   003   004   005   006   007 
## 14875 40497 24024  5421 19104 22891 25681 
## 
## $`ssc-miR-101`
##     001     002     003     004     005     006     007 
## 17862.5  7990.5  8471.5 12772.0  4261.0 11461.0 10299.0 
## 
## $`ssc-miR-103`
##    001    002    003    004    005    006    007 
## 2991.5 3343.0 2910.5 1727.0 2337.0 4359.5 4040.5 
## 
## $`ssc-miR-105-1`
## 001 002 003 004 005 006 007 
##   5   2   0   0   0   1   4 
## 
## $`ssc-miR-105-2`
## 001 002 003 004 005 006 007 
##   1   0   0   3   0   0   3 
## 
## $`ssc-miR-106a`
## 001 002 003 004 005 006 007 
##  44  41  27  28  19  44  63 
## 
## $`ssc-miR-107`
## 001 002 003 004 005 006 007 
## 605 618 444 277 392 755 634 
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
## 121871.0 120489.0  99531.5  70907.5  88434.0 184955.0 105337.5   4604.5 
##      009      010      011      012      013      014      015      016 
##  56074.0  74954.0  31466.5  62566.5  61593.0  50682.5  42404.5  66980.0 
##      017      018      019      020      021      022      023      024 
## 128323.0  44153.5  47459.0  75309.0  79047.0  99397.5  96866.5  96217.5 
## 
## $`ssc-let-7c`
##   001   002   003   004   005   006   007   008   009   010   011   012 
## 52830 77943 51631 28774 52407 90522 69993  6221 23344 35784 24846 36361 
##   013   014   015   016   017   018   019   020   021   022   023   024 
## 26465 24293 22798 39871 48683 32208 20438 39214 33516 32698 39353 42164 
## 
## $`ssc-let-7d-3p`
##  001  002  003  004  005  006  007  008  009  010  011  012  013  014  015 
##  489 1396 1213  535 1042 1594 1123 1166  290  340 1786  936  189  244  161 
##  016  017  018  019  020  021  022  023  024 
##  742  549 1740  195  671  412  393  413  635 
## 
## $`ssc-let-7d-5p`
##   001   002   003   004   005   006   007   008   009   010   011   012 
##  7371  6641  6233  4304  4989 11692  8054   334  3189  4127  2278  3911 
##   013   014   015   016   017   018   019   020   021   022   023   024 
##  3813  2523  2967  4139  7365  3184  2772  4662  5003  6280  5437  5834 
## 
## $`ssc-let-7e`
##  001  002  003  004  005  006  007  008  009  010  011  012  013  014  015 
## 3636 4248 3334 2285 3023 5894 3186  146 1995 2174  871 1929 2255 1625 1326 
##  016  017  018  019  020  021  022  023  024 
## 2430 4168 1566 1708 2362 2635 2883 3592 3155 
## 
## $`ssc-let-7f`
##     001     002     003     004     005     006     007     008     009 
## 79784.0 36097.5 38865.5 56123.5 24332.0 82877.0 53897.0  1089.5 33521.5 
##     010     011     012     013     014     015     016     017     018 
## 52903.5 18376.0 35595.5 41933.0 32378.0 28576.5 20309.0 85063.0 18608.0 
##     019     020     021     022     023     024 
## 32601.0 33097.5 54375.0 79049.0 63528.0 66253.5
```


### 3. Transform the list output back into a data.frame using the plyr package to prepare for gblup function:
Example: ldply(.data, .fun, .id)

id = name of the index column (used if data is a named list). Pass NULL to avoid creation
     of the index column. For compatibility, omit this argument or pass "NA" to avoid converting the index column
     to a factor; in this case, ".id" is used as column name.


```r
pcrdup.bioo24.dfmeanrc<-ldply(meanrc, fun=NULL, id=names(meanrc))

head(pcrdup.bioo24.dfmeanrc[1:8])
```

```
##             .id    001      002     003     004   005    006      007
## 1    ssc-let-7a 121871 120489.0 99531.5 70907.5 88434 184955 105337.5
## 2    ssc-let-7c  52830  77943.0 51631.0 28774.0 52407  90522  69993.0
## 3 ssc-let-7d-3p    489   1396.0  1213.0   535.0  1042   1594   1123.0
## 4 ssc-let-7d-5p   7371   6641.0  6233.0  4304.0  4989  11692   8054.0
## 5    ssc-let-7e   3636   4248.0  3334.0  2285.0  3023   5894   3186.0
## 6    ssc-let-7f  79784  36097.5 38865.5 56123.5 24332  82877  53897.0
```

```r
dim(pcrdup.bioo24.dfmeanrc)
```

```
## [1] 411  25
```

These dimensions are what would be expected, because there are 411 mature sus scrofa miRNA sequences in miRBase,
and there are 174 animals in the analysis, plus the miRNA column.

Check that the correct miRNA name went with the correct data:


```r
if (sum(names(meanrc)!=pcrdup.bioo24.dfmeanrc[,1]) != 0) stop ("miRNA names are not the same")

colnames(pcrdup.bioo24.dfmeanrc)[[1]]<-"miRNA"

if (sum(colnames(pcrdup.bioo24.dfmeanrc)!=colnames(mirquant)) != 0) stop ("animal order not the same")

head(pcrdup.bioo24.dfmeanrc[,1:10])
```

```
##           miRNA    001      002     003     004   005    006      007
## 1    ssc-let-7a 121871 120489.0 99531.5 70907.5 88434 184955 105337.5
## 2    ssc-let-7c  52830  77943.0 51631.0 28774.0 52407  90522  69993.0
## 3 ssc-let-7d-3p    489   1396.0  1213.0   535.0  1042   1594   1123.0
## 4 ssc-let-7d-5p   7371   6641.0  6233.0  4304.0  4989  11692   8054.0
## 5    ssc-let-7e   3636   4248.0  3334.0  2285.0  3023   5894   3186.0
## 6    ssc-let-7f  79784  36097.5 38865.5 56123.5 24332  82877  53897.0
##      008     009
## 1 4604.5 56074.0
## 2 6221.0 23344.0
## 3 1166.0   290.0
## 4  334.0  3189.0
## 5  146.0  1995.0
## 6 1089.5 33521.5
```


Set first column of pcrdup.bioo24.dfmeanrc (miRNA ids) as the row.names:


```r
rownames(pcrdup.bioo24.dfmeanrc)<-pcrdup.bioo24.dfmeanrc$miRNA
```

Eliminate column of row names:


```r
pcrdup.bioo24.dfmeanrc<-pcrdup.bioo24.dfmeanrc[,-c(1)]
head(pcrdup.bioo24.dfmeanrc[,1:10])
```

```
##                  001      002     003     004   005    006      007    008
## ssc-let-7a    121871 120489.0 99531.5 70907.5 88434 184955 105337.5 4604.5
## ssc-let-7c     52830  77943.0 51631.0 28774.0 52407  90522  69993.0 6221.0
## ssc-let-7d-3p    489   1396.0  1213.0   535.0  1042   1594   1123.0 1166.0
## ssc-let-7d-5p   7371   6641.0  6233.0  4304.0  4989  11692   8054.0  334.0
## ssc-let-7e      3636   4248.0  3334.0  2285.0  3023   5894   3186.0  146.0
## ssc-let-7f     79784  36097.5 38865.5 56123.5 24332  82877  53897.0 1089.5
##                   009     010
## ssc-let-7a    56074.0 74954.0
## ssc-let-7c    23344.0 35784.0
## ssc-let-7d-3p   290.0   340.0
## ssc-let-7d-5p  3189.0  4127.0
## ssc-let-7e     1995.0  2174.0
## ssc-let-7f    33521.5 52903.5
```

```r
dim(pcrdup.bioo24.dfmeanrc)
```

```
## [1] 411  24
```

### 4. Restore the pig IDs in place of the 3-digit codes as the column names of the data frame, for use with the gblup function of gwaR


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

The object pcrdup.bioo24.dfmeanrc has column names that need to be re-named. I have the config file which contains
the current column names and the desired column names. What I am doing in this code is re-ordering the config file 
based on where the config file "code" column matches the position of the pcrdup.bioo24.dfmeanrc object's column names, then having it return the corresponding value in column "pigid". 

So, when using match, need to have the first argument be the matrix/dataframe you want to change or match, and the second argument be what you want to index it by or match it against. 

"Where does [vector] match in [matrix]?" or "Match the column names of pcrdup.bioo24.dfmeanrc to the configfile "code" column, then return the corresponding pigid."


```r
configfile[match(colnames(pcrdup.bioo24.dfmeanrc),configfile$code),"pigid"]
```

```
##  [1] "1034" "1058" "1080" "1096" "1116" "1134" "1154" "1170" "1194" "1240"
## [11] "1278" "1300" "1426" "1434" "1458" "1484" "1502" "1512" "1534" "1580"
## [21] "1594" "1640" "1644" "1662"
```

Assign the column names using match:


```r
colnames(pcrdup.bioo24.dfmeanrc)<- configfile[match(colnames(pcrdup.bioo24.dfmeanrc),configfile$code),"pigid"]
head(pcrdup.bioo24.dfmeanrc[1:10])
```

```
##                 1034     1058    1080    1096  1116   1134     1154   1170
## ssc-let-7a    121871 120489.0 99531.5 70907.5 88434 184955 105337.5 4604.5
## ssc-let-7c     52830  77943.0 51631.0 28774.0 52407  90522  69993.0 6221.0
## ssc-let-7d-3p    489   1396.0  1213.0   535.0  1042   1594   1123.0 1166.0
## ssc-let-7d-5p   7371   6641.0  6233.0  4304.0  4989  11692   8054.0  334.0
## ssc-let-7e      3636   4248.0  3334.0  2285.0  3023   5894   3186.0  146.0
## ssc-let-7f     79784  36097.5 38865.5 56123.5 24332  82877  53897.0 1089.5
##                  1194    1240
## ssc-let-7a    56074.0 74954.0
## ssc-let-7c    23344.0 35784.0
## ssc-let-7d-3p   290.0   340.0
## ssc-let-7d-5p  3189.0  4127.0
## ssc-let-7e     1995.0  2174.0
## ssc-let-7f    33521.5 52903.5
```

```r
dim(pcrdup.bioo24.dfmeanrc)
```

```
## [1] 411  24
```

```r
if (sum(colnames(pcrdup.bioo24.dfmeanrc)!=(configfile$pigid))!=0) stop ("match function did not work correctly")
```

### 5. Round the mean read counts to the nearest integer for use with the voom function


```r
pcrdup.bioo24.dfmeanrc[1:10,1:6]
```

```
##                  1034     1058     1080      1096   1116    1134
## ssc-let-7a     121871 120489.0  99531.5   70907.5  88434  184955
## ssc-let-7c      52830  77943.0  51631.0   28774.0  52407   90522
## ssc-let-7d-3p     489   1396.0   1213.0     535.0   1042    1594
## ssc-let-7d-5p    7371   6641.0   6233.0    4304.0   4989   11692
## ssc-let-7e       3636   4248.0   3334.0    2285.0   3023    5894
## ssc-let-7f      79784  36097.5  38865.5   56123.5  24332   82877
## ssc-let-7g      33758  25786.0  26002.0   24050.0  16175   45276
## ssc-let-7i      16206  10927.0   9434.0    8780.0   6770   16013
## ssc-miR-1     3000633 560629.0 635504.0 3194495.0 350213 1521754
## ssc-miR-100     14875  40497.0  24024.0    5421.0  19104   22891
```

```r
pcrdup.bioo24.dfmeanrcround<-round(pcrdup.bioo24.dfmeanrc)

pcrdup.bioo24.dfmeanrcround[1:10,1:6]
```

```
##                  1034   1058   1080    1096   1116    1134
## ssc-let-7a     121871 120489  99532   70908  88434  184955
## ssc-let-7c      52830  77943  51631   28774  52407   90522
## ssc-let-7d-3p     489   1396   1213     535   1042    1594
## ssc-let-7d-5p    7371   6641   6233    4304   4989   11692
## ssc-let-7e       3636   4248   3334    2285   3023    5894
## ssc-let-7f      79784  36098  38866   56124  24332   82877
## ssc-let-7g      33758  25786  26002   24050  16175   45276
## ssc-let-7i      16206  10927   9434    8780   6770   16013
## ssc-miR-1     3000633 560629 635504 3194495 350213 1521754
## ssc-miR-100     14875  40497  24024    5421  19104   22891
```

The final matrix of rounded, mean read counts needs to have miRNA rownames and Animal ID colnames


```r
head(rownames(pcrdup.bioo24.dfmeanrcround))
```

```
## [1] "ssc-let-7a"    "ssc-let-7c"    "ssc-let-7d-3p" "ssc-let-7d-5p"
## [5] "ssc-let-7e"    "ssc-let-7f"
```

```r
head(colnames(pcrdup.bioo24.dfmeanrcround))
```

```
## [1] "1034" "1058" "1080" "1096" "1116" "1134"
```

```r
if (sum(colnames(pcrdup.bioo24.dfmeanrcround)!= configfile$pigid)!= 0) stop ("rownames do not match pigid")
```

## Save data
What I am saving here is the rounded, average read counts in an .Rdata object


```r
save(pcrdup.bioo24.dfmeanrcround, file = "../4_24bioosci_PCRDUP_mirna_expression_matrix/1_24bioosci_PCRDUP_rounded_mean_mature_mirna_expression.Rdata")
```

