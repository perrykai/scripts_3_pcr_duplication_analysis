**Script:** `1_extract_24bioosci_libraries_for_quantification.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/3_24bioosci_mirna_expression_characterization/scripts/`

**Date:**  8/17/16

**Input File Directory:**  

1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/2_create_expression_matrix_of_known_mirna/`
2. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/`
3. `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/6_prepped_fasta_24illumina_output/`

**Input File(s):** 

1. `1_exp_filtered_rounded_mean_mature_mirna_expression.Rdata`
2. `2_mature_mirna_annotation.Rdata`
3. `1_config_24illumina_mapper.txt`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/3_24bioosci_mirna_expression_characterization/`

**Output File(s):** `1_24bioosci_mirna_expression.Rdata`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Save data](#save-data)

## Objectives

The objective of this script is to extract the 24 libraries from the dataset of 174 Bioo Scientific-prepped libraries matching the 24 Illumina-prepped libraries for the PCR duplication analysis.

The libraries will be extracted from the expression matrix that was used to create the dge object for the 174 libraries, output from the miRNA quantification using miRDeep2. 
The read counts will then be normalized considering only the 24 libraries as a dataset. 

THIS ANALYSIS COMPLETED WITH R/3.2.0
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

setwd("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/3_24bioosci_mirna_expression_characterization/scripts/")
```

## Load data

Load the expression matrix from the 174 Bioo Scientific-prepped data:


```r
load("../../../2_mirna_characterization_expression/2_create_expression_matrix_of_known_mirna/1_exp_filtered_rounded_mean_mature_mirna_expression.Rdata")
```

Load the annotation file for the miRNAs in the dataset:


```r
load("../../../2_mirna_characterization_expression/3_build_dge_object_for_eqtl/2_mature_mirna_annotation.Rdata")
```

Load the config file from the 24 Illumina-prepped libraries:


```r
illuminaconfig<-read.table("../../1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/6_prepped_fasta_24illumina_output/1_config_24illumina_mapper.txt", colClasses=c("character","character"), col.names=c("inputfile", "id"))
head(illuminaconfig)
```

```
##              inputfile  id
## 1 1034_mapper_input.fa 001
## 2 1058_mapper_input.fa 002
## 3 1080_mapper_input.fa 003
## 4 1096_mapper_input.fa 004
## 5 1116_mapper_input.fa 005
## 6 1134_mapper_input.fa 006
```

## Analysis

Check the format of the expression matrix from the 174 Bioo Scientific libraries


```r
no.zero.dfmeanrcround[1:5,1:5]
```

```
##                1034  1036  1041  1049  1058
## ssc-let-7a    48132 23427 28448 29860 40758
## ssc-let-7c    32745 14987 18144 18681 34313
## ssc-let-7d-3p   381   192   198   269   778
## ssc-let-7d-5p  4925  1938  2511  3076  3472
## ssc-let-7e     2811  1302  1463  1690  2512
```

Check the format of the 24 Illumina libraries' config file:


```r
head(illuminaconfig)
```

```
##              inputfile  id
## 1 1034_mapper_input.fa 001
## 2 1058_mapper_input.fa 002
## 3 1080_mapper_input.fa 003
## 4 1096_mapper_input.fa 004
## 5 1116_mapper_input.fa 005
## 6 1134_mapper_input.fa 006
```

Remove the "_mapper_input.fa" from the inputfile column and add it as its own column


```r
filesplit<-strsplit(illuminaconfig$inputfile, "_")
illuminaconfig$filename<-sapply(filesplit, "[", 1)
head(illuminaconfig)
```

```
##              inputfile  id filename
## 1 1034_mapper_input.fa 001     1034
## 2 1058_mapper_input.fa 002     1058
## 3 1080_mapper_input.fa 003     1080
## 4 1096_mapper_input.fa 004     1096
## 5 1116_mapper_input.fa 005     1116
## 6 1134_mapper_input.fa 006     1134
```

```r
str(illuminaconfig)
```

```
## 'data.frame':	24 obs. of  3 variables:
##  $ inputfile: chr  "1034_mapper_input.fa" "1058_mapper_input.fa" "1080_mapper_input.fa" "1096_mapper_input.fa" ...
##  $ id       : chr  "001" "002" "003" "004" ...
##  $ filename : chr  "1034" "1058" "1080" "1096" ...
```

Use that column as an index to subset the no.zero.dfmeanrcround matrix to the 24 libraries.


```r
bioosci24libmirnaexp<-no.zero.dfmeanrcround[ ,illuminaconfig$filename]
dim(bioosci24libmirnaexp)
```

```
## [1] 335  24
```

```r
head(bioosci24libmirnaexp)
```

```
##                1034  1058  1080  1096  1116  1134  1154 1170  1194  1240
## ssc-let-7a    48132 40758 38799 35977 36678 59997 49982 2366 30864 35912
## ssc-let-7c    32745 34313 28022 20772 29022 45232 37824 3077 17962 24015
## ssc-let-7d-3p   381   778   774   440   710  1071   750  525   256   269
## ssc-let-7d-5p  4925  3472  3705  3259  3083  6661  5225  175  2588  3037
## ssc-let-7e     2811  2512  2229  1898  2069  3942  2481   75  1716  1752
## ssc-let-7f    35432 16512 19175 29076 13446 34402 29411  576 20800 27197
##                1278  1300  1426  1434  1458  1484  1502  1512  1534  1580
## ssc-let-7a    15876 28829 33168 28218 28524 32304 47522 20253 27358 35516
## ssc-let-7c    14181 21753 20260 18147 17704 24491 28812 17415 16277 25721
## ssc-let-7d-3p  1033   652   170   212   147   566   383   979   180   492
## ssc-let-7d-5p  1356  2545  3037  2058  2518  2799  4754  1831  2315  3304
## ssc-let-7e      557  1383  1968  1413  1225  1786  3047   996  1532  1880
## ssc-let-7f    10064 18572 24686 20088 20538 12570 35356 10026 20536 19231
##                1594  1640  1644  1662
## ssc-let-7a    36868 37179 38866 39304
## ssc-let-7c    22762 19671 24475 25909
## ssc-let-7d-3p   334   271   306   417
## ssc-let-7d-5p  3616  3814  3651  3820
## ssc-let-7e     2096  2017  2711  2262
## ssc-let-7f    27982 30385 27878 29883
```

```r
colnames(bioosci24libmirnaexp)
```

```
##  [1] "1034" "1058" "1080" "1096" "1116" "1134" "1154" "1170" "1194" "1240"
## [11] "1278" "1300" "1426" "1434" "1458" "1484" "1502" "1512" "1534" "1580"
## [21] "1594" "1640" "1644" "1662"
```

```r
sum(colnames(bioosci24libmirnaexp) != illuminaconfig$filename)
```

```
## [1] 0
```

```r
sum(rownames(bioosci24libmirnaexp) != rownames(no.zero.dfmeanrcround))
```

```
## [1] 0
```

Check that the data is the same between the subset and the no.zero.dfmeanrcround dataset:


```r
for (i in colnames(bioosci24libmirnaexp)){
	print(all.equal(bioosci24libmirnaexp[,i], no.zero.dfmeanrcround[,i]))
}
```

```
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
## [1] TRUE
```

## Save data


```r
save(bioosci24libmirnaexp, file = "../1_24bioosci_mirna_expression.Rdata")
```

