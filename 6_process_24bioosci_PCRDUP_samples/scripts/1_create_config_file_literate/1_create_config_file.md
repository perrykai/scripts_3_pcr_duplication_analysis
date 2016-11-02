**Script:** `1_create_config_file.R`

**Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/scripts`

**Date:**  `9/20/16`

**Input File Directory:**  `../../1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/6_prepped_fasta_24illumina_output/`

**Input File(s):** `1_config_24illumina_mapper.txt`

**Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/`

**Output File(s):** `1_config_24bioosci_PCRDUP_mapper.txt`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Save data](#save-data)

## Objectives

The objective of this script is to create the config file to feed the 24 PCR duplicate-containing Bioo Scientific - prepped libraries into the miRDeep2 mapper module.

## Install libraries

## Load data


```r
rm(list=ls())
setwd("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/scripts/")

fn<- read.table("../../1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/6_prepped_fasta_24illumina_output/1_config_24illumina_mapper.txt")
```

## Analysis


```r
fn<-sub("mapper_input.fa", "24bioosci_PCRDUP_mapper_input.fa", fn$V1)
head(fn)
```

```
## [1] "1034_24bioosci_PCRDUP_mapper_input.fa"
## [2] "1058_24bioosci_PCRDUP_mapper_input.fa"
## [3] "1080_24bioosci_PCRDUP_mapper_input.fa"
## [4] "1096_24bioosci_PCRDUP_mapper_input.fa"
## [5] "1116_24bioosci_PCRDUP_mapper_input.fa"
## [6] "1134_24bioosci_PCRDUP_mapper_input.fa"
```

```r
length(fn)
```

```
## [1] 24
```

```r
numf<-seq(1,length(fn),1)           #create a vector of numbers the same length as the number of file names
numf
```

```
##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
## [24] 24
```

```r
dig<-sprintf("%03d",numf)           #sprintf function allows for the creation of a 3-digit numeric identifier
dig
```

```
##  [1] "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011"
## [12] "012" "013" "014" "015" "016" "017" "018" "019" "020" "021" "022"
## [23] "023" "024"
```

```r
config<-cbind(fn,dig)               #bind the two columns together
head(config)
```

```
##      fn                                      dig  
## [1,] "1034_24bioosci_PCRDUP_mapper_input.fa" "001"
## [2,] "1058_24bioosci_PCRDUP_mapper_input.fa" "002"
## [3,] "1080_24bioosci_PCRDUP_mapper_input.fa" "003"
## [4,] "1096_24bioosci_PCRDUP_mapper_input.fa" "004"
## [5,] "1116_24bioosci_PCRDUP_mapper_input.fa" "005"
## [6,] "1134_24bioosci_PCRDUP_mapper_input.fa" "006"
```

## Save data


```r
write.table(config, file = "../1_mirdeep2_mapper_input_24bioosci_PCRDUP_samples/1_config_24bioosci_PCRDUP_mapper.txt", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE) #write the config object to a text file
```

