#' **Script:** `6_extract_24bioosci_PCRDUP_mirna_expression_df.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/scripts`
#' 
#' **Date:**  9/22/16
#' 
#' **Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/3_mirdeep2_core_output_24bioosci_PCRDUP_samples`
#' 
#' **Input File(s):** `miRNAs_expressed_all_samples_22_09_2016_t_14_48_56.csv`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/4_24bioosci_PCRDUP_mirna_expression_matrix/`
#' 
#' **Output File(s):** `1_24bioosci_PCRDUP_rounded_mean_mature_mirna_expression.Rdata`
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 5. [Save data](#save-data)
#' 
#' ## Objectives
#' The objective of this script is to create a matrix of miRNA expression data (read counts) for incorporation into the miRNA DE analysis, beginning with the known miRNA expression data output from the miRDeep2 core module.
#' Later, expression data of the putative novel miRNA candidates can also be included, also output from the miRDeep2 core module. 
#' To acheive one read count per miRNA per animal, I will need to take the average of the mature read counts in instances of multiple precursor sequences.
#' 
#' The result of this script will be a data frame containing one average read count per miRNA per animal.
#' 
#' So, what I need to do:
#' 
#' 1. Extract the columns of the mature miRNA read counts for each animal
#' 2. Use the 'by' function to apply the function colmeans to the entire data frame of read counts (What this will do is go down the columns looking at the index of grouped miRNA names and take the average of the read counts for that group of miRNA. The result of this will be a list containing the average read counts for each miRNA for each animal)
#' 3. Transform the list output back into a data.frame using the plyr package to prepare for gblup function
#' 4. Restore the pig IDs in place of the 3-digit codes as the column names of the data frame, for use with the gblup function of gwaR
#' 5. Round the mean read counts to the nearest integer for use with the voom function
#' 6. Save the final data.frame.

#' ## Install libraries
rm(list=ls())
library(plyr)

#' ## Load data
setwd("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/scripts")

#' ### 1. Read in the csv of the miRNA expression profiles, using check.names default to correct column names with non-ASCII characters
rc<-read.csv("../3_mirdeep2_core_output_24bioosci_PCRDUP_samples/miRNAs_expressed_all_samples_22_09_2016_t_14_48_56.csv", sep = "\t", header = TRUE, row.names=NULL)

#' Set the name of the first column to "miRNA":
colnames(rc)[[1]]<-"miRNA"

#' Remove the "X" character from the beginning of each column:
colnames(rc)<-gsub("X", "", colnames(rc))
colnames(rc)

#' View the data set:
head(rc[1:8])

#' Do a for loop to check the class of each column in the data.frame

colclass<-NULL

for (i in colnames(rc)) {
   colclass<-c(colclass,class(rc[,i]))
}

table(colclass)
head(colclass)


#' Notice here that the mature miRNA names and the miRNA precursor names are considered factors (columns 1 and 3).

#' ### 2. Read in the config file, maintaining the characters in the 3-digit code names
configfile<-read.table("../1_mirdeep2_mapper_input_24bioosci_PCRDUP_samples/1_config_24bioosci_PCRDUP_mapper.txt", header = FALSE, sep = " ", row.names=NULL, colClasses = c('character','character'))
head(configfile)

#' Remove the "_express.fa" from each file name to leave the pig id:
configfile$V1<-gsub("_24bioosci_PCRDUP_mapper_input.fa", "", configfile$V1)

#' Make filenames more informative:
colnames(configfile)<-c("pigid","code")
colnames(configfile)
head(configfile)

#' ## Analysis

#' ### 1. Extract the columns of mature read counts for each miRNA for each animal
mirquant<-rc[,c(1,5:28)]
colnames(mirquant)

dim(mirquant)

head(mirquant[1:8])


#' Take a subset of this data.frame for testing:
test<-mirquant[1:20,1:8]

test[1:10,]

#' ### 2. Use the 'by' function to apply the function colMeans to the entire data frame of read counts:
#' (What this will do is go down the columns looking at the index of grouped miRNA names and take the average of the read counts for that miRNA)
#' The result of this will be a list containing the average read counts for each miRNA for each animal.
#' 
#' Example: by(data, index, function)
bytst<-by(test[,2:ncol(test)], test[,1], colMeans)
bytst[1:25]
#' Notice here that the rest of the miRNA names remain since the miRNA name is a factor, but since there is no data for them they are filled with NULL.
#' 

#' Apply the by function to the full dataframe:
meanrc<-by(mirquant[,2:ncol(mirquant)], mirquant[,1], colMeans)

#' This should be 411 (the number of mature pig miRNAs in miRBase), meaning we have one expression profile for each mature miRNA:
length(meanrc)

head(meanrc)

#' 
#' ### 3. Transform the list output back into a data.frame using the plyr package to prepare for gblup function:

#' Example: ldply(.data, .fun, .id)
#' 
#' id = name of the index column (used if data is a named list). Pass NULL to avoid creation
#'      of the index column. For compatibility, omit this argument or pass "NA" to avoid converting the index column
#'      to a factor; in this case, ".id" is used as column name.

pcrdup.bioo24.dfmeanrc<-ldply(meanrc, fun=NULL, id=names(meanrc))

head(pcrdup.bioo24.dfmeanrc[1:8])

dim(pcrdup.bioo24.dfmeanrc)
#' These dimensions are what would be expected, because there are 411 mature sus scrofa miRNA sequences in miRBase,
#' and there are 174 animals in the analysis, plus the miRNA column.
#' 
#' Check that the correct miRNA name went with the correct data:
if (sum(names(meanrc)!=pcrdup.bioo24.dfmeanrc[,1]) != 0) stop ("miRNA names are not the same")

colnames(pcrdup.bioo24.dfmeanrc)[[1]]<-"miRNA"

if (sum(colnames(pcrdup.bioo24.dfmeanrc)!=colnames(mirquant)) != 0) stop ("animal order not the same")

head(pcrdup.bioo24.dfmeanrc[,1:10])

#' 
#' Set first column of pcrdup.bioo24.dfmeanrc (miRNA ids) as the row.names:
rownames(pcrdup.bioo24.dfmeanrc)<-pcrdup.bioo24.dfmeanrc$miRNA

#' Eliminate column of row names:
pcrdup.bioo24.dfmeanrc<-pcrdup.bioo24.dfmeanrc[,-c(1)]
head(pcrdup.bioo24.dfmeanrc[,1:10])
dim(pcrdup.bioo24.dfmeanrc)



#' ### 4. Restore the pig IDs in place of the 3-digit codes as the column names of the data frame, for use with the gblup function of gwaR

head(configfile)

#' Now I need to substitute the 3-digit code with the pig IDs, ensuring the names stay in the correct order:
#' 
#' Use match function to find positional index and match column names:
#' 
#' The object pcrdup.bioo24.dfmeanrc has column names that need to be re-named. I have the config file which contains
#' the current column names and the desired column names. What I am doing in this code is re-ordering the config file 
#' based on where the config file "code" column matches the position of the pcrdup.bioo24.dfmeanrc object's column names, then having it return the corresponding value in column "pigid". 
#' 
#' So, when using match, need to have the first argument be the matrix/dataframe you want to change or match, and the second argument be what you want to index it by or match it against. 
#' 
#' "Where does [vector] match in [matrix]?" or "Match the column names of pcrdup.bioo24.dfmeanrc to the configfile "code" column, then return the corresponding pigid."
configfile[match(colnames(pcrdup.bioo24.dfmeanrc),configfile$code),"pigid"]

#' Assign the column names using match:
colnames(pcrdup.bioo24.dfmeanrc)<- configfile[match(colnames(pcrdup.bioo24.dfmeanrc),configfile$code),"pigid"]
head(pcrdup.bioo24.dfmeanrc[1:10])
dim(pcrdup.bioo24.dfmeanrc)

if (sum(colnames(pcrdup.bioo24.dfmeanrc)!=(configfile$pigid))!=0) stop ("match function did not work correctly")

#' ### 5. Round the mean read counts to the nearest integer for use with the voom function
pcrdup.bioo24.dfmeanrc[1:10,1:6]

pcrdup.bioo24.dfmeanrcround<-round(pcrdup.bioo24.dfmeanrc)

pcrdup.bioo24.dfmeanrcround[1:10,1:6]

#' The final matrix of rounded, mean read counts needs to have miRNA rownames and Animal ID colnames
head(rownames(pcrdup.bioo24.dfmeanrcround))

head(colnames(pcrdup.bioo24.dfmeanrcround))

if (sum(colnames(pcrdup.bioo24.dfmeanrcround)!= configfile$pigid)!= 0) stop ("rownames do not match pigid")

#' ## Save data
#' What I am saving here is the rounded, average read counts in an .Rdata object
save(pcrdup.bioo24.dfmeanrcround, file = "../4_24bioosci_PCRDUP_mirna_expression_matrix/1_24bioosci_PCRDUP_rounded_mean_mature_mirna_expression.Rdata")
