#' **Script:** `1_create_24illumina_mirna_expression_df.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/2_24illumina_mirna_expression_characterization/scripts`
#' 
#' **Date:**  8/17/16 #UPDATED 8/25/16
#' 
#' **Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/9_mirdeep2_core_quantify_predict_24illumina_output`
#' 
#' **Input File(s):** `miRNAs_expressed_all_samples_17_08_2016_t_11_33_33.csv`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/2_24illumina_mirna_expression_characterization`
#' 
#' **Output File(s):** `1_24illumina_rounded_mean_mature_mirna_exp.Rdata`
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 6. [Save data](#save-data)
#' 
#' ## Objectives
#' 
#' The objective of this script is to create a matrix of miRNA expression data (read counts) for incorporation into the PCR duplicate analysis, beginning with the known miRNA expression data output from the miRDeep2 core module.
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
#' 5. Round the mean read counts to the nearest integer for use with the edgeR package.

#' ## Install libraries
library(plyr)


setwd("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/2_24illumina_mirna_expression_characterization/scripts")
rm(list=ls())

#' ## Load data

#' ### 1. Read in the raw read counts output from miRDeep2 quantifier module
rawcounts<-read.csv("../../1_preprocess_24_illumina_samples/9_mirdeep2_core_quantify_predict_24illumina_output/miRNAs_expressed_all_samples_17_08_2016_t_11_33_33.csv", sep = "\t", header = TRUE, row.names=NULL)
head(rawcounts)
dim(rawcounts)

#' Set the name of the first column to 'miRNA'
colnames(rawcounts)[[1]] <- "miRNA"

#' Remove the "X" character from the beginning of each column:
colnames(rawcounts)<-gsub("X", "", colnames(rawcounts))
colnames(rawcounts)

#' View the data set:
head(rawcounts[1:8])

colclass<-NULL

for (i in colnames(rawcounts)) {
   colclass<-c(colclass,class(rawcounts[,i]))
}

table(colclass)
head(colclass)

#' Notice columns 1 and 3 are factors - the mature miRNA names and the miRNA precursor names.

#' ### 2. Read in the config file, maintaining the characters in the 3-digit code names
configfile<-read.table("../../1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/6_prepped_fasta_24illumina_output/1_config_24illumina_mapper.txt", header = FALSE, sep = " ", row.names=NULL, colClasses = c('character','character'))
head(configfile)

#' Remove the "_mapper_input.fa" from each file name to leave the pig id:
configfile$V1<-gsub("_mapper_input.fa", "", configfile$V1)
head(configfile)

colnames(configfile)<-c("pigid", "code")
head(configfile)

#' ## Analysis

#' ### 1. Extract the columns of mature read counts for each miRNA for each animal
mirquant<-rawcounts[,c(1,5:28)]
colnames(mirquant)

dim(mirquant)

head(mirquant[1:8])

#' Take a subset of this data.frame for testing:
test <- mirquant[1:20,1:8]

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

illumina24.dfmeanrc<-ldply(meanrc, fun=NULL, id=names(meanrc))

head(illumina24.dfmeanrc[1:8])

dim(illumina24.dfmeanrc)
#' These dimensions are what would be expected, because there are 411 mature sus scrofa miRNA sequences in miRBase,
#' and there are 24 animals in the analysis, plus the miRNA column.
#' 
#' Check that the correct miRNA name went with the correct data:
if (sum(names(meanrc)!=illumina24.dfmeanrc[,1]) != 0) stop ("miRNA names are not the same")

colnames(illumina24.dfmeanrc)[[1]]<-"miRNA"

if (sum(colnames(illumina24.dfmeanrc)!=colnames(mirquant)) != 0) stop ("animal order not the same")

head(illumina24.dfmeanrc[,1:10])

#' Set first column of dfmeanrc (miRNA ids) as the row.names:
rownames(illumina24.dfmeanrc)<-illumina24.dfmeanrc$miRNA

#' Eliminate column of row names:
illumina24.dfmeanrc<-illumina24.dfmeanrc[,-c(1)]
head(illumina24.dfmeanrc[,1:10])
dim(illumina24.dfmeanrc)

#' 
#' How many miRNAs are not expressed (total read count across all libraries = 0)?
head(rowSums(illumina24.dfmeanrc))
tail(rowSums(illumina24.dfmeanrc))
table(rowSums(illumina24.dfmeanrc)==0)

#' So, 86 miRNA profiles contain 0 read counts total, meaning 0 expression
#' 

#' ### 4. Restore the pig IDs in place of the 3-digit codes as the column names of the data frame, for use with the gblup function of gwaR

head(configfile)

#' Now I need to substitute the 3-digit code with the pig IDs, ensuring the names stay in the correct order:
#' 
#' Use match function to find positional index and match column names:
#' 
#' The object dfmeanrc has column names that need to be re-named. I have the config file which contains
#' the current column names and the desired column names. What I am doing in this code is re-ordering the config file 
#' based on where the config file "code" column matches the position of the dfmeanrc object's column names, then having it return the corresponding value in column "pigid". 
#' 
#' So, when using match, need to have the first argument be the matrix/dataframe you want to change or match, and the second argument be what you want to index it by or match it against. 
#' 
#' "Where does [vector] match in [matrix]?" or "Match the column names of illumina24.dfmeanrc to the configfile "code" column, then return the corresponding pigid."
configfile[match(colnames(illumina24.dfmeanrc),configfile$code),"pigid"]

#' Assign the column names using match:
colnames(illumina24.dfmeanrc)<- configfile[match(colnames(illumina24.dfmeanrc),configfile$code),"pigid"]
head(illumina24.dfmeanrc[1:10])
dim(illumina24.dfmeanrc)

if (sum(colnames(illumina24.dfmeanrc)!=(configfile$pigid))!=0) stop ("match function did not work correctly")

#' ### 5. Round the mean read counts to the nearest integer for use with the voom function
illumina24.dfmeanrc[1:10,1:6]

illumina24.dfmeanrcround<-round(illumina24.dfmeanrc)

illumina24.dfmeanrcround[1:10,1:6]

#' The final matrix of rounded, mean read counts needs to have miRNA rownames and Animal ID colnames
head(rownames(illumina24.dfmeanrcround))

head(colnames(illumina24.dfmeanrcround))

if (sum(colnames(illumina24.dfmeanrcround)!= configfile$pigid)!= 0) stop ("rownames do not match pigid")


#' ## Save data
#' What I am saving here is the rounded, average read counts in an .Rdata object
save(illumina24.dfmeanrcround, file = "../1_24illumina_rounded_mean_mature_mirna_exp.Rdata")
