#' **Script:** `1_extract_24bioosci_libraries_for_quantification.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/3_24bioosci_mirna_expression_characterization/scripts/`
#' 
#' **Date:**  8/17/16
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/2_create_expression_matrix_of_known_mirna/`
#' 2. `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/3_build_dge_object_for_eqtl/`
#' 3. `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/6_prepped_fasta_24illumina_output/`
#' 
#' **Input File(s):** 
#' 
#' 1. `1_exp_filtered_rounded_mean_mature_mirna_expression.Rdata`
#' 2. `2_mature_mirna_annotation.Rdata`
#' 3. `1_config_24illumina_mapper.txt`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/3_24bioosci_mirna_expression_characterization/`
#' 
#' **Output File(s):** `1_24bioosci_mirna_expression.Rdata`
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
#' 
#' The objective of this script is to extract the 24 libraries from the dataset of 174 Bioo Scientific-prepped libraries matching the 24 Illumina-prepped libraries for the PCR duplication analysis.
#' 
#' The libraries will be extracted from the expression matrix that was used to create the dge object for the 174 libraries, output from the miRNA quantification using miRDeep2. 
#' The read counts will then be normalized considering only the 24 libraries as a dataset. 
#' 
#' THIS ANALYSIS COMPLETED WITH R/3.2.0


#' ## Install libraries
library(limma)
library(edgeR)

rm(list=ls())

setwd("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/3_24bioosci_mirna_expression_characterization/scripts/")
#' ## Load data
#' 
#' Load the expression matrix from the 174 Bioo Scientific-prepped data:
load("../../../2_mirna_characterization_expression/2_create_expression_matrix_of_known_mirna/1_exp_filtered_rounded_mean_mature_mirna_expression.Rdata")

#' Load the annotation file for the miRNAs in the dataset:
load("../../../2_mirna_characterization_expression/3_build_dge_object_for_eqtl/2_mature_mirna_annotation.Rdata")

#' Load the config file from the 24 Illumina-prepped libraries:
illuminaconfig<-read.table("../../1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/6_prepped_fasta_24illumina_output/1_config_24illumina_mapper.txt", colClasses=c("character","character"), col.names=c("inputfile", "id"))
head(illuminaconfig)

#' ## Analysis
#' 
#' Check the format of the expression matrix from the 174 Bioo Scientific libraries
no.zero.dfmeanrcround[1:5,1:5]

#' Check the format of the 24 Illumina libraries' config file:
head(illuminaconfig)

#' Remove the "_mapper_input.fa" from the inputfile column and add it as its own column
filesplit<-strsplit(illuminaconfig$inputfile, "_")
illuminaconfig$filename<-sapply(filesplit, "[", 1)
head(illuminaconfig)
str(illuminaconfig)

#' Use that column as an index to subset the no.zero.dfmeanrcround matrix to the 24 libraries.
bioosci24libmirnaexp<-no.zero.dfmeanrcround[ ,illuminaconfig$filename]
dim(bioosci24libmirnaexp)
head(bioosci24libmirnaexp)
colnames(bioosci24libmirnaexp)
sum(colnames(bioosci24libmirnaexp) != illuminaconfig$filename)
sum(rownames(bioosci24libmirnaexp) != rownames(no.zero.dfmeanrcround))

#' Check that the data is the same between the subset and the no.zero.dfmeanrcround dataset:
for (i in colnames(bioosci24libmirnaexp)){
	print(all.equal(bioosci24libmirnaexp[,i], no.zero.dfmeanrcround[,i]))
}

#' ## Save data
save(bioosci24libmirnaexp, file = "../1_24bioosci_mirna_expression.Rdata")
