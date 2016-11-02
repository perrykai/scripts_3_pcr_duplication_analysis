#' **Script:** `1_create_config_file.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/scripts`
#' 
#' **Date:**  `9/20/16`
#' 
#' **Input File Directory:**  `../../1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/6_prepped_fasta_24illumina_output/`
#' 
#' **Input File(s):** `1_config_24illumina_mapper.txt`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/`
#' 
#' **Output File(s):** `1_config_24bioosci_PCRDUP_mapper.txt`
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
#' The objective of this script is to create the config file to feed the 24 PCR duplicate-containing Bioo Scientific - prepped libraries into the miRDeep2 mapper module.
#' 
#' ## Install libraries
#' 
#' ## Load data
rm(list=ls())
setwd("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/scripts/")

fn<- read.table("../../1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/6_prepped_fasta_24illumina_output/1_config_24illumina_mapper.txt")

#' ## Analysis
fn<-sub("mapper_input.fa", "24bioosci_PCRDUP_mapper_input.fa", fn$V1)
head(fn)

length(fn)

numf<-seq(1,length(fn),1)           #create a vector of numbers the same length as the number of file names
numf

dig<-sprintf("%03d",numf)           #sprintf function allows for the creation of a 3-digit numeric identifier
dig

config<-cbind(fn,dig)               #bind the two columns together
head(config)

#' ## Save data
write.table(config, file = "../1_mirdeep2_mapper_input_24bioosci_PCRDUP_samples/1_config_24bioosci_PCRDUP_mapper.txt", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE) #write the config object to a text file