#' **Script:** `3_24illumina_create_dge_normalize_mirna_exp.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/2_24illumina_mirna_expression_characterization/scripts`
#' 
#' **Date:**  8/18/16
#' 
#' **Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/2_24illumina_mirna_expression_characterization/`
#' 
#' **Input File(s):** 
#' 1. `1_24illumina_filtered_rounded_mean_mature_mirna_exp.Rdata`
#' 2. `2_24illumina_mature_mirna_annotation.Rdata`
#'
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/2_24illumina_mirna_expression_characterization/`
#' 
#' **Output File(s):** `4_24illumina_dge_normalized_mirna_expression.Rdata`
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
#' The objective of this script is to create a dge object of the subset of 24 Illumina libraries. 
#' Additionally, the read counts of the miRNAs will be normalized using the cpm function of edgeR, filtered for expression (rough filter: < 1cpm in > 6 libraries removed) 
#' and the calcNormFactors and estimateCommonDisp functions of edgeR will be applied to the non-cpm read counts prior to the final cpm normalization of the read counts. 
#' 
#' THIS ANALYSIS COMPLETED USING R/3.2.0

#' ## Install libraries
library(limma)
library(edgeR)

rm(list=ls())

setwd("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/2_24illumina_mirna_expression_characterization/scripts/")

#' ## Load data
#' Load the miRDeep2 read count data for the 24 Illumina libraries
load("../1_24illumina_filtered_rounded_mean_mature_mirna_exp.Rdata")


#' Load the annotation file from the mature miRNAs
load("../2_24illumina_mature_mirna_annotation.Rdata")

#' 
#' ## Analysis
#' 
#' ### Create the dge object and transform the read counts: log-cpm, then filter genes by expression (1 cpm in 6 or more samples retained):
#' 

#' Create the dge object:
dge<-DGEList(counts=illumina24.no.zero.dfmeanrcround,genes=illumina24.total.mature.annot2)
dim(dge)
dge[1:5,1:5]

#' Calculate the read counts per million in order to filter miRNAs by normalized expression:
cpm.dge<-cpm(dge)
dim(cpm.dge)
cpm.dge[1:5,1:5]

if (sum(rownames(illumina24.no.zero.dfmeanrcround)!=rownames(cpm.dge))!=0) stop ("miRNAs not the same between read counts and cpm")
if (sum(colnames(illumina24.no.zero.dfmeanrcround)!=colnames(cpm.dge))!=0) stop ("animal ids not the same between read counts and cpm")
if (sum(rownames(illumina24.no.zero.dfmeanrcround)!=rownames(illumina24.total.mature.annot2))!=0) stop ("miRNAs not the same between read counts and annotation")

#' Filter miRNAs with at least 1 cpm in at least 1/4 of the samples (24/4=6)
filtercpm<-rowSums(cpm.dge>=1)>=6
sum(filtercpm)

nrow(cpm.dge) - sum(filtercpm)

#' We are removing 59 miRNA profiles from the analysis
#' 
#' So, keep the miRNA profiles in dge based on those retained in the cpm-filtering step:
#' 
#' This retains the rounded, filtered mean read counts, not the cpm (this will be done later):

dge<-dge[filtercpm,]
names(dge)
dge[1:5,1:5]
dim(dge$counts)

if (sum(colnames(dge)!=colnames(cpm.dge))!=0) stop ("colnames not the same between dge and cpm.dge")

#' Apply the TMM normalization:
dge<-calcNormFactors(dge)
head(dge$samples)
hist(dge$samples$norm.factors)

#' This function (estimateCommonDisp) applies normalization factors, caluclates normalized expression based on robust count of normalized reads.
dge<-estimateCommonDisp(dge,verbose=TRUE)
dge$common.dispersion

#' The cpm conversion is used to normalize the read counts between libraries (accounts for differences in sequencing depth)
illumina24.dge.filter.cpm<-cpm(dge)
head(illumina24.dge.filter.cpm)
dim(illumina24.dge.filter.cpm)


#' ## Save data
save(illumina24.dge.filter.cpm, file="../4_24illumina_dge_normalized_mirna_expression.Rdata")