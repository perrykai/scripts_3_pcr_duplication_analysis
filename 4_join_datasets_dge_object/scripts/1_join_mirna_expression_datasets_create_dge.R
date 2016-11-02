#' **Script:** `1_join_mirna_expression_datasets_create_dge.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/4_join_datasets_dge_object/scripts`
#' 
#' **Date:**  8/29/16
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/2_24illumina_mirna_expression_characterization`
#' 2. `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/4_join_datasets_dge_object`
#' 
#' **Input File(s):**
#' 
#' 1. `1_24illumina_rounded_mean_mature_mirna_exp.Rdata`
#' 2. `1_24bioosci_rounded_mean_mature_mirna_expression.Rdata`
#' 3. `1_mature_mirna_annotation.Rdata`
#' 
#' **Output File Directory:** ``
#' 
#' **Output File(s):** 

#' 1. `3_joint_mirna_expression_df.Rdata`
#' 2. `4_joint_mirna_dge_object.Rdata`
#' 3. `5_joint_mirna_deanalysis_results.Rdata`
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

#' The objective of this script is to (1) join the two miRNA expression datasets into one and (2) make them into a dge object. 
#' This dge object will then (3) be filtered for miRNA expression: first, the read counts of the miRNAs will be normalized using the cpm function of edgeR, then filtered for expression (rough filter: < 1cpm in > 6 libraries removed) 
#' and the calcNormFactors and estimateCommonDisp functions of edgeR will be applied to the non-cpm read counts prior to the differential expression analysis of the read counts. 
#' 
#' THIS ANALYSIS COMPLETED WITH R/3.2.0

#' ## Install libraries
library(methods)
library(limma)
library(edgeR)
library(statmod)

rm(list=ls())

setwd("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/4_join_datasets_dge_object/scripts")

#' ## Load data
#' Load the miRDeep2 read count data for the 24 Illumina libraries
load("../../2_24illumina_mirna_expression_characterization/1_24illumina_rounded_mean_mature_mirna_exp.Rdata")

#' Load the miRDeep2 read count data for the 24 Bioo Scientific libraries
load("../../3_24bioosci_mirna_expression_characterization/1_24bioosci_rounded_mean_mature_mirna_expression.Rdata")

#' Load the annotation file for the 24 libraries 
load("../1_mature_mirna_annotation.Rdata")
ls()

#' ## Analysis
#' 



#' ### 1. Join the two expression datasets (rounded, mean mature miRNA expression)
#' 
#' Make a matrix of miRNA read counts from the Illumina prepped libraries
illumina24.dfmeanrcround<-as.matrix(illumina24.dfmeanrcround)
dim(illumina24.dfmeanrcround)
#' Make the sample IDs more descriptive (ITS = Illumina TruSeq)
colnames(illumina24.dfmeanrcround)<-paste(colnames(illumina24.dfmeanrcround), ".ITS", sep = "")
illumina24.dfmeanrcround[1:5,1:5]

#' Define which miRNAs have zero expression in the Bioo prepped libraries:
illumina24zeromir<-rownames(illumina24.dfmeanrcround[rowSums(illumina24.dfmeanrcround)==0,])
length(illumina24zeromir)
illumina24zeromir


#' Make a matrix of miRNA read counts from the Bioo prepped libraries
bioosci24libmirnaexp<-as.matrix(bioosci24libmirnaexp)
dim(bioosci24libmirnaexp)
#' Make the sample IDs more descriptive (BNF = Bioo Scientific Next Flex)
colnames(bioosci24libmirnaexp)<-paste(colnames(bioosci24libmirnaexp), ".BNF", sep = "")
bioosci24libmirnaexp[1:5,1:5]
#' Define which miRNAs have zero expression in the Bioo prepped libraries:
bioozeromir<-rownames(bioosci24libmirnaexp[rowSums(bioosci24libmirnaexp)==0,])
length(bioozeromir)
bioozeromir

#' In total, 78 miRNAs are not expressed in either dataset:
sum(illumina24zeromir %in% bioozeromir)

#' miRNAs not expressed in either dataset:
illumina24zeromir[illumina24zeromir %in% bioozeromir]

#' (8) miRNAs not expressed only in the illumina24 dataset:
illumina24zeromir[!illumina24zeromir %in% bioozeromir]

#' (7) miRNAs not expressed only in the bioo dataset:
bioozeromir[!bioozeromir %in% illumina24zeromir]

#' Check that the rownames are equal between the datasets
if (sum(rownames(bioosci24libmirnaexp) != rownames(illumina24.dfmeanrcround)) != 0) stop ("rownames not equal between datasets")

#' Join the two datasets using cbind
jointmirnaexp<-cbind(illumina24.dfmeanrcround,bioosci24libmirnaexp)
jointmirnaexp[1:10,1:30]
dim(jointmirnaexp)
str(jointmirnaexp)

#' Check that the data remained the same between the datasets (no shuffling took place):
for (i in colnames(illumina24.dfmeanrcround)){
        print(sum(illumina24.dfmeanrcround[,i] != jointmirnaexp[,i]))
}

for (i in colnames(bioosci24libmirnaexp)){
        print(sum(bioosci24libmirnaexp[,i] != jointmirnaexp[,i]))
}

if (sum(rownames(illumina24.dfmeanrcround) != rownames(jointmirnaexp)) != 0) stop ("rownames not equal between datasets")
if (sum(rownames(bioosci24libmirnaexp) != rownames(jointmirnaexp)) != 0) stop ("rownames not equal between datasets")




#' ### 2. Create dge object for the joint miRNA expression dataset, including a group option with the library preparation kit determining the groups.
jointmirnadge<-DGEList(counts=jointmirnaexp, genes=illumina24.total.mature.annot2)
names(jointmirnadge)
dim(jointmirnadge$counts)
head(jointmirnadge$counts)
jointmirnadge$samples


#' Check that the data remained the same between the datasets (no shuffling took place):
if (sum(rownames(illumina24.dfmeanrcround) != rownames(jointmirnadge)) != 0) stop ("rownames not equal between datasets")
if (sum(rownames(bioosci24libmirnaexp) != rownames(jointmirnadge)) != 0) stop ("rownames not equal between datasets")

for (i in colnames(illumina24.dfmeanrcround)){
        print(sum(illumina24.dfmeanrcround[,i] != jointmirnadge$counts[,i]))
}

for (i in colnames(bioosci24libmirnaexp)){
        print(sum(bioosci24libmirnaexp[,i] != jointmirnadge$counts[,i]))
}




#' ### 3. Filter the joint expression dataset (dge object)

#' First, eliminate miRNAs whose total expression is 0 among the datasets
sum(rowSums(jointmirnadge$counts) == 0)

#' So, 78 miRNAs have 0 expression among the datasets

jointmirnadge<-jointmirnadge[rowSums(jointmirnadge$counts)>0,]
dim(jointmirnadge)

#' Calculate the read counts per million in order to filter miRNAs by normalized expression:
cpm.jointmirnadge<-cpm(jointmirnadge)
dim(cpm.jointmirnadge)
cpm.jointmirnadge[1:5,1:5]

if (sum(colnames(illumina24.dfmeanrcround)!=colnames(cpm.jointmirnadge[1:24]))!=0) stop ("animal ids not the same between read counts and cpm")
if (sum(rownames(illumina24.dfmeanrcround)!=rownames(illumina24.total.mature.annot2))!=0) stop ("miRNAs not the same between read counts and annotation")

if (sum(colnames(bioosci24libmirnaexp)!=colnames(cpm.jointmirnadge[25:ncol(cpm.jointmirnadge)]))!=0) stop ("animal ids not the same between read counts and cpm")
if (sum(rownames(bioosci24libmirnaexp)!=rownames(illumina24.total.mature.annot2))!=0) stop ("miRNAs not the same between read counts and annotation")

#' Filter miRNAs with at least 1 cpm in at least 1/4 of the samples (24/4=6)
filtercpm<-rowSums(cpm.jointmirnadge>=1)>=6
sum(filtercpm)

nrow(cpm.jointmirnadge) - sum(filtercpm)

#' We are removing 34 miRNA profiles from the analysis
#' 
#' So, keep the miRNA profiles in dge based on those retained in the cpm-filtering step:
#' 
#' This retains the rounded, mean read counts, not the cpm (this will be done later):

jointmirnadge<-jointmirnadge[filtercpm,]
names(jointmirnadge)
jointmirnadge[1:5,1:5]
dim(jointmirnadge$counts)

if (sum(colnames(jointmirnadge)!=colnames(cpm.jointmirnadge))!=0) stop ("colnames not the same between dge and cpm.jointmirnadge")

#' Apply the TMM normalization (normalizes for RNA composition by finding a set of scaling factors for the library sizes that minimize the log-fold changes between the samples for the most genes):
#' 
#' The result is the effective library size, which is equal to the product of the original library size and the scaling factor. Effective library size is what is used in downstream analyses.
jointmirnadge<-calcNormFactors(jointmirnadge)
jointmirnadge$samples
hist(jointmirnadge$samples$norm.factors)
jointmirnadge$samples$norm.factors
#' This function (estimateCommonDisp) applies normalization factors, caluclates normalized expression based on robust count of normalized reads.
#' 
#' estimateCommonDisp estimates the quantile-adjusted conditional maximum likelihood (qCML) common dispersion.
#' 
#' This creates a matrix of pseudo-counts, used internally to speed up computation of conditional likelihood used for dispersion estimation & exact tests in the limma pipeline.
#' Pseudo-counts represent the equivalent counts that would have been observed had the library sizes all been equal, assuming the fitted model.
#' DO NOT INTERPRET PSEUDO-COUNTS AS GENERAL-PURPOSE NORMALIZED COUNTS.
jointmirnadge<-estimateCommonDisp(jointmirnadge,verbose=TRUE)
jointmirnadge$common.dispersion


#' ### 4. Perform the DE analysis (following the case study in edgeR User Guide: "edgeR: differential expression analysis of digital gene expression data")
#' 
#' Examine the samples for outliers and for other relationships by using the function plotMDS():
#' 
#' This function draws a multi-dimensional scaling plot of the RNA samples in which distances correspond to leading log-fold-changes between each pair of RNA samples.
#' The leading log-fold-change is the average (root-mean-square) of the largest absolute log-fold-changes between each pair of samples. Can be viewed as unsupervised clustering.

plotMDS(jointmirnadge)

#' What this shows is that dimension 1 separates the ITS from the BNF samples, while dimension 2 roughly corresponds to sample number. 
#' This shows the paired nature of the samples.

#' 
#' Now to create the design matrix for the proposed hypothesis:
#' Testing for "differential expression" between 1. (PCR-duplciate-removed) BNF and 2. ITS kits, adjusting for differences between individuals.
#'
#' Statistically, this is an additive linear model with individual as a blocking factor

#' 
#' First, create factors of individual and kit to create the design matrix:
indiv<-factor(substr(colnames(jointmirnadge$counts), 1, 4))
head(indiv)
str(indiv)

kit<-factor(rep(c("ITS", "BNF"), c(24,24)))
kit
str(kit)

#' View this as a data.frame to check for correctness:
data.frame(Sample=colnames(jointmirnadge$counts), indiv, kit)

#' Use the model.matrix() command to create the design matrix for this analysis:
design <- model.matrix(~indiv+kit)
rownames(design)<-colnames(jointmirnadge$counts)
design

#' Estimate the dispersion factors for this dataset: 
jointmirnadge<-estimateDisp(jointmirnadge, design=design, robust=TRUE)
jointmirnadge$common.dispersion

#' The square root of the common.dispersion is the coefficient of variation of biological variation.
sqrt(jointmirnadge$common.dispersion)

#' View the dispersion estimates:
plotBCV(jointmirnadge, pch=19, cex=0.6)

#' Now proceed with DE analysis. Fit genewise glms:

fit <- glmFit(jointmirnadge,design)

#' Conduct likelihoo ratio test for ITS vs BNF kits and show the top miRNA
lrtfit<-glmLRT(fit)
topTags(lrtfit)

#' glmLRT() conducted a test for the last coefficient in the linear model, which we can see is kitITS vs kitBFN
colnames(design)

#' Here, the miRNA-wise tests are for ITS vs BFN kit differential expression, adjusting for baseline differences between individuals.
#' The tests can be viewed as analogous to paired t-tests. The top DE tags have very small p-values and FDR values
#' 
#' Take a look at the counts-per-million in individual samples for the top genes:
o <- order(lrtfit$table$PValue)
cpm(jointmirnadge)[o[1:10],]

#' Notice that all the top miRNAs have consistent differences for the individuals.
#' 
#' Total number of DE miRNAs at 1%FDR is given by:
summary(de <- decideTestsDGE(lrtfit, adjust.method="BH", p.value=0.01))

#' Plot log-fold change against log-cpm with DE miRNAs highlighted
detags<-rownames(jointmirnadge)[as.logical(de)]
plotSmear(lrtfit, de.tags=detags, cex=0.6)
abline(h=c(-1,1),col="blue")
#' The blue line indicates 2-fold change
#' 

#' Check the histogram of p-values from the LRT:
y<-lrtfit$table$PValue
hist(y, xlab="LRT p-values")

#' ## Save data
save(jointmirnaexp, file="../3_joint_mirna_expression_df.Rdata")
save(jointmirnadge, file="../4_joint_mirna_dge_object.Rdata")
save(fit, lrtfit, file="../5_joint_mirna_deanalysis_results.Rdata")
