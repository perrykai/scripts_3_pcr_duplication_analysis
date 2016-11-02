#' **Script:** `1_join_24Bioosci_PCRDUP_noPCRdup_datasets_dge_DE_analysis.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/8_DE_analysis_BiooNextFlex_PCRDup_noPCRDup/scripts`
#' 
#' **Date:**  9/27/16
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/3_24bioosci_mirna_expression_characterization`
#' 2. `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/4_24bioosci_PCRDUP_mirna_expression_matrix`
#' 3. `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/4_join_datasets_dge_object`
#' 
#' **Input File(s):**
#' 
#' 1. `1_24bioosci_rounded_mean_mature_mirna_expression.Rdata`
#' 2. `1_24bioosci_PCRDUP_rounded_mean_mature_mirna_expression.Rdata`
#' 3. `1_mature_mirna_annotation.Rdata`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/8_DE_analysis_BiooNextFlex_PCRDup_noPCRDup`
#' 
#' **Output File(s):** 

#' 1. `1_24bioosci_PCRDUP_noPCRDup_joint_mirna_expression_df.Rdata`
#' 2. `2_24bioosci_PCRDUP_noPCRDup_joint_mirna_dge_object.Rdata`
#' 3. `3_24bioosci_PCRDUP_noPCRDup_joint_mirna_voom_deanalysis_results.Rdata`
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
#' and the calcNormFactors and estimateCommonDisp functions of edgeR will be applied to the non-cpm read counts prior to the differential expression analysis of the read counts using the voom function. 
#' 
#' THIS ANALYSIS COMPLETED WITH R/3.2.0

#' ## Install libraries
library(methods)
library(limma)
library(edgeR)
library(statmod)
library(qvalue)

rm(list=ls())

setwd("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/8_DE_analysis_BiooNextFlex_PCRDup_noPCRDup/scripts")

#' ## Load data
#' Load the miRDeep2 read count data for the 24 Bioo Scientific libraries without PCR duplicates
load("../../3_24bioosci_mirna_expression_characterization/1_24bioosci_rounded_mean_mature_mirna_expression.Rdata")

#' Load the miRDeep2 read count data for the 24 Bioo Scientific libraries containing PCR DUPLICATES
load("../../6_process_24bioosci_PCRDUP_samples/4_24bioosci_PCRDUP_mirna_expression_matrix/1_24bioosci_PCRDUP_rounded_mean_mature_mirna_expression.Rdata")

#' Load the annotation file for the 24 Illumina libraries 
load("../../4_join_datasets_dge_object/1_mature_mirna_annotation.Rdata")
ls()

#' ## Analysis
#' 
#' ### 1. Join the two expression datasets (rounded, mean mature miRNA expression)
#' 
#' Make a matrix of miRNA read counts from the Bioo PCR duplicate-removed libraries
nopcrdup.dfmeanrcround<-as.matrix(bioosci24libmirnaexp)
dim(nopcrdup.dfmeanrcround)
#' Make the sample IDs more descriptive (BNFnoPCRDup = PCR Duplicate-removed libraries)
colnames(nopcrdup.dfmeanrcround)<-paste(colnames(nopcrdup.dfmeanrcround), ".BNFnoPCRDup", sep = "")
nopcrdup.dfmeanrcround[1:5,1:5]

#' Define which miRNAs have zero expression in the Bioo PCR-dup-removed libraries:
nopcrdupzeromir<-rownames(nopcrdup.dfmeanrcround[rowSums(nopcrdup.dfmeanrcround)==0,])
length(nopcrdupzeromir)
nopcrdupzeromir


#' Make a matrix of miRNA read counts from the Bioo prepped libraries
pcrdup.bioosci24libmirnaexp<-as.matrix(pcrdup.bioo24.dfmeanrcround)
dim(pcrdup.bioosci24libmirnaexp)
#' Make the sample IDs more descriptive (BNFPCRDup = Bioo Scientific Next Flex WITH PCR duplicates)
colnames(pcrdup.bioosci24libmirnaexp)<-paste(colnames(pcrdup.bioosci24libmirnaexp), ".BNFPCRDup", sep = "")
pcrdup.bioosci24libmirnaexp[1:5,1:5]
#' Define which miRNAs have zero expression in the Bioo prepped libraries:
pcrdup.bioozeromir<-rownames(pcrdup.bioosci24libmirnaexp[rowSums(pcrdup.bioosci24libmirnaexp)==0,])
length(pcrdup.bioozeromir)
pcrdup.bioozeromir

#' In total, 85 miRNAs are not expressed in either dataset:
sum(nopcrdupzeromir %in% pcrdup.bioozeromir)

#' miRNAs not expressed in either dataset:
nopcrdupzeromir[nopcrdupzeromir %in% pcrdup.bioozeromir]

#' Check that the rownames are equal between the datasets
if (sum(rownames(pcrdup.bioosci24libmirnaexp) != rownames(nopcrdup.dfmeanrcround)) != 0) stop ("rownames not equal between datasets")

#' Join the two datasets using cbind
pcrdup.jointmirnaexp<-cbind(nopcrdup.dfmeanrcround,pcrdup.bioosci24libmirnaexp)
pcrdup.jointmirnaexp[1:10,1:30]
dim(pcrdup.jointmirnaexp)
str(pcrdup.jointmirnaexp)

#' Check that the data remained the same between the datasets (no shuffling took place):
for (i in colnames(nopcrdup.dfmeanrcround)){
        print(sum(nopcrdup.dfmeanrcround[,i] != pcrdup.jointmirnaexp[,i]))
}

for (i in colnames(pcrdup.bioosci24libmirnaexp)){
        print(sum(pcrdup.bioosci24libmirnaexp[,i] != pcrdup.jointmirnaexp[,i]))
}

if (sum(rownames(nopcrdup.dfmeanrcround) != rownames(pcrdup.jointmirnaexp)) != 0) stop ("rownames not equal between datasets")
if (sum(rownames(pcrdup.bioosci24libmirnaexp) != rownames(pcrdup.jointmirnaexp)) != 0) stop ("rownames not equal between datasets")




#' ### 2. Create dge object for the joint miRNA expression dataset, including a group option with the library preparation kit determining the groups.
pcrdup.jointmirnadge<-DGEList(counts=pcrdup.jointmirnaexp, genes=illumina24.total.mature.annot2)
names(pcrdup.jointmirnadge)
dim(pcrdup.jointmirnadge$counts)
head(pcrdup.jointmirnadge$counts)
pcrdup.jointmirnadge$samples


#' Check that the data remained the same between the datasets (no shuffling took place):
if (sum(rownames(nopcrdup.dfmeanrcround) != rownames(pcrdup.jointmirnadge)) != 0) stop ("rownames not equal between datasets")
if (sum(rownames(pcrdup.bioosci24libmirnaexp) != rownames(pcrdup.jointmirnadge)) != 0) stop ("rownames not equal between datasets")

for (i in colnames(nopcrdup.dfmeanrcround)){
        print(sum(nopcrdup.dfmeanrcround[,i] != pcrdup.jointmirnadge$counts[,i]))
}

for (i in colnames(pcrdup.bioosci24libmirnaexp)){
        print(sum(pcrdup.bioosci24libmirnaexp[,i] != pcrdup.jointmirnadge$counts[,i]))
}


#' ### 3. Filter the joint expression dataset (dge object)

#' First, eliminate miRNAs whose total expression is 0 among the datasets
sum(rowSums(pcrdup.jointmirnadge$counts) == 0)

#' So, 85 miRNAs have 0 expression among the datasets

pcrdup.jointmirnadge<-pcrdup.jointmirnadge[rowSums(pcrdup.jointmirnadge$counts)>0,]
dim(pcrdup.jointmirnadge)

#' Calculate the read counts per million in order to filter miRNAs by normalized expression:
cpm.pcrdup.jointmirnadge<-cpm(pcrdup.jointmirnadge)
dim(cpm.pcrdup.jointmirnadge)
cpm.pcrdup.jointmirnadge[1:5,1:5]

if (sum(colnames(nopcrdup.dfmeanrcround)!=colnames(cpm.pcrdup.jointmirnadge[1:24]))!=0) stop ("animal ids not the same between read counts and cpm")
if (sum(rownames(nopcrdup.dfmeanrcround)!=rownames(illumina24.total.mature.annot2))!=0) stop ("miRNAs not the same between read counts and annotation")

if (sum(colnames(pcrdup.bioosci24libmirnaexp)!=colnames(cpm.pcrdup.jointmirnadge[25:ncol(cpm.pcrdup.jointmirnadge)]))!=0) stop ("animal ids not the same between read counts and cpm")
if (sum(rownames(pcrdup.bioosci24libmirnaexp)!=rownames(illumina24.total.mature.annot2))!=0) stop ("miRNAs not the same between read counts and annotation")

#' Filter miRNAs with at least 1 cpm in at least 1/4 of the samples (24/4=6)
filtercpm<-rowSums(cpm.pcrdup.jointmirnadge>=1)>=6
sum(filtercpm)

nrow(cpm.pcrdup.jointmirnadge) - sum(filtercpm)

#' We are removing 26 miRNA profiles from the analysis
#' 
#' So, keep the miRNA profiles in dge based on those retained in the cpm-filtering step:
#' 
#' This retains the rounded, mean read counts, not the cpm (this will be done later):

pcrdup.jointmirnadge<-pcrdup.jointmirnadge[filtercpm,]
names(pcrdup.jointmirnadge)
pcrdup.jointmirnadge[1:5,1:5]
dim(pcrdup.jointmirnadge$counts)

if (sum(colnames(pcrdup.jointmirnadge)!=colnames(cpm.pcrdup.jointmirnadge))!=0) stop ("colnames not the same between dge and cpm.pcrdup.jointmirnadge")

#' Apply the TMM normalization (normalizes for RNA composition by finding a set of scaling factors for the library sizes that minimize the log-fold changes between the samples for the most genes):
#' 
#' The result is the effective library size, which is equal to the product of the original library size and the scaling factor. Effective library size is what is used in downstream analyses.
pcrdup.jointmirnadge<-calcNormFactors(pcrdup.jointmirnadge)
pcrdup.jointmirnadge$samples
hist(pcrdup.jointmirnadge$samples$norm.factors)
pcrdup.jointmirnadge$samples$norm.factors
#' This function (estimateCommonDisp) applies normalization factors, caluclates normalized expression based on robust count of normalized reads.
#' 
#' estimateCommonDisp estimates the quantile-adjusted conditional maximum likelihood (qCML) common dispersion.
#' 
#' This creates a matrix of pseudo-counts, used internally to speed up computation of conditional likelihood used for dispersion estimation & exact tests in the limma pipeline.
#' Pseudo-counts represent the equivalent counts that would have been observed had the library sizes all been equal, assuming the fitted model.
#' DO NOT INTERPRET PSEUDO-COUNTS AS GENERAL-PURPOSE NORMALIZED COUNTS.
pcrdup.jointmirnadge<-estimateCommonDisp(pcrdup.jointmirnadge,verbose=TRUE)
pcrdup.jointmirnadge$common.dispersion


#' ### 4. Perform the DE analysis (following the case study in edgeR User Guide: "edgeR: differential expression analysis of digital gene expression data")
#' 
#' Apply the voom transformation and build the design matrix
#' 
#' First, create factors of individual and kit to create the design matrix:
indiv<-factor(substr(colnames(pcrdup.jointmirnadge$counts), 1, 4))
head(indiv)
str(indiv)

kit<-factor(rep(c("BNFnoPCRDup", "BNFPCRDup"), c(24,24)))
kit
str(kit)

#' View this as a data.frame to check for correctness:
data.frame(Sample=colnames(pcrdup.jointmirnadge$counts), indiv, kit)

#' Use the model.matrix() command to create the design matrix for this analysis:
design <- model.matrix(~indiv+kit)
rownames(design)<-colnames(pcrdup.jointmirnadge$counts)
design

#' Apply the voom transformation to the dataset
v <- voom(pcrdup.jointmirnadge, design=design, plot=TRUE)
names(v)
dim(v$E)
v$E[1:5,1:5]
dim(v$genes)
v$genes[1:5,1:5]
dim(v$weights)
v$weights[1:5,1:5]
dim(v$targets)
v$targets[1:5,]
dim(v$design)
v$design[1:5,]

#' Now be super paranoid that no NAs or missing data exists...
sum(is.na(v$E))
sum(is.na(v$weights))
sum(v$E == "NA")
sum(v$weights == "NA")
sum(v$E == "NaN")
sum(v$weights == "NaN")
#' ... And we're good! Carry on.

#' Examine the samples for outliers and for other relationships by using the function plotMDS():
#' 
#' This function draws a multi-dimensional scaling plot of the RNA samples in which distances correspond to leading log-fold-changes between each pair of RNA samples.
#' The leading log-fold-change is the average (root-mean-square) of the largest absolute log-fold-changes between each pair of samples. Can be viewed as unsupervised clustering.
plotMDS(v)

#' ### 4. Perform the DE analysis
fit<-lmFit(v, design=design)
names(fit)

#' Residual standard deviations for each gene:
head(fit$sigma)
summary(fit$sigma)

#' Compute moderated t-statistics, moderated F-statistics, and log-odds of DE by empirical Bayes moderation of standard errors towards a common value:
fitmodt<-eBayes(fit)
names(fitmodt)

#' Summarize the differentially expressed genes between the sequencing kits using tobTable.
#' The argument "adjust.method = "BH" means using the FDR method to correct for multiple testing.
toptab<-(topTable(fitmodt, coef=ncol(design), n=Inf, adjust.method="BH"))
dim(toptab)
head(toptab)
rownames(toptab)

sum(toptab$adj.P.Val <= 0.05)
sum(toptab$adj.P.Val <= 0.01)

#' The eBayes function is giving a warning message that the var.prior cannot be estimated. 
#' Investigate manually calculating the unmoderated ("ordinary") t-statistic from the lmFit cor.results and compare the cor.results to the eBayes moderated t-statistics p-values (var.prior = default value)
#' 
#' Manual calculation of t-statistics from lmFit:
ordt<-fit$coef/fit$stdev.unscaled/fit$sigma
rownames(ordt)
head(ordt)

#' Sort the t-statistics to match the order of the toptab rows and extract the kitBNFPCRDup column
ordtkitBNFPCRDup<-ordt[rownames(toptab),"kitBNFPCRDup"]
head(ordtkitBNFPCRDup)

#' Make sure the rownames match between the two datasets:
sum(rownames(ordtkitBNFPCRDup) != rownames(toptab))

#' Compute the p-value using the tail area probability function (make sure to multiply by 2 to get both tails, df= n-p = 48-25 = 23)
#' 
#' To get the correct p-values for all t-statistics doing a two-sided test, take the negative of the absolute value of the t-statistics
ordtpval<-2*pt(-abs(ordtkitBNFPCRDup), df=23)
head(ordtpval)

#' Make sure the rownames match between the two datasets:
sum(rownames(ordtpval) != rownames(ordtkitBNFPCRDup))
sum(rownames(ordtpval) != rownames(toptab))

#' Join the ordinary t-statistic and their associated p-values to the toptab object.
toptab$ord.t.BNFPCRDup<-ordtkitBNFPCRDup
toptab$ord.t.BNFPCRDup.pval<-ordtpval

head(toptab)

#' Plot the -log10 p-values from the test statistics, both moderated and unmoderated
plot(-log(toptab$ord.t.BNFPCRDup.pval), -log(toptab$P.Value))
abline(0,1)

#' The cor.results are highly correlated
cor.test(-log(toptab$ord.t.BNFPCRDup.pval), -log(toptab$P.Value))

#' So, since the two p-values are so similar, either method can be used. 
#' 
#' For this analysis, we will continue using the eBayes method. 


#' ----------------------------------------------

#' The next step is to compute the gene-wise correlation analysis for the two kits.
#' For this analysis, I will use the voom log2 cpm for each library
#' 
#' Using the cor.test function on the voom log2 cpm objects to obtain gene-wise correlations between samples in the ITS and BFN kits. 
#' 
#' We will use the alternative hypothesis of "greater than", as we expect there would be positive association between the two gene expression datasets.

head(v$E)

itscpm<-v$E[,1:24]
dim(itscpm)
bnfcpm<-v$E[,25:ncol(v$E)]
dim(bnfcpm)

#' Make sure the samples are in the same order in each subset
sum(substr(colnames(itscpm), 1, 4) != substr(colnames(bnfcpm), 1, 4))

#' Create an empty matrix to fill with the correlation results:
cor.results<-matrix(NA, nrow(itscpm), 3)

#' This for loop conducts the gene-wise correlation analysis for the 299 miRNAs between the two kits:
for (i in 1:nrow(itscpm)){
	x<-cor.test(itscpm[i,], bnfcpm[i,], alternative="greater")
	cor.results[i,1]<-x$estimate
	cor.results[i,2]<-x$statistic
	cor.results[i,3]<-x$p.value
}

rownames(cor.results)<-rownames(itscpm)
colnames(cor.results)<-c("cor.estimate", "t.statistic", "p.value")

dim(cor.results)
head(cor.results)
summary(cor.results[,"cor.estimate"])

#' Check a histogram of the correlation coefficients:
hist(cor.results[,"cor.estimate"])

colnames(cor.results)

#' Adjust the p.values from the correlation analysis for multiple tests:
qcor.results<-qvalue(cor.results[,"p.value"])
names(qcor.results)

#' Check the summary of the q-value results to check pi0 (want it close to 0) and the number of significant calls (q<0.05)
summary(qcor.results)

#' Reorder the cor.results dataframe to match the toptable rownames for plotting:
cor.results<-cor.results[rownames(toptab),]
dim(cor.results)
head(cor.results)
#' Plot the average expression of the gene (from the toptable object) vs the -log10 pvalue of the correlation to see the relationship between expression and correlation:
plot(toptab$AveExpr, -log(cor.results[,"p.value"]))
abline(0,1)


#' Plot the expression of each miRNA in the two datasets:
plot(v$E[,1:24], v$E[,25:48], ylab="BNFnoPCRDup", xlab="BNFPCRDup")
abline(0,1)


#' ## Save data
save(pcrdup.jointmirnaexp, file="../1_24illumina_24bioosci_PCRDUP_joint_mirna_expression_df.Rdata")
save(pcrdup.jointmirnadge, file="../2_24illumina_24bioosci_PCRDUP_joint_mirna_dge_object.Rdata")
save(v, fit, toptab, file="../3_24illumina_24bioosci_PCRDUP_joint_mirna_voom_deanalysis_results.Rdata")
save(cor.results, file="../4_24illumina_24bioosci_PCRDUP_joint_mirna_genewise_cor_results.Rdata")