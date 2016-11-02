#' **Script:** `1_voom_de_analysis.R`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/5_voom_de_analysis/scripts`
#' 
#' **Date:**  09/08/16
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/4_join_datasets_dge_object`
#' 
#' **Input File(s):**
#' 
#' 1. `3_joint_mirna_expression_df.Rdata`
#' 2. `1_mature_mirna_annotation.Rdata`
#' 
#' **Output File Directory:** 
#' 
#' 1. `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/5_voom_de_analysis/`
#' 
#' **Output File(s):** 
#' 
#' 1. `1_voom_de_analysis_cor.results.Rdata`
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 5. [Save data](#save-data)

#' ## Objectives
#' The objective of this script is to conduct a differential expression analysis of the 24 libraries sequenced with both the Illumina TruSeq kit and the Bioo Scientific NEXTFlex v2 library prep kit.
#' This analysis will utilize voom and the limma pipeline to compare cor.results with those obtained using edgeR.
#' 
#' ## Install libraries
library(methods)
library(limma)
library(edgeR)
library(qvalue)

sessionInfo()

#' ## Load data
rm(list=ls())

setwd("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/5_voom_de_analysis/scripts")

load("../../4_join_datasets_dge_object/3_joint_mirna_expression_df.Rdata")

load("../../4_join_datasets_dge_object/1_mature_mirna_annotation.Rdata")

#' ## Analysis
dim(jointmirnaexp)
head(jointmirnaexp)

#' ### 1. Create the dge object
jointmirnadge<-DGEList(counts=jointmirnaexp, genes=illumina24.total.mature.annot2)
names(jointmirnadge)
dim(jointmirnadge$counts)
head(jointmirnadge$counts)
jointmirnadge$samples

#' ### 2. Filter the joint expression dataset (dge object)

#' First, eliminate miRNAs whose total expression is 0 among the datasets
sum(rowSums(jointmirnadge$counts) == 0)

#' So, 78 miRNAs have 0 expression among the datasets

jointmirnadge<-jointmirnadge[rowSums(jointmirnadge$counts)>0,]
dim(jointmirnadge)

#' Calculate the read counts per million in order to filter miRNAs by normalized expression:
cpm.jointmirnadge<-cpm(jointmirnadge)
dim(cpm.jointmirnadge)
cpm.jointmirnadge[1:5,1:5]

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
summary(jointmirnadge$samples$norm.factors)

#' This function (estimateCommonDisp) applies normalization factors, caluclates normalized expression based on robust count of normalized reads.
#' 
#' estimateCommonDisp estimates the quantile-adjusted conditional maximum likelihood (qCML) common dispersion.
#' 
#' This creates a matrix of pseudo-counts, used internally to speed up computation of conditional likelihood used for dispersion estimation & exact tests in the limma pipeline.
#' Pseudo-counts represent the equivalent counts that would have been observed had the library sizes all been equal, assuming the fitted model.
#' DO NOT INTERPRET PSEUDO-COUNTS AS GENERAL-PURPOSE NORMALIZED COUNTS.
jointmirnadge<-estimateCommonDisp(jointmirnadge,verbose=TRUE)
jointmirnadge$common.dispersion


newcpm.jointmirnadge<-cpm(jointmirnadge)
head(newcpm.jointmirnadge)



#' ### 3. Apply the voom transformation and build the design matrix
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

#' Apply the voom transformation to the dataset
v <- voom(jointmirnadge, design=design, plot=TRUE)
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
rownames(toptab)[1:5]

sum(toptab$adj.P.Val <= 0.05)
sum(toptab$adj.P.Val <= 0.01)

#' The eBayes function is giving a warning message that the var.prior cannot be estimated. 
#' Investigate manually calculating the unmoderated ("ordinary") t-statistic from the lmFit cor.results and compare the cor.results to the eBayes moderated t-statistics p-values (var.prior = default value)
#' 
#' Manual calculation of t-statistics from lmFit:
ordt<-fit$coef/fit$stdev.unscaled/fit$sigma
rownames(ordt)[1:5]
head(ordt)

#' Sort the t-statistics to match the order of the toptab rows and extract the kitITS column
ordtkitITS<-ordt[rownames(toptab),"kitITS"]
head(ordtkitITS)

#' Make sure the rownames match between the two datasets:
sum(rownames(ordtkitITS) != rownames(toptab))

#' Compute the p-value using the tail area probability function (make sure to multiply by 2 to get both tails, df= n-p = 48-25 = 23)
#' 
#' To get the correct p-values for all t-statistics doing a two-sided test, take the negative of the absolute value of the t-statistics
ordtpval<-2*pt(-abs(ordtkitITS), df=23)
head(ordtpval)

#' Make sure the rownames match between the two datasets:
sum(rownames(ordtpval) != rownames(ordtkitITS))
sum(rownames(ordtpval) != rownames(toptab))

#' Join the ordinary t-statistic and their associated p-values to the toptab object.
toptab$ord.t.ITS<-ordtkitITS
toptab$ord.t.ITS.pval<-ordtpval

head(toptab)

#' Plot the -log10 p-values from the test statistics, both moderated and unmoderated
plot(-log(toptab$ord.t.ITS.pval), -log(toptab$P.Value))
abline(0,1)

#' The cor.results are highly correlated
cor.test(-log(toptab$ord.t.ITS.pval), -log(toptab$P.Value))

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


#' This for loop conducts the individual-wise correlation analysis for the 24 individuals between the two kits:
corind<-NULL

for (i in 1:ncol(itscpm)){
	corind[i]<-cor(itscpm[,i], bnfcpm[,i])
}
corind
summary(corind)

#' Create the vectors of correlation between individuals in the different kits:
#' 
#' 1034.ITS vs all of ITS:
cor1ITS<-cor(itscpm[,"1034.ITS"], itscpm[,2:ncol(itscpm)])
cor1ITS
mean(cor1ITS)
range(cor1ITS)

#' 1034.ITS vs all of BNF:
cor1ITSBNF<-cor(itscpm[,"1034.ITS"], bnfcpm[,2:ncol(bnfcpm)])
cor1ITSBNF
mean(cor1ITSBNF)
range(cor1ITSBNF)

#' 1034.BNF vs all of BNF:
cor1BNF<-cor(bnfcpm[,"1034.BNF"], bnfcpm[,2:ncol(bnfcpm)])
cor1BNF
mean(cor1BNF)
range(cor1BNF)

#' 1034.BNF vs all of ITS:
cor1BNFITS<-cor(bnfcpm[,"1034.BNF"], itscpm[,2:ncol(itscpm)])
cor1BNFITS
mean(cor1BNFITS)
range(cor1BNFITS)

#' Plot the expression of each miRNA in the two datasets:
plot(v$E[,1:24], v$E[,25:48], xlab="kitITS", ylab="kitBNF")
abline(0,1)

#' ----------------------------------------------
#' 
#' Reorder the cpm.jointmirnadge object to compute the correlations:
newcpm.jointmirnadge<-newcpm.jointmirnadge[rownames(toptab),]
nrow(newcpm.jointmirnadge)
nrow(cor.results)
sum(rownames(newcpm.jointmirnadge) != rownames(cor.results))

rowMeans(newcpm.jointmirnadge)[1:5]
rowSums(newcpm.jointmirnadge)[1:5]

rowSums(newcpm.jointmirnadge)[1:5]/ncol(newcpm.jointmirnadge)

head(newcpm.jointmirnadge)[1:5, 1:5]
newitscpm<-newcpm.jointmirnadge[,1:24]
newbnfcpm<-newcpm.jointmirnadge[,25:ncol(newcpm.jointmirnadge)]

#' Compute the correlations similarly to those above, but with the cpm counts and not the voom transformed values:
#' 
#' First, going column-wise (individuals) between the two kits:
corindcpm<-NULL

for (i in 1:ncol(newitscpm)){
	corindcpm[i]<-cor(newitscpm[,i], newbnfcpm[,i])
}
corindcpm
summary(corindcpm)

#' Create the vectors of correlation between individuals in the different kits:
#' 
#' 1034.ITS vs all of ITS:
cor1ITScpm<-cor(newitscpm[,"1034.ITS"], newitscpm[,2:ncol(newitscpm)])
cor1ITScpm
mean(cor1ITScpm)
range(cor1ITScpm)

#' 1034.ITS vs all of BNF:
cor1ITSBNFcpm<-cor(newitscpm[,"1034.ITS"], newbnfcpm[,2:ncol(newbnfcpm)])
cor1ITSBNFcpm
mean(cor1ITSBNFcpm)
range(cor1ITSBNFcpm)

#' 1034.BNF vs all of BNF:
cor1BNFcpm<-cor(newbnfcpm[,"1034.BNF"], newbnfcpm[,2:ncol(newbnfcpm)])
cor1BNFcpm
mean(cor1BNFcpm)
range(cor1BNFcpm)

#' 1034.BNF vs all of ITS:
cor1BNFITScpm<-cor(newbnfcpm[,"1034.BNF"], newitscpm[,2:ncol(newitscpm)])
cor1BNFITScpm
mean(cor1BNFITScpm)
range(cor1BNFITScpm)

#' Re-do the gene-wise correlation analysis with the cpm counts as well:
#' 
#' Create an empty matrix to fill with the correlation results:
cpmcor.results<-NULL

#' This for loop conducts the gene-wise correlation analysis for the 299 miRNAs between the two kits:

for (i in 1:nrow(newitscpm)){
    cpmcor.results[i]<-cor(newitscpm[i,], newbnfcpm[i,])
}

names(cpmcor.results)<-rownames(newitscpm)
summary(cpmcor.results)
hist(cpmcor.results)


#' ----------------------------------------------
#' 
#' Plot the average expression of the gene vs the -log10 pvalue of the correlation to see the relationship between expression and correlation:
#' 
#' Plotting just the raw rowMeans of the cpm jointmirnadge object:
plot(rowMeans(newcpm.jointmirnadge), -log(cor.results[,"p.value"]))
abline(0,1)

#' Plotting the log of the rowMeans of the cpm jointmirnadge object:
plot(log(rowMeans(newcpm.jointmirnadge)), -log(cor.results[,"p.value"]))
abline(0,1)


#' ## Save data
save(v, fit, toptab, cor.results, qcor.results, file="../1_voom_de_analysis_cor.results.Rdata")
save(jointmirnadge, design, file="../2_dge_object_design_mx.Rdata")