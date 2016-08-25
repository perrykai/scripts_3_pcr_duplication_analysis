#============================================
#  File:  6_librarysize_scatterplot_24illumina.R
#  Directory Code:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/scripts
#  Date:  7/14/15 #UPDATED 8/15/16
#  Description:  This code creates a scatterplot plotting the total processed library size of each pig
#                Pig ID is on the X axis, and Library Size (Reads) is on the Y axis.
#============================================
#  Input File Directory:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/5_trimfilter_statistics_table_24illumina_output

#  Input File(s):  2_fastx_quality_filter_stats_table_24illumina.txt

#  Output File Directory: /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/6_librarysize_scatterplot_24illumina_output

#  Output File(s): 1_librarysize_scatterplot_24illumina.pdf

#============================================
#  Module Used: R/3.1.0
#============================================
#  Example Output (If Applicable):
#============================================
R
setwd("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/5_trimfilter_statistics_table_24illumina_output")
getwd()
# [1] "/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/5_trimfilter_statistics_table_24illumina_output"
rm(list=ls())
ls()
# character(0)
 
#==============================================================================
#Now, pull in the preliminary 24 samples 
qcstats24<-read.csv("2_fastx_quality_filter_stats_table_24illumina.txt", sep="", header=TRUE)
head(qcstats24)
#       FileName   Input Discarded  Output PropReadsRemaining RawReadsCutadapt
# 1 1034_report: 2301816     26926 2274890          0.9883023          4879521
# 2 1058_report: 4566844     58605 4508239          0.9871673          8058674
# 3 1080_report: 2787328     40900 2746428          0.9853264          5692069
# 4 1096_report: 4523585     61873 4461712          0.9863221          7675196
# 5 1116_report: 4500833     57917 4442916          0.9871319          6483344
# 6 1134_report: 5385000    105956 5279044          0.9803239          8388301
#   PropRawReadsCutadapt PropProcessedLibraryTotal
# 1            0.4662117                0.01831635
# 2            0.5594269                0.03629822
# 3            0.4825008                0.02211295
# 4            0.5813157                0.03592361
# 5            0.6852815                0.03577227
# 6            0.6293341                0.04250438

qcstats24$Filenamenum24<-as.numeric(sub("_report:", "", qcstats24$FileName))
head(qcstats24)
#       FileName   Input Discarded  Output PropReadsRemaining RawReadsCutadapt
# 1 1034_report: 2301816     26926 2274890          0.9883023          4879521
# 2 1058_report: 4566844     58605 4508239          0.9871673          8058674
# 3 1080_report: 2787328     40900 2746428          0.9853264          5692069
# 4 1096_report: 4523585     61873 4461712          0.9863221          7675196
# 5 1116_report: 4500833     57917 4442916          0.9871319          6483344
# 6 1134_report: 5385000    105956 5279044          0.9803239          8388301
#   PropRawReadsCutadapt PropProcessedLibraryTotal Filenamenum24
# 1            0.4662117                0.01831635          1034
# 2            0.5594269                0.03629822          1058
# 3            0.4825008                0.02211295          1080
# 4            0.5813157                0.03592361          1096
# 5            0.6852815                0.03577227          1116
# 6            0.6293341                0.04250438          1134

sum(qcstats24$Output)
# [1] 124199994

#Plot the results to see what they look like
setwd("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/6_librarysize_scatterplot_24illumina_output/")
pdf("1_librarysize_scatterplot_24illumina.pdf")
plot(qcstats24$Filenamenum24, qcstats24$Output, xlim=c(1034,1662), ylim=c(0,22000000), 
     xlab="Pig ID", ylab ="Total Processed Reads", 
     main="Total Processed Reads per Pig",type="p", pch=19, col="blue",
     cex.main=1.8,cex.lab=1.5,cex.axis=1.2)
dev.off()
q()
n