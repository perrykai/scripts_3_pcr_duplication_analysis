#============================================
#  File:  3_quality_summary_plot_24illumina.R
#  Directory Code:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/scripts
#  Date:  2/26/15 #UPDATED 8/12/16
#  Description:  This code creates one graph plotting the Q1, median, and Q3 Phred scores
#                for each position of each read, for the preliminary 24 samples.
#============================================
#  Input File Directory:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/3_trimfilter_statistics_24illumina_output/2_trimmed_filtered_fastq_quality_assessment

#  Input File(s):  xxxx_posttrimqstats.txt (One output file from the preprocessing steps for each sample)

#  Output File Directory: /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/3_trimfilter_statistics_24illumina_output/3_quality_summary_24illumina_plot_output

#  Output File(s): 1_posttrimfilter_QCsummary_plot_24illumina.pdf
#============================================
#  Module Used: R/3.1.0 (run interactively)
#============================================
#  Example Output (If Applicable):
#============================================
rm(list=ls())
setwd("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/3_trimfilter_statistics_24illumina_output/2_trimmed_filtered_fastq_quality_assessment/")

files<-list.files(pattern="*.txt")
id<-unlist(lapply(strsplit(files,"_"), function(x) x[1]))
lane<-unlist(lapply(strsplit(files,"_"), function(x) x[3]))
pe<-unlist(lapply(strsplit(files,"_"), function(x) x[4]))
nm<-paste(id,lane,pe,sep="_")

tb<-read.table(files[1],header=T)
pdf("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/3_quality_summary_24illumina_plot_output/1_posttrimfilter_QCsummary_plot_24illumina.pdf")
plot(tb$column,tb$med,ylim=c(0,45),type="n",xlab="read position",ylab="Phred score",main="Quality Scores Preliminary 24 Samples")
for (i in files){
  tb<-read.table(i,header=T)
  points(tb$column,tb$med,type="l",lwd=2,lty="dotted")
  points(tb$column,tb$Q1,type="l",lty="dashed",col="red")
  points(tb$column,tb$Q3,type="l",lty="dashed",col="blue")
  legend("bottomright", c("Median","Q1","Q3"), col=c("black","red","blue"),lty=c("dotted","dashed","dashed"),lwd=c(2,1,1))
}
dev.off()
q()
n