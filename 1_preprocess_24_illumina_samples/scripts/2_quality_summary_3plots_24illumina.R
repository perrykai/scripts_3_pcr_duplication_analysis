#============================================
#  File:  2_quality_summary_3plots_24illumina.R
#  Directory Code:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/scripts
#  Date:  2/18/15 #UPDATED 8/12/16
#  Description:  This code creates three graphs for each sample, displaying quality control information including:
#                One graph of the Phred score at each read position, one graph of the average base composition at each read position, 
#                and one graph of the number of bases called at each read position. 
#============================================
#  Input File Directory:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/3_trimfilter_statistics_24illumina_output/2_trimmed_filtered_fastq_quality_assessment

#  Input File(s):  xxxx_posttrimqstats.txt (One output file from the preprocessing steps for each sample)

#  Output File Directory: /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/2_quality_summary_24illumina_3plots_output

#  Output File(s): 1_posttrimfilter_QC3plots_24illumina.pdf

#============================================
#  Example Output (If Applicable):
#  This script was run interactively with R/3.1.0
#============================================

setwd ("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/3_trimfilter_statistics_24illumina_output/2_trimmed_filtered_fastq_quality_assessment/")

pdf("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/2_quality_summary_24illumina_3plots_output/1_posttrimfilter_QC3plots_24illumina.pdf")
files<-list.files(pattern="*.txt")
for (i in files){
  tb<-read.table(i,header=T)
  totcount<-rowSums(tb[,c("A_Count","C_Count","G_Count","T_Count")])
  prop<-sweep(tb[,c("A_Count","C_Count","G_Count","T_Count")],1,"/",STATS=totcount)
  nf <- layout(matrix(c(1,2,3,1,2,3,1,2,3), 3, 3, byrow=TRUE), respect=TRUE)
  plot(tb$column,tb$med,ylim=c(0,45),type="l",lwd=2,xlab="read position",ylab="Phred score",main=i)
  points(tb$column,tb$Q1,type="l",lty="dashed")
  points(tb$column,tb$Q3,type="l",lty="dashed")
  
  plot(tb$column,prop[,1],ylim=c(0,max(prop)),type="l",xlab="read position",
       ylab="average base composition",,main=i)
  for (j in 2:4){
    points(tb$column,prop[,j],type="l")
  }
  
  plot(tb$column,totcount,ylim=c(0,max(totcount)),type="l",xlab="read position",
       ylab="bases called",main=i)
}
dev.off()
q()
n