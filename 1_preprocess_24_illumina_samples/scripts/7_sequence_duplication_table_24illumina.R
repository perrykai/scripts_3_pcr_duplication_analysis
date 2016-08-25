#============================================
#  File:  SeqDupInfoUpdate_prelim24.R
#  Directory Code:  /mnt/research/pigeqtl/analyses/microRNA/Code
#  Date:  8/4/15 #UPDATED 8/15/16
#  Description:  This code creates an updated matrix containing PCR duplication information
#                The columns contain: FileName: The name of each sample (Pig ID)
#                UniqFrag: The number of unique fragments in the collapsed fasta file, before random adaptor(R.A.) trimming
#                PropExp>x: The proportion of unique fragments expressed greater than 1, 10, 100, or 1000 times
#                MeanExpPre: The average expression of each unique fragment in the collapsed fasta file before R.A. trimming
#============================================
#  Input File Directory:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/5_collapsed_fasta_24illumina_output/

#  Input File(s):  xxxx_collapsed.fa (One collapsed fasta file from the preprocessing steps for each sample)

#  Output File Directory: /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/7_sequence_duplication_table_24illumina_output

#  Output File(s): 1_sequence_duplication_information_24illumina.txt
#============================================
#  Module Used: R/3.1.0
#============================================
#  Example Output (If Applicable):

#      FileName    UniqFrag        PropExp>1               PropExp>10              PropExp>100             PropExp>1000            MeanExpPre             UniqSeq  MeanExpPost
# 1        1036    1210210         0.240726816007139       0.0243775873608713      0.001274985333124       1.07419373497162e-05    2.29977937713288        102037  11.8605015827592
# 2        1041    1393661         0.233143497593748       0.0263571987735898      0.00167903098386193     1.07630191273201e-05    2.42857481123457        134835  10.3360477620796
# 3        1049    1594434         0.262079835226795       0.028314749936341       0.00206217378706174     3.70037267143074e-05    2.66241437400357        162399  9.81800380544215

#============================================
R
setwd("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/5_collapsed_fasta_24illumina_output")
rm(list=ls())

library("ShortRead")

files<-list.files(pattern="_collapsed.fa")
seqdupinfo<-matrix(NA, nrow=length(files), ncol=7)

  for (i in 1:length(files)){
    files1<-files[[i]]
    reads <- readFasta(files1) #Read in the file as a fasta file using ShortRead package
    seqs <- sread(reads) #Isolate the sequences from the collapsed fasta file
    uniqfrag<-length(seqs) #The length of this file is the number of unique fragments, prior to random adaptor trimming
    iddf<-as.data.frame(id(reads)) #Create a data.frame of the sequence ids (this value contains the read counts of each fragment)
    iddf$expression <- as.numeric(sub("^.*-", "", iddf$x)) #isolate the read count data
    expmx<-as.matrix(iddf$expression)
    expg1<-apply(expmx, 2, function(x) (sum(x > 1))/length(expmx))      #Calculates proportion of unique fragments expressed more than once in the sample
    expg10<-apply(expmx, 2, function(x) (sum(x > 10))/length(expmx))    #Calculates proportion of unique fragments expressed more than 10x in the sample
    expg100<-apply(expmx, 2, function(x) (sum(x > 100))/length(expmx))  #Calculates proportion of unique fragments expressed more than 100x in the sample
    expg1000<-apply(expmx, 2, function(x) (sum(x > 1000))/length(expmx))#Calculates proportion of unique fragments expressed more than 1000x in the sample
    mnuniqfrag<-mean(iddf$expression) #Take the average read counts of all unique fragments
    filenmcp<-sub("_collapsed.fa", "", files1) #Create the file name (take the name of the file, remove the _collapsed.fa to obtain the pig ID)
    
    seqdupinfo[i,1]<-as.numeric(filenmcp) #assign columns of matrix
    seqdupinfo[i,2]<-as.numeric(uniqfrag)
    seqdupinfo[i,3]<-as.numeric(expg1)
    seqdupinfo[i,4]<-as.numeric(expg10)
    seqdupinfo[i,5]<-as.numeric(expg100)
    seqdupinfo[i,6]<-as.numeric(expg1000)
    seqdupinfo[i,7]<-as.numeric(mnuniqfrag)
    colnames(seqdupinfo)<-c("FileName","UniqFrag", "PropExp>1", "PropExp>10", "PropExp>100", "PropExp>1000", "MeanExpPre")
  }

write.table(seqdupinfo, file="/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/7_sequence_duplication_table_24illumina_output/1_sequence_duplication_information_24illumina.txt", quote = FALSE, sep = "\t ", col.names = TRUE)
q()
n