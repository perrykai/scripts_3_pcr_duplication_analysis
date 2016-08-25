#============================================
#  File:  4a_readlengthdist_plot_table_24illumina.R
#  Directory Code:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/scripts
#  Date:  2/19/15 #UPDATED 8/12/16
#  Description:  This code creates a .pdf plot of the read length distributions of the 24 preliminary smallRNA libraries.
#                Additionally, it creates a .txt file of a table containing columns "File", "ReadLength", "NumReads", and "PropTotal",
#                describing the most common read length in each library, the number of reads of that length, and the proportion of total 
#                reads that are that length.
#============================================
#  Input File Directory:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/4_readlength_distribution_24illumina_output

#  Input File(s):  xxxx_readlengthdistfull.txt 

#  Output File Directory: /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/4_readlengthdist_histogram_24illumina_plot_output/

#  Output File(s): 1_readlengthdist_fullplot_24illumina.pdf (visual plot of read length distribution)
#                  2_readlength_frequency_24illumina.txt  (Table containing read length distribution info)

#============================================
#  Example Output (If Applicable):

# "File" "ReadLength" "NumReads" "PropTotal"
#  "1" "1034_readlengthdistfull.txt" 22 1256308 0.552249998901046
#  "2" "1058_readlengthdistfull.txt" 22 2570618 0.570204463427959
#  "3" "1080_readlengthdistfull.txt" 22 1098600 0.400010486348086
#  "4" "1096_readlengthdistfull.txt" 22 2297557 0.514949642648382
#  "5" "1116_readlengthdistfull.txt" 22 2257224 0.50805011843573
#  "6" "1134_readlengthdistfull.txt" 22 2714353 0.514175104431787

#============================================
#Run visualization R code: MODIFIED 2/19/2015 with help of Dr. Steibel

rm(list=ls())
setwd("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/4_readlength_distribution_24illumina_output/")

pdf("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/4_readlengthdist_histogram_24illumina_plot_output/1_readlengthdist_fullplot_24illumina.pdf")

files<-list.files(pattern="*_readlengthdistfull.txt")
mx<-matrix(0,length(files),3)
reads<-read.csv(file=files[1], sep="", header=FALSE)
plot (reads$V2,reads$V1,type="l", main="Read Length Dist. 24 Libraries", xlab="Read Length",ylab="Occurences",ylim=c(0,9000000))

#The mx matrix is to view the read length with the highest frequency of reads for all 24 samples.
#First, extract the read length that has the highest frequency into the first column:
mx[1,1]<-reads$V2[which.max(reads$V1)]
#The second column is the number of reads in the specified read length:
mx[1,2]<-max(reads$V1)
#The third column is the maximum reads divided by the total number of reads in the sample.
#This determines the proportion of reads that are the specified read length:
mx[1,3]<-max(reads$V1)/sum(reads$V1)
j<-1
for (i in files[2:length(files)]) {
  reads2<-read.csv(file = i, sep = "", header = FALSE)
  points(reads2$V2,reads2$V1,type="l")
  j<-j+1
  mx[j,1]<-reads2$V2[which.max(reads2$V1)]
  mx[j,2]<-max(reads2$V1)
  mx[j,3]<-max(reads2$V1)/sum(reads2$V1)
}
mx<-data.frame(files,mx)
colnames(mx)<- c("File","ReadLength","NumReads","PropTotal")
dev.off()
print (mx)
write.table(mx, file="/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/4_readlengthdist_histogram_24illumina_plot_output/2_readlength_frequency_24illumina.txt", sep=" ", col.names=TRUE)
q()
n