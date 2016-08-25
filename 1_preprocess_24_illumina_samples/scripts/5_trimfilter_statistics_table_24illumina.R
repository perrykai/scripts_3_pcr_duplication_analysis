#============================================
#  File:  5_trimfilter_statistics_table_24illumina.R
#  Directory Code:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/scripts
#  Date:  3/3/15 #UPDATED 8/15/16
#  Description:  Building a table of summary statistics from QC data obtained from cutadapt and fastq_quality_filter output.
#                This analysis is for the 24 illumina samples.
#                First, create a data.frame from the cutadapt adaptor trim step.
#                Columns: FileName | RawReads | TooShort | TooLong | Remaining | Proportion
#                All can be obtained from the Cutadapt Summary (mirsummary), and the remaining and proportion can be calculated in R.
#                Then, create a data.frame from the Fastq_qual_filter report (qualreport).
#                Columns: FileName | Input | Discarded | Output | PropInput | PropRaw
#                Where PropInput is the proportion of Output compared to Input reads
#                and PropRaw is the Proportion of Output compared to Raw Reads

#============================================
#  Input File Directory:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/1_cutadapt_adaptor_trim_size_filter_24illumina_output/mirsummary
#                         /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/2_fastx_quality_filter_24illumina_output/qualityreport
#
#  Input File(s):  xxxx_summary report for each input .fastq file
#                  xxxx_report for each .fastq file
#
#  Output File Directory: /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/5_trimfilter_statistics_table_24illumina_output
#
#  Output File(s): CutadaptStats.txt
#                  CutadaptStats2.txt (Contains final column of each library's proportion of total reads)
#                  FastxqfStats.txt
#                  FastxqfStats2.txt (Contains final column of each library's proportion of total reads)

#============================================
#  Module used: R/3.1.0
#============================================
#  Example Output (If Applicable):
#     CutadaptStats.txt:
#          FileName        RawReads     TooShort   TooLong         Remaining       Prop
# 1        1034_summary:   4879521         61546   2516159         2301816         0.471729909554647
# 2        1058_summary:   8058674         185654  3306176         4566844         0.566699186491475

#     CutadaptStats2.txt:
#          FileName        RawReads      ooShort   TooLong         Remaining       Prop                    PropTotal
# 1        1034_summary:   4879521         61546   2516159         2301816         0.471729909554647       0.0224718072944734
# 2        1058_summary:   8058674         185654  3306176         4566844         0.566699186491475       0.0371128578352225

#     FastxqfStats.txt:
#          FileName        Input       Discarded   Output          PropInput            RawReadsCA         PropRaw
# 1        1034_report:    2301816         26926   2274890         0.988302279591418       4879521         0.466211744964311
# 2        1058_report:    4566844         58605   4508239         0.987167286642592       8058674         0.559426898271353

#     FastxqfStats.txt:
#          FileName          Input     Discarded    Output         PropInput            RawReadsCA         PropRaw                 PropProcess
# 1        1034_report:    2301816         26926   2274890         0.988302279591418       4879521         0.466211744964311       0.0183163454903226
# 2        1058_report:    4566844         58605   4508239         0.987167286642592       8058674         0.559426898271353       0.0362982223654536

#============================================
#To extract the Raw Reads, use the processed reads from cutadapt report (mirsummary):
cd /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/1_cutadapt_adaptor_trim_size_filter_24illumina_output/mirsummary/
  grep "Processed reads" *> processed_reads.txt

#Extract the TooShort Reads:
  grep "Too short reads" * > tooshort_reads.txt

#Extract the TooLong Reads:
  grep "Too long reads" * > toolong_reads.txt

#All the above files saved in the directory /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/1_cutadapt_adaptor_trim_size_filter_24illumina_output/mirsummary

#====================================================================
#Create a data.frame of summary statistics for Cutadapt output:
#====================================================================
#Now, in R, reading the data in as data.frames to compile into one large data.frame.
R
rm(list=ls())
setwd("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/1_cutadapt_adaptor_trim_size_filter_24illumina_output/mirsummary")
input <- read.table("processed_reads.txt")
head(input)
#              V1        V2     V3      V4
# 1 1034_summary: Processed reads: 4879521
# 2 1058_summary: Processed reads: 8058674
# 3 1080_summary: Processed reads: 5692069
# 4 1096_summary: Processed reads: 7675196
# 5 1116_summary: Processed reads: 6483344
# 6 1134_summary: Processed reads: 8388301

tooshort <- read.table("tooshort_reads.txt")
head(tooshort)
#              V1  V2    V3     V4     V5     V6 V7        V8     V9
# 1 1034_summary: Too short reads:  61546  (1.3% of processed reads)
# 2 1058_summary: Too short reads: 185654  (2.3% of processed reads)
# 3 1080_summary: Too short reads: 505152  (8.9% of processed reads)
# 4 1096_summary: Too short reads: 299097  (3.9% of processed reads)
# 5 1116_summary: Too short reads: 644174  (9.9% of processed reads)
# 6 1134_summary: Too short reads: 923734 (11.0% of processed reads)

toolong <- read.table("toolong_reads.txt")
head(toolong)
#              V1  V2   V3     V4      V5     V6 V7        V8     V9
# 1 1034_summary: Too long reads: 2516159 (51.6% of processed reads)
# 2 1058_summary: Too long reads: 3306176 (41.0% of processed reads)
# 3 1080_summary: Too long reads: 2399589 (42.2% of processed reads)
# 4 1096_summary: Too long reads: 2852514 (37.2% of processed reads)
# 5 1116_summary: Too long reads: 1338337 (20.6% of processed reads)
# 6 1134_summary: Too long reads: 2079567 (24.8% of processed reads)

#Extract correct columns from data.frame:
input2<-input[,c("V1","V4")]
head(input2)
#              V1      V4
# 1 1034_summary: 4879521
# 2 1058_summary: 8058674
# 3 1080_summary: 5692069
# 4 1096_summary: 7675196
# 5 1116_summary: 6483344
# 6 1134_summary: 8388301

#Rename columns of data.frame:
colnames(input2)<- c("FileName","RawReads")
head(input2)
#        FileName RawReads
# 1 1034_summary:  4879521
# 2 1058_summary:  8058674
# 3 1080_summary:  5692069
# 4 1096_summary:  7675196
# 5 1116_summary:  6483344
# 6 1134_summary:  8388301

#Now, if I set up the other two data.frames similarly, I can merge them all together,
#Double checking that the order of the files remains the same throughout.
#Preparing TooShort data.frame:
tooshort2<-tooshort[,c("V1","V5")]
head(tooshort2)
#              V1     V5
# 1 1034_summary:  61546
# 2 1058_summary: 185654
# 3 1080_summary: 505152
# 4 1096_summary: 299097
# 5 1116_summary: 644174
# 6 1134_summary: 923734

colnames(tooshort2)<- c("FileName","TooShort")
head(tooshort2)
#        FileName TooShort
# 1 1034_summary:    61546
# 2 1058_summary:   185654
# 3 1080_summary:   505152
# 4 1096_summary:   299097
# 5 1116_summary:   644174
# 6 1134_summary:   923734

#Preparing TooLong data.frame:
toolong2<-toolong[,c("V1","V5")]
head(toolong2)
#              V1      V5
# 1 1034_summary: 2516159
# 2 1058_summary: 3306176
# 3 1080_summary: 2399589
# 4 1096_summary: 2852514
# 5 1116_summary: 1338337
# 6 1134_summary: 2079567

colnames(toolong2)<- c("FileName","TooLong")
head(toolong2)
#        FileName TooLong
# 1 1034_summary: 2516159
# 2 1058_summary: 3306176
# 3 1080_summary: 2399589
# 4 1096_summary: 2852514
# 5 1116_summary: 1338337
# 6 1134_summary: 2079567

#Now, merge the three data.frames, ensuring the file names stay in the correct order
merge1<-merge(input2,tooshort2,by="FileName")
head(merge1)
#        FileName RawReads TooShort
# 1 1034_summary:  4879521    61546
# 2 1058_summary:  8058674   185654
# 3 1080_summary:  5692069   505152
# 4 1096_summary:  7675196   299097
# 5 1116_summary:  6483344   644174
# 6 1134_summary:  8388301   923734

merge2<-merge(merge1,toolong2,by="FileName")
head(merge2)
#        FileName RawReads TooShort TooLong
# 1 1034_summary:  4879521    61546 2516159
# 2 1058_summary:  8058674   185654 3306176
# 3 1080_summary:  5692069   505152 2399589
# 4 1096_summary:  7675196   299097 2852514
# 5 1116_summary:  6483344   644174 1338337
# 6 1134_summary:  8388301   923734 2079567

#Then, calculate the Remaining reads (RawReads - TooShort - TooLong) column
merge2$Remaining<-(merge2$RawReads-merge2$TooShort-merge2$TooLong)
head(merge2)
#        FileName RawReads TooShort TooLong Remaining
# 1 1034_summary:  4879521    61546 2516159   2301816
# 2 1058_summary:  8058674   185654 3306176   4566844
# 3 1080_summary:  5692069   505152 2399589   2787328
# 4 1096_summary:  7675196   299097 2852514   4523585
# 5 1116_summary:  6483344   644174 1338337   4500833
# 6 1134_summary:  8388301   923734 2079567   5385000

#Now, calculate the proportion of RemainingReads compared to RawReads (Remaining/RawReads)
merge2$PropRemaining<-(merge2$Remaining/merge2$RawReads)
head(merge2)
#        FileName RawReads TooShort TooLong Remaining PropRemaining
# 1 1034_summary:  4879521    61546 2516159   2301816     0.4717299
# 2 1058_summary:  8058674   185654 3306176   4566844     0.5666992
# 3 1080_summary:  5692069   505152 2399589   2787328     0.4896863
# 4 1096_summary:  7675196   299097 2852514   4523585     0.5893771
# 5 1116_summary:  6483344   644174 1338337   4500833     0.6942147
# 6 1134_summary:  8388301   923734 2079567   5385000     0.6419655

#Range of proportion of reads kept after Adaptor Trimming:
range(merge2$PropRemaining)
#[1] 0.4505851 0.7057477

mean(merge2$PropRemaining)
# [1] 0.6132788

sd(merge2$PropRemaining)
# [1] 0.07058544

#==========================================================
#Creating a column of the proportion of each library vs total number of reads
#This column will be added to both CutadaptStats and FastxqfStats
#Fastxqfstats having a column of Proportion of Processed Reads vs Total Processed Reads
#==========================================================

#In the summary report, also want the proportion of raw reads in each library/total raw reads!
sum(merge2$RawReads)
# [1] 217139678

sum(merge2$Remaining)
# [1] 125975358

merge2$PropRawLibraryTotal<-(merge2$RawReads/(sum(merge2$RawReads)))
head(merge2)
#        FileName RawReads TooShort TooLong Remaining PropRemaining
# 1 1034_summary:  4879521    61546 2516159   2301816     0.4717299
# 2 1058_summary:  8058674   185654 3306176   4566844     0.5666992
# 3 1080_summary:  5692069   505152 2399589   2787328     0.4896863
# 4 1096_summary:  7675196   299097 2852514   4523585     0.5893771
# 5 1116_summary:  6483344   644174 1338337   4500833     0.6942147
# 6 1134_summary:  8388301   923734 2079567   5385000     0.6419655
#   PropRawLibraryTotal
# 1          0.02247181
# 2          0.03711286
# 3          0.02621386
# 4          0.03534681
# 5          0.02985794
# 6          0.03863090

merge2$PropRemainingLibraryTotal<-(merge2$Remaining/(sum(merge2$Remaining)))
head(merge2)
#        FileName RawReads TooShort TooLong Remaining PropRemaining
# 1 1034_summary:  4879521    61546 2516159   2301816     0.4717299
# 2 1058_summary:  8058674   185654 3306176   4566844     0.5666992
# 3 1080_summary:  5692069   505152 2399589   2787328     0.4896863
# 4 1096_summary:  7675196   299097 2852514   4523585     0.5893771
# 5 1116_summary:  6483344   644174 1338337   4500833     0.6942147
# 6 1134_summary:  8388301   923734 2079567   5385000     0.6419655
#   PropRawLibraryTotal PropRemainingLibraryTotal
# 1          0.02247181                0.01827195
# 2          0.03711286                0.03625188
# 3          0.02621386                0.02212598
# 4          0.03534681                0.03590849
# 5          0.02985794                0.03572788
# 6          0.03863090                0.04274646

#Save this file:
write.table(merge2, file = "/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/5_trimfilter_statistics_table_24illumina_output/1_cutadapt_adaptor_trim_stats_table_24illumina.txt", quote = FALSE, sep = "\t ", col.names = TRUE)

#This .txt file saved to /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/5_trimfilter_statistics_table_24illumina_output/

q()
n

#====================================================================
#Create a summary statistics table for Fastq_qual_Filter Data:
#====================================================================
cd /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/2_fastx_quality_filter_24illumina_output/qualityreport/

#Now, to follow the same procedure with the Fastq_qual_filter report to obtain the statistics there:
grep Output * >output_reads_fastx_quality_filter.txt
grep Input * >input_reads_fastx_quality_filter.txt
grep discarded * >discarded_reads_fastx_quality_filter.txt

#Now, build the table in R:
#Input/Read files:
R
rm(list=ls())
setwd("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/2_fastx_quality_filter_24illumina_output/qualityreport/")

inputqf<-read.table("input_reads_fastx_quality_filter.txt")
head(inputqf)
#                   V1      V2     V3
# 1 1034_report:Input: 2301816 reads.
# 2 1058_report:Input: 4566844 reads.
# 3 1080_report:Input: 2787328 reads.
# 4 1096_report:Input: 4523585 reads.
# 5 1116_report:Input: 4500833 reads.
# 6 1134_report:Input: 5385000 reads.

discardqf<-read.table("discarded_reads_fastx_quality_filter.txt")
head(discardqf)
#                      V1     V2   V3          V4     V5
# 1 1034_report:discarded  26926 (1%) low-quality reads.
# 2 1058_report:discarded  58605 (1%) low-quality reads.
# 3 1080_report:discarded  40900 (1%) low-quality reads.
# 4 1096_report:discarded  61873 (1%) low-quality reads.
# 5 1116_report:discarded  57917 (1%) low-quality reads.
# 6 1134_report:discarded 105956 (1%) low-quality reads.

outputqf<-read.table("output_reads_fastx_quality_filter.txt")
head(outputqf)
#                    V1      V2     V3
# 1 1034_report:Output: 2274890 reads.
# 2 1058_report:Output: 4508239 reads.
# 3 1080_report:Output: 2746428 reads.
# 4 1096_report:Output: 4461712 reads.
# 5 1116_report:Output: 4442916 reads.
# 6 1134_report:Output: 5279044 reads.

#Trim needed columns and names of files for merging:

#For this analysis, I only need $V1 and $V2 from each file. The rest can be discarded.
#Input file:
inputqf2<-inputqf[,c("V1","V2")]
head(inputqf2)
#                   V1      V2
# 1 1034_report:Input: 2301816
# 2 1058_report:Input: 4566844
# 3 1080_report:Input: 2787328
# 4 1096_report:Input: 4523585
# 5 1116_report:Input: 4500833
# 6 1134_report:Input: 5385000

colnames(inputqf2)<- c("FileName","Input")
head(inputqf2)
#             FileName   Input
# 1 1034_report:Input: 2301816
# 2 1058_report:Input: 4566844
# 3 1080_report:Input: 2787328
# 4 1096_report:Input: 4523585
# 5 1116_report:Input: 4500833
# 6 1134_report:Input: 5385000

#Discard File:
discardqf2<-discardqf[,c("V1","V2")]
head(discardqf2)
#                      V1     V2
# 1 1034_report:discarded  26926
# 2 1058_report:discarded  58605
# 3 1080_report:discarded  40900
# 4 1096_report:discarded  61873
# 5 1116_report:discarded  57917
# 6 1134_report:discarded 105956

colnames(discardqf2)<- c("FileName","Discarded")
head(discardqf2)
#                FileName Discarded
# 1 1034_report:discarded     26926
# 2 1058_report:discarded     58605
# 3 1080_report:discarded     40900
# 4 1096_report:discarded     61873
# 5 1116_report:discarded     57917
# 6 1134_report:discarded    105956

#Output File:
outputqf2<-outputqf[,c("V1","V2")]
head(outputqf2)
#                    V1      V2
# 1 1034_report:Output: 2274890
# 2 1058_report:Output: 4508239
# 3 1080_report:Output: 2746428
# 4 1096_report:Output: 4461712
# 5 1116_report:Output: 4442916
# 6 1134_report:Output: 5279044

colnames(outputqf2)<- c("FileName","Output")
head(outputqf2)
#              FileName  Output
# 1 1034_report:Output: 2274890
# 2 1058_report:Output: 4508239
# 3 1080_report:Output: 2746428
# 4 1096_report:Output: 4461712
# 5 1116_report:Output: 4442916
# 6 1134_report:Output: 5279044

#Now I need to trim the text in the file name after "_report:" to merge the files
#Example using outputqf2 file:
outputqf3<-as.data.frame(sapply(outputqf2,gsub,pattern="Output:",replacement=""))
head(outputqf3)
#       FileName  Output
# 1 1034_report: 2274890
# 2 1058_report: 4508239
# 3 1080_report: 2746428
# 4 1096_report: 4461712
# 5 1116_report: 4442916
# 6 1134_report: 5279044

#Repeat this step with the other two files:
inputqf3<-as.data.frame(sapply(inputqf2,gsub,pattern="Input:",replacement=""))
head(inputqf3)
#       FileName   Input
# 1 1034_report: 2301816
# 2 1058_report: 4566844
# 3 1080_report: 2787328
# 4 1096_report: 4523585
# 5 1116_report: 4500833
# 6 1134_report: 5385000

discardqf3<-as.data.frame(sapply(discardqf2,gsub,pattern="discarded",replacement=""))
head(discardqf3)
#       FileName Discarded
# 1 1034_report:     26926
# 2 1058_report:     58605
# 3 1080_report:     40900
# 4 1096_report:     61873
# 5 1116_report:     57917
# 6 1134_report:    105956

#Merge the files to create one data.frame:
#Now that the FileNames match, I can merge the three files into one data.frame
mergeqf1<-merge(inputqf3,discardqf3,by="FileName")
head(mergeqf1)
#       FileName   Input Discarded
# 1 1034_report: 2301816     26926
# 2 1058_report: 4566844     58605
# 3 1080_report: 2787328     40900
# 4 1096_report: 4523585     61873
# 5 1116_report: 4500833     57917
# 6 1134_report: 5385000    105956

mergeqf2<-merge(mergeqf1,outputqf3,by="FileName")
head(mergeqf2)
#       FileName   Input Discarded  Output
# 1 1034_report: 2301816     26926 2274890
# 2 1058_report: 4566844     58605 4508239
# 3 1080_report: 2787328     40900 2746428
# 4 1096_report: 4523585     61873 4461712
# 5 1116_report: 4500833     57917 4442916
# 6 1134_report: 5385000    105956 5279044

#Add additional columns of information (Proportions):
#Next, add a column of PropInput, the Proportion of the Output compared to the Proportion of Reads Input
#Need to convert columns 2-4 to numeric, by first converting to character
mergeqf2[,2]<-as.numeric(as.character(mergeqf2[,2]))
mergeqf2[,3]<-as.numeric(as.character(mergeqf2[,3]))
mergeqf2[,4]<-as.numeric(as.character(mergeqf2[,4]))
head(mergeqf2)
#       FileName   Input Discarded  Output
# 1 1034_report: 2301816     26926 2274890
# 2 1058_report: 4566844     58605 4508239
# 3 1080_report: 2787328     40900 2746428
# 4 1096_report: 4523585     61873 4461712
# 5 1116_report: 4500833     57917 4442916
# 6 1134_report: 5385000    105956 5279044

mergeqf2$PropReadsRemaining<-(mergeqf2$Output/mergeqf2$Input)
head(mergeqf2)
#       FileName   Input Discarded  Output PropReadsRemaining
# 1 1034_report: 2301816     26926 2274890          0.9883023
# 2 1058_report: 4566844     58605 4508239          0.9871673
# 3 1080_report: 2787328     40900 2746428          0.9853264
# 4 1096_report: 4523585     61873 4461712          0.9863221
# 5 1116_report: 4500833     57917 4442916          0.9871319
# 6 1134_report: 5385000    105956 5279044          0.9803239

mean(mergeqf2$Output)
# [1] 5175000
range(mergeqf2$Output)
# [1]  1402261 20300167

#I believe I can then calculate the Proportion of Raw Reads
#By dividing Output here by "RawReads" column from CutadaptStats.txt
cutadaptstats<-read.table("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/5_trimfilter_statistics_table_24illumina_output/1_cutadapt_adaptor_trim_stats_table_24illumina.txt")
head(cutadaptstats)
#        FileName RawReads TooShort TooLong Remaining PropRemaining
# 1 1034_summary:  4879521    61546 2516159   2301816     0.4717299
# 2 1058_summary:  8058674   185654 3306176   4566844     0.5666992
# 3 1080_summary:  5692069   505152 2399589   2787328     0.4896863
# 4 1096_summary:  7675196   299097 2852514   4523585     0.5893771
# 5 1116_summary:  6483344   644174 1338337   4500833     0.6942147
# 6 1134_summary:  8388301   923734 2079567   5385000     0.6419655
#   PropRawLibraryTotal PropRemainingLibraryTotal
# 1          0.02247181                0.01827195
# 2          0.03711286                0.03625188
# 3          0.02621386                0.02212598
# 4          0.03534681                0.03590849
# 5          0.02985794                0.03572788
# 6          0.03863090                0.04274646

#I could include a column of RawReadsCutadapt after the PropReadsRemaining, merging the files to ensure the files match correctly.
mergeqf2$RawReadsCutadapt<-cutadaptstats$RawReads
head(mergeqf2)
#       FileName   Input Discarded  Output PropReadsRemaining RawReadsCutadapt
# 1 1034_report: 2301816     26926 2274890          0.9883023          4879521
# 2 1058_report: 4566844     58605 4508239          0.9871673          8058674
# 3 1080_report: 2787328     40900 2746428          0.9853264          5692069
# 4 1096_report: 4523585     61873 4461712          0.9863221          7675196
# 5 1116_report: 4500833     57917 4442916          0.9871319          6483344
# 6 1134_report: 5385000    105956 5279044          0.9803239          8388301

#Calculate the proportion of reads remaining versus the raw reads input (before adaptor trimming)
mergeqf2$PropRawReadsCutadapt<-(mergeqf2$Output/mergeqf2$RawReadsCutadapt)
head(mergeqf2)
#       FileName   Input Discarded  Output PropReadsRemaining RawReadsCutadapt
# 1 1034_report: 2301816     26926 2274890          0.9883023          4879521
# 2 1058_report: 4566844     58605 4508239          0.9871673          8058674
# 3 1080_report: 2787328     40900 2746428          0.9853264          5692069
# 4 1096_report: 4523585     61873 4461712          0.9863221          7675196
# 5 1116_report: 4500833     57917 4442916          0.9871319          6483344
# 6 1134_report: 5385000    105956 5279044          0.9803239          8388301
#   PropRawReadsCutadapt
# 1            0.4662117
# 2            0.5594269
# 3            0.4825008
# 4            0.5813157
# 5            0.6852815
# 6            0.6293341

#Can also check that "Remaining" from CutadaptStats.txt matches "Input" in this data.frame
(mergeqf2$Input-cutadaptstats$Remaining)
 # [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

#Range of proportion of reads kept after both processing steps:
range(mergeqf2$PropRawReadsCutadapt)
# [1] 0.4438810 0.6952216

mean(mergeqf2$PropRawReadsCutadapt)
# [1] 0.6047061

sd(mergeqf2$PropRawReadsCutadapt)
# [1] 0.06955766

#Proportion of Processed Reads per Library vs Total Processed Reads:
#Sum of reads output from fastq_quality_filter step
sum(mergeqf2$Output)
# [1] 124199994

mergeqf2$PropProcessedLibraryTotal<-(mergeqf2$Output/(sum(mergeqf2$Output)))
head(mergeqf2)
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

#Save this file:
getwd()
# [1] "/mnt/research/pigeqtl/analyses/microRNA/1_preprocess_fastq_files/1_preprocess_fastq_files_output/2_fastx_quality_filter_output/qualityreport" 
write.table(mergeqf2, file = "/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/5_trimfilter_statistics_table_24illumina_output/2_fastx_quality_filter_stats_table_24illumina.txt", quote = FALSE, sep = "\t ", col.names = TRUE)

range(mergeqf2$PropProcessedLibraryTotal)
# [1] 0.01129035 0.16344741
mean(mergeqf2$PropProcessedLibraryTotal)
# [1] 0.04166667
sd(mergeqf2$PropProcessedLibraryTotal)
# [1] 0.03130403
q()
n