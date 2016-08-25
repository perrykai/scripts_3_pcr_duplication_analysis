#============================================
#  File:  1_preprocess_24illumina.sh
#  Directory Code:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/scripts
#  Date:  2/11/15 #UPDATED 8/12/16
#  Description:  This code completes all pre-processing steps for the preliminary 24 samples of smallRNAseq data. The steps are as follows:
#                1. Load modules cutadapt and FASTX
#                2. Copy the raw data from the datasets directory to the 3_pcr_duplication_analysis/1_preprocess_24_illumina_samples directory and unzip the files
#                3. Move to the microRNA directory to begin the analysis
#                4. Use fastx_quality_stats to assess the quality of the preliminary 24 samples prior to analysis
#                5. Use cutadapt to trim the universal adaptor sequence and the length of the reads to between 18-30nt
#                   (This step creates multiple output files useful for further analysis; see output files below)
#                6. Use fastq_quality_filter to filter out the low-quality reads (reads with >50% of reads having Phred <30 are eliminated)
#                7. Use fastx_quality_stats to assess the quality of the preliminary 24 samples following analysis
#                8. Use awk to extract the second line from each group of four lines (the sequence line in a .fastq file),
#                   then measure the length of the sequence corresponding to that line (read length), incrementing the array cell
#                   corresponding to that length. When all the lines have been read, it loops over the array to print its content.
#                   (Awk code obtained from Biostars: https://www.biostars.org/p/72433/)
#                9. Use fastx_collapser to collapse the adaptor trimmed and quality filtered .fastq files into .fa files,
#                   while combining duplicated reads into a unique sequence ID with an expression tag.
#                   It also creates a verbose output file containing information on number of input and output sequences
#                   and the number of reads represented. This code is used in the first step of investigating
#                   PCR Duplication in the smallRNA seq data sequenced 1/2015.
#               10. After submitting the job, go back to the "3_pcr_duplication_analysis/1_preprocess_24_illumina_samples" directory and remove the unzipped raw data file. 
#============================================
#  Input File Directory:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/

#  Input File(s):  One raw .fastq file for each sample.

#  Output File Directory: /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/1_cutadapt_adaptor_trim_size_filter_24illumina_output/mirinfo
#                         /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/1_cutadapt_adaptor_trim_size_filter_24illumina_output/mirrest
#                         /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/1_cutadapt_adaptor_trim_size_filter_24illumina_output/mirtoolong
#                         /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/1_cutadapt_adaptor_trim_size_filter_24illumina_output/mirsummary
#                         /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/1_cutadapt_adaptor_trim_size_filter_24illumina_output
#                         /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/2_fastx_quality_filter_24illumina_output
#                         /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/2_fastx_quality_filter_24illumina_output/qualityreport
#                         /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/3_trimfilter_statistics_24illumina_output/1_raw_fastq_quality_assessment
#                         /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/3_trimfilter_statistics_24illumina_output/2_trimmed_filtered_fastq_quality_assessment
#                         /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/4_readlength_distribution_24illumina_output
#                         /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/5_collapsed_fasta_24illumina_output
#                         /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/5_collapsed_fasta_24illumina_output/verbose

#  Output File(s): xxxx_pretrimqstats.txt    (output of pre-trim fastx_quality_stats)
#                  xxxx_info                 (output from cutadapt, info file containing information about each read and its adaptor matches)
#                  xxxx_rest                 (output from cutadapt, when the adaptor matches in the middle of a read, the rest of the read is written here)
#                  xxxx_toolong              (output from cutadapt, when the read is too long (>30nt), it is written to this file)
#                  xxxx_cutadaptout.fastq    (output from cutadapt, output .fastq file of the remaining adaptor-trimmed reads meeting length requirements)
#                  xxxx_summary              (output from cutadapt, verbose output for each sample saved to individual files)
#                  xxxx_fastqfilterout.fastq (output from fastq_quality_filter)
#                  xxxx_report               (output from fastq_quality_filter, verbose output file for each sample)
#                  xxxx_posttrimqstats.txt   (output of post-trim and quality filter fastx_quality_stats)
#                  xxxx_collapsed.fa         (One output file of collapsed reads for each sample)
#                  xxxx_verbose.txt          (A verbose output file for each collapsed sample)

#============================================
#Final Code for Preprocessing All 24 Samples:

f1=(`ls /mnt/research/pigeqtl/datasets/microrna/ -1|grep gz`)
f2=(`ls /mnt/research/pigeqtl/datasets/microrna/ -1|grep fastq|cut -f1 -d_`)
echo $f1
echo $f2
cd /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/scripts
  
for (( i = 0 ; i < ${#f1[@]} ; i++ )) do
 echo '#!/bin/sh  -login' > base.sh
 echo '#PBS -l nodes=1:ppn=1:intel14,walltime=2:00:00,mem=2Gb'  >> base.sh
 echo '#PBS -N illumina_sample_'${f2[$i]}  >> base.sh
 echo '#PBS -j oe'  >> base.sh
 echo '#PBS -o /mnt/research/pigeqtl/analyses/microRNA/OutputsErrors/3_preprocess_24illumina_fastq/'  >> base.sh
 echo '#PBS -m a'  >> base.sh
 echo '#PBS -M perrykai@msu.edu'  >> base.sh
 
 echo 'module load cutadapt/1.4.1'  >> base.sh
 echo 'module load FASTX/0.0.14'  >> base.sh

 echo 'cd /mnt/research/pigeqtl/datasets/microrna/'  >> base.sh
 echo 'cp' ${f1[$i]} '/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/'  >> base.sh
 echo 'cd /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/'  >> base.sh
 echo 'gunzip' ${f1[$i]} >> base.sh   
 
 echo 'cd /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output'  >> base.sh   

 echo 'fastx_quality_stats -i /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/'`echo ${f1[$i]}|cut -d. -f1,2` '-o /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/3_trimfilter_statistics_24illumina_output/1_raw_fastq_quality_assessment/'${f2[$i]}'_pretrimqstats.txt -Q33'  >> base.sh

 echo 'cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC -m 18 -M 30 --info-file=/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/1_cutadapt_adaptor_trim_size_filter_24illumina_output/mirinfo/'${f2[$i]}'_info --rest-file=/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/1_cutadapt_adaptor_trim_size_filter_24illumina_output/mirrest/'${f2[$i]}'_rest --too-long-output=/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/1_cutadapt_adaptor_trim_size_filter_24illumina_output/mirtoolong/'${f2[$i]}'_toolong /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/'`echo ${f1[$i]}|cut -d. -f1,2` '-o /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/1_cutadapt_adaptor_trim_size_filter_24illumina_output/'${f2[$i]}'_cutadaptout.fastq > /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/1_cutadapt_adaptor_trim_size_filter_24illumina_output/mirsummary/'${f2[$i]}'_summary'  >> base.sh

 echo 'fastq_quality_filter -q 30 -p 50 -i /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/1_cutadapt_adaptor_trim_size_filter_24illumina_output/'${f2[$i]}'_cutadaptout.fastq -o /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/2_fastx_quality_filter_24illumina_output/'${f2[$i]}'_fastqfilterout.fastq -v > /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/2_fastx_quality_filter_24illumina_output/qualityreport/'${f2[$i]}'_report -Q33'  >> base.sh

 echo 'fastx_quality_stats -i /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/2_fastx_quality_filter_24illumina_output/'${f2[$i]}'_fastqfilterout.fastq -o /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/3_trimfilter_statistics_24illumina_output/2_trimmed_filtered_fastq_quality_assessment/'${f2[$i]}'_posttrimqstats.txt -Q33'  >> base.sh

 echo 'cat /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/2_fastx_quality_filter_24illumina_output/'${f2[$i]}'_fastqfilterout.fastq | awk '"'{if(NR%4==2) print length("'$1'")}'"' | sort -n | uniq -c > /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/4_readlength_distribution_24illumina_output/'${f2[$i]}'_readlengthdistfull.txt' >> base.sh
    
 echo 'fastx_collapser -Q33 -i /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/2_fastx_quality_filter_24illumina_output/'${f2[$i]}'_fastqfilterout.fastq -o /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/5_collapsed_fasta_24illumina_output/'${f2[$i]}'_collapsed.fa -v > /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/5_collapsed_fasta_24illumina_output/verbose/'${f2[$i]}'_verbose.txt' >> base.sh

 echo 'qstat -f ${PBS_JOBID}'  >> base.sh
 
 
 echo 'cd /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/' >> base.sh
 echo 'rm' `echo ${f1[$i]}|cut -d. -f1,2`  >> base.sh
 
 qsub base.sh
 
 done
