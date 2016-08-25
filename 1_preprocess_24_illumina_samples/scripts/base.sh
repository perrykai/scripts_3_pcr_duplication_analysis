#!/bin/sh  -login
#PBS -l nodes=1:ppn=1:intel14,walltime=2:00:00,mem=2Gb
#PBS -N illumina_sample_1662
#PBS -j oe
#PBS -o /mnt/research/pigeqtl/analyses/microRNA/OutputsErrors/3_preprocess_24illumina_fastq/
#PBS -m a
#PBS -M perrykai@msu.edu
module load cutadapt/1.4.1
module load FASTX/0.0.14
cd /mnt/research/pigeqtl/datasets/microrna/
cp 1662_CCAACA_L007_R1_001.fastq.gz /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/
cd /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/
gunzip 1662_CCAACA_L007_R1_001.fastq.gz
cd /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output
fastx_quality_stats -i /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1662_CCAACA_L007_R1_001.fastq -o /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/3_trimfilter_statistics_24illumina_output/1_raw_fastq_quality_assessment/1662_pretrimqstats.txt -Q33
cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC -m 18 -M 30 --info-file=/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/1_cutadapt_adaptor_trim_size_filter_24illumina_output/mirinfo/1662_info --rest-file=/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/1_cutadapt_adaptor_trim_size_filter_24illumina_output/mirrest/1662_rest --too-long-output=/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/1_cutadapt_adaptor_trim_size_filter_24illumina_output/mirtoolong/1662_toolong /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1662_CCAACA_L007_R1_001.fastq -o /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/1_cutadapt_adaptor_trim_size_filter_24illumina_output/1662_cutadaptout.fastq > /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/1_cutadapt_adaptor_trim_size_filter_24illumina_output/mirsummary/1662_summary
fastq_quality_filter -q 30 -p 50 -i /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/1_cutadapt_adaptor_trim_size_filter_24illumina_output/1662_cutadaptout.fastq -o /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/2_fastx_quality_filter_24illumina_output/1662_fastqfilterout.fastq -v > /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/2_fastx_quality_filter_24illumina_output/qualityreport/1662_report -Q33
fastx_quality_stats -i /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/2_fastx_quality_filter_24illumina_output/1662_fastqfilterout.fastq -o /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/3_trimfilter_statistics_24illumina_output/2_trimmed_filtered_fastq_quality_assessment/1662_posttrimqstats.txt -Q33
cat /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/2_fastx_quality_filter_24illumina_output/1662_fastqfilterout.fastq | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/4_readlength_distribution_24illumina_output/1662_readlengthdistfull.txt
fastx_collapser -Q33 -i /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/2_fastx_quality_filter_24illumina_output/1662_fastqfilterout.fastq -o /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/5_collapsed_fasta_24illumina_output/1662_collapsed.fa -v > /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/5_collapsed_fasta_24illumina_output/verbose/1662_verbose.txt
qstat -f ${PBS_JOBID}
cd /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/
rm 1662_CCAACA_L007_R1_001.fastq
