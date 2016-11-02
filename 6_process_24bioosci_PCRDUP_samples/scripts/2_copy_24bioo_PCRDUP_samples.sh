#' **Script:** `2_copy_24bioo_PCRDUP_samples.sh`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/scripts`
#' 
#' **Date:**  9/21/16
#' 
#' **Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/1_preprocess_fastq_files/1_preprocess_fastq_files_output/5_collapsed_fasta_output`
#' 
#' **Input File(s):** `xxxx_collapsed.fa`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/`
#' 
#' **Output File(s):** `xxxx_24bioosci_collapsed.fa`
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Analysis](#analysis)
#' 3. [Save data](#save-data)
#' 
#' ## Objectives
#' 
#' The objective of this script is to copy the 24 xxxx_collapsed.fa files from the 174 library analysis to this directory for analysis of the PCR-duplicate containing libraries.
#' 
#' ## Analysis

cd /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/scripts

f1=(`ls /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/5_collapsed_fasta_24illumina_output -1|grep '_collapsed.fa' |cut -f1,2 -d'.'` )
f2=(`ls /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/5_collapsed_fasta_24illumina_output -1|grep '_collapsed.fa' |cut -f1 -d_`)

echo '#!/bin/sh  -login' > base.sh
echo '#PBS -l nodes=1:ppn=1,walltime=00:10:00,mem=1Gb' >> base.sh
echo '#PBS -N 2_copy_24bioo_PCRDUP_samples' >> base.sh
echo '#PBS -j oe' >> base.sh
echo '#PBS -o /mnt/research/pigeqtl/analyses/microRNA/OutputsErrors/4_24bioo_PCRDUP_samples/' >> base.sh
echo '#PBS -m a' >> base.sh
echo '#PBS -M perrykai@msu.edu' >> base.sh

for (( i = 0 ; i < ${#f1[@]} ; i++ )) do
 echo 'cp /mnt/research/pigeqtl/analyses/microRNA/1_preprocess_fastq_files/1_preprocess_fastq_files_output/5_collapsed_fasta_output/'${f1[$i]}' /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/'${f2[$i]}'_24bioosci_collapsed.fa' >> base.sh
done

qsub base.sh

rm base.sh
