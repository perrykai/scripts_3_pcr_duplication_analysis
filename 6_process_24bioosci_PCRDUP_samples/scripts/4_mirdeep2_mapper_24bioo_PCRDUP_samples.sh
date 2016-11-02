#!/bin/sh  -login
#PBS -l nodes=1:ppn=1:intel16,walltime=2:00:00,mem=3Gb
#PBS -N 4_mirdeep2_mapper_24bioo_PCRDUP_samples
#PBS -j oe
#PBS -o /mnt/research/pigeqtl/analyses/microRNA/OutputsErrors/4_24bioo_PCRDUP_samples/
#PBS -m a
#PBS -M perrykai@msu.edu

#' **Script:** `4_mirdeep2_mapper_24bioo_PCRDUP_samples.sh`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/scripts`
#' 
#' **Date:**  9/22/16
#' 
#' **Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/1_mirdeep2_mapper_input_24bioosci_PCRDUP_samples/`
#' 
#' **Input File(s):** `1_config_24bioosci_PCRDUP_mapper.txt`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/2_mirdeep2_mapper_24bioosci_PCRDUP_output`
#' 
#' **Output File(s):** 
#' 
#' 1. `1_24bioosci_PCRDUP_reads_processed.fa`
#' 2. `2_24bioosci_PCRDUP_reads_mapped.fa`
#' 3. `mapper.log` (located in input directory)
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 
#' ## Objectives
#' 
#' Description:  This code maps the 24 smallRNA Bioo Scientific-prepped libraries containing PCR duplciates to the sus scrofa reference genome 10.2.79 using miRDeep2 mapper module.
#' 
#' 1. Uses config file to feed all 174 libraries to mapper module
#' 2. Uses bowtie-index of reference genome to map reads using: bowtie –f –n 0 –e 80 –l 18 –a –m 5 –best –strata
#' 3. Outputs a file of processed reads, and an .arf format alignment file for use with miRDeep2 core module
#' 
#' ## Install libraries
#' 
#' Module used: miRDeep2/0.0.5
#' 
#' ## Analysis

cd /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/1_mirdeep2_mapper_input_24bioosci_PCRDUP_samples/

module load miRDeep2/0.0.5

mapper.pl /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/1_mirdeep2_mapper_input_24bioosci_PCRDUP_samples/1_config_24bioosci_PCRDUP_mapper.txt -d -c -p /mnt/research/pigeqtl/indexes/allDNA/Sus_scrofa10.2.79/bowtie1_index/ssc_10_2_79_bowtie1_ref -s /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/2_mirdeep2_mapper_24bioosci_PCRDUP_output/1_24bioosci_PCRDUP_reads_processed.fa -t /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/2_mirdeep2_mapper_24bioosci_PCRDUP_output/2_24bioosci_PCRDUP_reads_mapped.arf -m -n -v

qstat -f ${PBS_JOBID}

# This is a comment to see if the stinking save thing works