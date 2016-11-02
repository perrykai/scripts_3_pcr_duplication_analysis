#!/bin/sh  -login
#PBS -l nodes=1:ppn=1:intel16,walltime=12:00:00,mem=10Gb
#PBS -N 5_24bioosci_PCRDUP_miRDeep2_core
#PBS -j oe
#PBS -o /mnt/research/pigeqtl/analyses/microRNA/OutputsErrors/4_24bioo_PCRDUP_samples/
#PBS -m abe
#PBS -M perrykai@msu.edu

#' **Script:** `5_mirdeep2_core_24bioosci_PCRDUP_samples.sh`
#' 
#' **Directory of Code:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/scripts`
#' 
#' **Date:**  9/21/16
#' 
#' **Input File Directory:**  `/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/2_mirdeep2_mapper_24bioosci_PCRDUP_output`
#' 
#' **Input File(s):** 
#' 
#' 1. `2_24bioosci_PCRDUP_reads_mapped.arf`
#' 2. `1_24bioosci_PCRDUP_reads_processed.fa
#' 
#' **Output File Directory:** `3_mirdeep2_core_output_24bioosci_PCRDUP_samples`
#' 
#' **Output File(s):** 
#' 
#' 1. `miRNAs_expressed_all_samples_[timestamp].csv`
#' 2. `expression_[timestamp].html`
#' 3. `expression_analyses directory containing quantifier module outputs`
#' 4. `dir_prepare_signature_timestamp directory`
#' 5. `result_[timestamp].csv, .html, and .bed files`
#' 6. `pdfs directory`  
#' 7. `1_report_24bioosci_PCRDUP_mirdeep2_core.log`
#' 8. `error.log`
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 5. [Save data](#save-data)
#' 
#' ## Objectives
#'
#' Description:  This code conducts the quantification and prediction steps for the 24 illumina smallRNA libraries using miRDeep2 core module.
#'                1. Runs mirdeep2 quantifier module on known miRNAs, using miRBase mature and precursor sequences as reference
#'                2. Excises potential miRNA precursor sequences from genome, using read mappings as a guide (see Friedlander et al 2012 for details)
#'                3. Uses bowtie to map sequence reads to the bowtie index of excised precursors
#'                4. Uses RNAfold to predict the secondary structures of the precursors, including Randfold p-value calculation
#'                5. Core algorithm evaluates structure and significance of potential miRNA precursors

#' ## Install libraries
#'
#' Module used: miRDeep2/0.0.5
#' 
#' ## Analysis

module load miRDeep2/0.0.5

cd /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/3_mirdeep2_core_output_24bioosci_PCRDUP_samples

miRDeep2.pl /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/2_mirdeep2_mapper_24bioosci_PCRDUP_output/1_24bioosci_PCRDUP_reads_processed.fa /mnt/research/pigeqtl/analyses/microRNA/reference_sequences/Sscrofa79ref.dna.fa /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/2_mirdeep2_mapper_24bioosci_PCRDUP_output/2_24bioosci_PCRDUP_reads_mapped.arf /mnt/research/pigeqtl/analyses/microRNA/reference_sequences/ssc_mature_mir.fa /mnt/research/pigeqtl/analyses/microRNA/reference_sequences/hsa_mature_mir.fa /mnt/research/pigeqtl/analyses/microRNA/reference_sequences/ssc_hairpin2.fa 2> 1_report_24bioosci_PCRDUP_mirdeep2_core.log

qstat -f ${PBS_JOBID}

