#============================================
#  File: 11_mirdeep2_24illumina_core_quantify_predict.sh
#  Directory Code:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/scripts
#  Date:  8/16/16 #UPDATE: Using miRDeep2/0.0.5 due to errors received when trying to run miRDeep2/0.0.7. 
#  Description:  This code conducts the quantification and prediction steps for the 24 illumina smallRNA libraries using miRDeep2 core module.
#                1. Runs mirdeep2 quantifier module on known miRNAs, using miRBase mature and precursor sequences as reference 
#                2. Excises potential miRNA precursor sequences from genome, using read mappings as a guide (see Friedlander et al 2012 for details)
#                3. Uses bowtie to map sequence reads to the bowtie index of excised precursors
#                4. Uses RNAfold to predict the secondary structures of the precursors, including Randfold p-value calculation
#                5. Core algorithm evaluates structure and significance of potential miRNA precursors
#============================================
#  Input File Directory:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/8_mirdeep2_genome_mapper_24illumina_output
#                         /mnt/research/pigeqtl/analyses/microRNA/reference_sequences/

#  Input File(s):   8_mirdeep2_genome_mapper_24illumina_output/1_24illumina_reads_processed.fa
#                   reference_sequences/Sscrofa79ref.dna.fa   <-- THIS REFERENCE SEQUENCE CREATED 10/5/15 by removing labels after the first space from the reference genome downloaded by DV in /mnt/research/pigeqtl/references/allDNA/Sus_scrofa10.2.79 directory. See script "greprefseq.sh" in /mnt/research/pigeqtl/analyses/microRNA/reference_sequences directory
#                   8_mirdeep2_genome_mapper_24illumina_output/2_24illumina_reads_mapped.arf
#                   reference_sequences/ssc_mature_mir.fa
#                   reference_sequences/hsa_mature_mir.fa
#                   reference_sequences/ssc_hairpin2.fa

#  Output File Directory: /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/9_mirdeep2_core_quantify_predict_24illumina_output

#  Output File(s): miRNAs_expressed_all_samples_[timestamp].csv
#                  expression_[timestamp].html
#                  expression_analyses directory containing quantifier module outputs
#                  dir_prepare_signature[timestamp] directory
#                  result_[timestamp].csv, .html, and .bed files
#                  pdfs directory
#                  1_report_24illumina_mirdeep2_core.log
#                  error.log

#============================================
#  Module used: miRDeep2/0.0.5
#============================================

#!/bin/sh  -login
#PBS -l nodes=1:ppn=1:intel16,walltime=10:00:00,mem=10Gb
#PBS -N 11_24illumina_miRDeep2_core
#PBS -j oe
#PBS -o /mnt/research/pigeqtl/analyses/microRNA/OutputsErrors/3_preprocess_24illumina_fastq/
#PBS -m abe
#PBS -M perrykai@msu.edu

module load miRDeep2/0.0.5

cd /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/9_mirdeep2_core_quantify_predict_24illumina_output

miRDeep2.pl /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/8_mirdeep2_genome_mapper_24illumina_output/1_24illumina_reads_processed.fa /mnt/research/pigeqtl/analyses/microRNA/reference_sequences/Sscrofa79ref.dna.fa /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/8_mirdeep2_genome_mapper_24illumina_output/2_24illumina_reads_mapped.arf /mnt/research/pigeqtl/analyses/microRNA/reference_sequences/ssc_mature_mir.fa /mnt/research/pigeqtl/analyses/microRNA/reference_sequences/hsa_mature_mir.fa /mnt/research/pigeqtl/analyses/microRNA/reference_sequences/ssc_hairpin2.fa 2> 1_report_24illumina_mirdeep2_core.log

qstat -f ${PBS_JOBID}
