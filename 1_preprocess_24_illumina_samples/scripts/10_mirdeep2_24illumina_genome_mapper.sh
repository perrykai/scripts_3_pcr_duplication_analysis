#============================================
#  File: 10_mirdeep2_24illumina_genome_mapper.sh
#  Directory Code:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/scripts
#  Date:  8/16/16 #UPDATE: Using miRDeep2/0.0.5 due to errors received when trying to run miRDeep2/0.0.7
#  Description:  This code maps the 24 smallRNA libraries to the sus scrofa reference genome 10.2.79 using miRDeep2 mapper module.
#                1. Uses config file to feed all 24 libraries to mapper module 
#                2. Uses bowtie-index of reference genome to map reads using: bowtie –f –n 0 –e 80 –l 18 –a –m 5 –best –strata
#                3. Outputs a file of processed reads, and an .arf format alignment file for use with miRDeep2 core module
#============================================
#  Input File Directory:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/6_prepped_fasta_24illumina_output
#                         /mnt/research/pigeqtl/indexes/allDNA/Sus_scrofa10.2.79/bowtie1_index/

#  Input File(s):   1_config_24illumina_mapper.txt
#                   (6) genome index files, each using the prefix ssc_10_2_79_bowtie1_ref

#  Output File Directory: /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/6_prepped_fasta_24illumina_output
#                         /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/8_mirdeep2_genome_mapper_24illumina_output/

#  Output File(s): 6_prepped_fasta_24illumina_output/mapper.log
#                  6_prepped_fasta_24illumina_output/bowtie.log
#                  8_mirdeep2_genome_mapper_24illumina_output/1_24illumina_reads_processed.fa
#                  8_mirdeep2_genome_mapper_24illumina_output/2_24illumina_reads_mapped.arf

#============================================
#  Module used: miRDeep2/0.0.5
#============================================

#!/bin/sh  -login
#PBS -l nodes=1:ppn=1:intel16,walltime=1:00:00,mem=2Gb
#PBS -N 9_miRDeep2_24illumina_genome_mapper
#PBS -j oe
#PBS -o /mnt/research/pigeqtl/analyses/microRNA/OutputsErrors/3_preprocess_24illumina_fastq/
#PBS -m a
#PBS -M perrykai@msu.edu

cd /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/6_prepped_fasta_24illumina_output/

module load miRDeep2/0.0.5

mapper.pl /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/6_prepped_fasta_24illumina_output/1_config_24illumina_mapper.txt -d -c -p /mnt/research/pigeqtl/indexes/allDNA/Sus_scrofa10.2.79/bowtie1_index/ssc_10_2_79_bowtie1_ref -s /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/8_mirdeep2_genome_mapper_24illumina_output/1_24illumina_reads_processed.fa -t /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/8_mirdeep2_genome_mapper_24illumina_output/2_24illumina_reads_mapped.arf -m -n -v

qstat -f ${PBS_JOBID}

