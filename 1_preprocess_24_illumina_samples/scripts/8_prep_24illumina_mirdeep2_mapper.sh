#!/bin/sh  -login
#PBS -l nodes=1:ppn=1:intel16,walltime=00:10:00,mem=2Gb
#PBS -N 8_prep_24illumina_mirdeep2_mapper
#PBS -j oe
#PBS -o /mnt/research/pigeqtl/analyses/microRNA/OutputsErrors/3_preprocess_24illumina_fastq/
#PBS -m a
#PBS -M perrykai@msu.edu

#==============================================================================
#============================================
#  File:  8_prep_24illumina_mirdeep2_mapper.sh
#  Directory Code:  /mnt/research/pigeqtl/analyses/microRNA/1_preprocess_fastq_files/scripts/
#  Date:  1/9/15
#  Description:  This code edits the sequence ids for the 24 illumina libraries for compatability with the miRDeep2 mapper module.
#                Makes the sequence reads follow the format >seq_id_xreadcount on the header lines, then saves them to the directory to be used for the mapper module and zips the non-prepped fasta file.
#                The output .fa files will be used for the mapping and miRNA quantification steps, and the PCR duplication analysis.
#============================================
#  Input File Directory:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/5_collapsed_fasta_24illumina_output

#  Input File(s):   xxxx_collapsed.fa
 
#  Output File Directory: /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/6_shortread_prepped_fasta_24illumina_output/

#  Output File(s): xxxx_mapper_input.fa (One output file for each sample)
#============================================
#  Example Output (If Applicable):

#  >seq_1_x54014
#  TGGAATGTAAAGAAGTATGTAT
#  >seq_2_x33188
#  TGGAATGTAAAGAAGTATGTAC
#  >seq_3_x29438
#  TGGAATGTAAGGAAGTGTGTGA
#  >seq_4_x23873
#  TGGAATGTAAAGAAGTATGTATT
#  >seq_5_x18417
#  TGGAATGTAAAGAAGTATGT

#============================================
#Softwares used: sed and gzip (command line)
#============================================

cd /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/5_collapsed_fasta_24illumina_output/

f1=(`ls *_collapsed.fa`)
f2=(`ls *_collapsed.fa | cut -d '_' -f1`)

# The first sed command, sed 's/^>/&'seq_'/' inserts the seq_ needed at the beginning of the header line.
# The second sed command substitutes the hyphon for the '_x' needed for the header line

for (( i = 0 ; i < ${#f1[@]} ; i++ )) do

	sed 's/^>/&'seq_'/' ${f1[$i]} | sed -r 's/'-'/'_x'/g' > /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/6_prepped_fasta_24illumina_output/"${f2[i]}"_mapper_input.fa

	gzip ${f1[$i]}

done

qstat -f ${PBS_JOBID}