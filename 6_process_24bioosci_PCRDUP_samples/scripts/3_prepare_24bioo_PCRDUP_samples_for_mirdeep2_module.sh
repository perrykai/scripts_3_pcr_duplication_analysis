#!/bin/sh  -login
#PBS -l nodes=1:ppn=1:intel16,walltime=00:10:00,mem=1Gb
#PBS -N 3_prepare_24bioo_PCRDUP_samples_for_mirdeep2_module
#PBS -j oe
#PBS -o /mnt/research/pigeqtl/analyses/microRNA/OutputsErrors/4_24bioo_PCRDUP_samples/
#PBS -m a
#PBS -M perrykai@msu.edu

#==============================================================================
#============================================
#  File:  3_prepare_24bioo_PCRDUP_samples_for_mirdeep2_module.sh
#  Directory Code:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/scripts/
#  Date:  9/21/16
#  Description:  This code edits the sequence ids for the 24 Bioo Scientific-prepped libraries containing PCR duplicates for compatability with the miRDeep2 mapper module and trims the first and last 4 nt from each sequence.
#                Makes the sequence reads follow the format >seq_id_xreadcount on the header lines, then saves them to the directory to be used for the mapper module and zips the non-prepped fasta file.
#                The output .fa files will be used for the mapping and miRNA quantification steps for the PCR duplication analysis.
#============================================
#  Input File Directory:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples

#  Input File(s):   xxxx_24bioosci_collapsed.fa
 
#  Output File Directory: /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/1_mirdeep2_mapper_input_24bioosci_PCRDUP_samples/

#  Output File(s): xxxx__24bioosci_PCRDUP_mapper_input.fa (One output file for each sample)
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

cd /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples

f1=(`ls *_24bioosci_collapsed.fa`)
f2=(`ls *_24bioosci_collapsed.fa | cut -d '_' -f1`)

# The first sed command removes the first 4 characters from strings beginning with A-Z characters
# The second sed command removes the last 4 characters from strings ending with A-Z characters
# The third sed command, sed 's/^>/&'seq_'/' inserts the seq_ needed at the beginning of the header line.
# The last sed command substitutes the hyphon for the '_x' needed for the header line


for (( i = 0 ; i < ${#f1[@]} ; i++ )) do

	sed 's/^[A-Z]\{4\}//' ${f1[$i]} | sed 's/.\{3\}[A-Z]$//' | sed 's/^>/&'seq_'/' | sed -r 's/'-'/'_x'/g' > /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/6_process_24bioosci_PCRDUP_samples/1_mirdeep2_mapper_input_24bioosci_PCRDUP_samples/"${f2[i]}"_24bioosci_PCRDUP_mapper_input.fa

	gzip ${f1[$i]}

done

qstat -f ${PBS_JOBID}