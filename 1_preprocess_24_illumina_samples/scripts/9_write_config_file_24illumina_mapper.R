#============================================
#  File: 8_write_config_file_24illumina_mapper.R
#  Directory Code:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/scripts
#  Date:  8/15/16
#  Description:  This code writes the config file for the 24 samples being analyzed by miRDeep2 mapper and core modules
#============================================
#  Input File Directory:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/6_prepped_fasta_24illumina_output/

#  Input File(s):  xxxx_mapper_input.fa file for each sample

#  Output File Directory: /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/6_prepped_fasta_24illumina_output/

#  Output File(s): 1_config_24illumina_mapper.txt
#============================================
#  Module used: R/3.1.0
#============================================
R
setwd("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/6_prepped_fasta_24illumina_output/")
fn<-as.vector(list.files(pattern="_mapper_input.fa"))               #create a vector of the file names

fn
#  [1] "1034_mapper_input.fa" "1058_mapper_input.fa" "1080_mapper_input.fa"
#  [4] "1096_mapper_input.fa" "1116_mapper_input.fa" "1134_mapper_input.fa"
#  [7] "1154_mapper_input.fa" "1170_mapper_input.fa" "1194_mapper_input.fa"
# [10] "1240_mapper_input.fa" "1278_mapper_input.fa" "1300_mapper_input.fa"
# [13] "1426_mapper_input.fa" "1434_mapper_input.fa" "1458_mapper_input.fa"
# [16] "1484_mapper_input.fa" "1502_mapper_input.fa" "1512_mapper_input.fa"
# [19] "1534_mapper_input.fa" "1580_mapper_input.fa" "1594_mapper_input.fa"
# [22] "1640_mapper_input.fa" "1644_mapper_input.fa" "1662_mapper_input.fa"

length(fn)
# [1] 24

numf<-seq(1,length(fn),1)           #create a vector of numbers the same length as the number of file names
numf
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24

dig<-sprintf("%03d",numf)           #sprintf function allows for the creation of a 3-digit numeric identifier
dig
# [1] "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012"
# [13] "013" "014" "015" "016" "017" "018" "019" "020" "021" "022" "023" "024"

config<-cbind(fn,dig)               #bind the two columns together
head(config)
#      fn                     dig
# [1,] "1034_mapper_input.fa" "001"
# [2,] "1058_mapper_input.fa" "002"
# [3,] "1080_mapper_input.fa" "003"
# [4,] "1096_mapper_input.fa" "004"
# [5,] "1116_mapper_input.fa" "005"
# [6,] "1134_mapper_input.fa" "006"

write.table(config, file = "1_config_24illumina_mapper.txt", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE) #write the config object to a text file
q()
n
