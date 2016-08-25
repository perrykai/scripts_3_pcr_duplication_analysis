    #============================================
    #  File:  4_readlengthdist_histogram.R
    #  Directory Code:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/scripts
    #  Date:  1/8/15 #UPDATED 8/12/16
    #  Description:  This code creates a histogram of the read length distributions of the preliminary 24 smallRNA libraries.
    #                
    #============================================
    #  Input File Directory:  /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/4_readlength_distribution_24illumina_output
    
    #  Input File(s):  xxxx_readlengthdistfull.txt
    
    #  Output File Directory: /mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/4_readlengthdist_histogram_24illumina_plot_output/
    
    #  Output File(s): 3_readlengthdist_histogram_24illumina.pdf (visual plot of read length distribution)

    #  Run interactively in R/3.1.0
    #============================================    
    
    rm(list=ls())
    
    setwd("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/1_preprocess_24illumina_output/4_readlength_distribution_24illumina_output/")
    
    filelist<-grep("readlengthdistfull.txt",list.files(),value=T)
    head(filelist)
    
    #Create the empty vector
    counts<-rep(0,50)
    
    #Create a for loop to run all 24 files
    for (i in filelist){
      cnt<-read.table(i)
      id<-cnt[,2]
      fr<-cnt[,1]
      counts[id]<-counts[id]+fr
    }
    
    
    #Now, counts is the total frequency of each read length for the 24 files.
    counts
    #  [1]        0        0        0        0        0        0        0        0
    #  [9]        0        0        0        0        0        0        0        0
    # [17]        0  1948256  1991326  4799054 11897766 63719774 22537908  6154642
    # [25]  1650429  1283827  1465194  1502118  2082391  3167309        0        0
    # [33]        0        0        0        0        0        0        0        0
    # [41]        0        0        0        0        0        0        0        0
    # [49]        0        0
    
    counts2<-counts/sum(counts)
    counts2 
    #  [1] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #  [7] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    # [13] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.01568644
    # [19] 0.01603322 0.03863973 0.09579522 0.51304168 0.18146465 0.04955429
    # [25] 0.01328848 0.01033677 0.01179705 0.01209435 0.01676643 0.02550168
    # [31] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    # [37] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    # [43] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    # [49] 0.00000000 0.00000000
    
    #Then, use barplot to visually represent the distribution of read lengths in this complete dataset.
    pdf("/mnt/research/pigeqtl/analyses/microRNA/3_pcr_duplication_analysis/1_preprocess_24_illumina_samples/4_readlengthdist_histogram_24illumina_plot_output/3_readlengthdist_histogram_24illumina.pdf")
    par(oma=c(2,2,2,2))
    barplot(counts2[17:31],names.arg=17:31,space=0,axes=T,
            main="Read Length Distribution 24 Libraries",
            xlab="Read Length (nt)",ylab="Proportional Frequency", ylim=c(0,0.7),
            col="deepskyblue",border=TRUE,
            cex.main=1.8,cex.lab=1.5,cex.axis=1.2)
    dev.off()
    q()
    n