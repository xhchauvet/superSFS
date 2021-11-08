# superSFS
This is a tool for speculation of ancestral allel, calculation of sfs and drawing its bar plot. It is easy-to-use and runing fast. What you should prepare is the phased vcf file containg the data of populations you intrested and the outgroup, the outgroup name file, and the annotation file. Enjoy it!!!

It has four models:
    0：Using all function, from original vcf data to sfs barplot
    1: Only speculate the ancestral allel and output new vcf file using speculated allel as reference
    2: Only count the frequency of derived allel in each snp of each population
    3: Only draw bar polt of sfs using data generated from the results of calutation of sfs

Example:

    Model 0: python superSFS 0 ogdir threshold vcfdir annodir modir coutdir plotdir group
    Model 1: python superSFS 1 ogdir threshold vcfdir outdir
    Model 2: python superSFS 2 annodir modir coutdir
    Model 3: python superSFS 3 coutdir plotdir group

Explation for each parameter:
    ogdir: direction of outgroup names file
    threshold: a number that if the sum of variant allel in outpgroup greater than it,the variant allel will be counted as ancestral allel
    vcfdir: direction of vcf data
    vannodir: direction of annotation file with sample names in first column and group name in second colum. This file should has header in first row
    vmodir: assign the output direction of generated vcf file using speculated allel as reference
    countdir: assign the output direction of calculation of derived allels for each snp in each group
    plotdir: assign the output direction of bar plot of sfs
    group: the group that you want to analysis
