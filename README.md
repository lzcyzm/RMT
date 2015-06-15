# RMT
RMT is now part of exomePeak package

# Install the current exomePeak package from Bioconductor 
# to ensure all the dependencies are properly installed
source("http://bioconductor.org/biocLite.R")
biocLite("GenomicAlignments")
biocLite("exomePeak")

# Install devtools
install.packages("devtools")

# Install the developmental version of exomePeak
library(devtools)
install_github('lzcyzm/exomePeak')

# the package should be ready to go now
library(exomePeak)
?RMT # please check the description of RMT 
example(RMT) # an example

# Please pay attention to the following part
# Alternatively, exomePeak package can automatically download the complete transcriptome.
# And then scan the entire transcriptome for RNA methylation sites
# It will take a long time
# RMT(INPUT_BAM=input_bam, IP_BAM=ip_bam, INPUT2IP=input2bam, GENOME="hg19")
