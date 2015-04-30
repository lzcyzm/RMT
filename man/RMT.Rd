\name{RMT}
\alias{RMT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This is the main function of RNA modification tools R-package.
}
\description{
This is the main function of RNA modification tools R-package.
Affinity-based epitranscriptome sequencing data: exome-based peaks and compute reads counts, rpkm and fold change. 

1. Accept and statistically supports multiple biological replicates

2. Remove possible PCR artifacts and mapping ambiguity caused by multi-reads (reads that can be mapped to multiple genomic locations)

3. exome-based peaks and compute reads counts, rpkm and fold change

4. Optionally output the intermediate results in Rdata format
}
\usage{
RMT(INPUT_BAM, IP_BAM, 
          INPUT2IP = NA, 
          EXOME_OUTPUT_DIR = NA, 
          GO_OUTPUT_DIR = NA, 
          GO_EXPERIMENT_NAME = "GO_Result")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{INPUT_BAM}{ a vector of strings, which specifies a number of Input control samples from the untreated condition in bam format
%%     ~~Describe \code{INPUT_BAM} here~~
}
  \item{IP_BAM}{ a vector of strings, which specifies a number of IP samples from the untreated condition in bam format
%%     ~~Describe \code{IP_BAM} here~~
}
  \item{INPUT2IP}{ a list of number, the default is NA, if it is NA, the length between input and ip is same, and each input bam corresponds one ip bam; otherwise, the length of INPUT2IP equals the number of IP_BAM, and each row specifies the input bams which correspondence the ip bam.
%%     ~~Describe \code{INPUT2IP} here~~ 
}
  \item{EXOME_OUTPUT_DIR}{ a string, which specifies folder name generated in the output directory that contains all the results, default: "exomePeak_output"
%%     ~~Describe \code{EXOME_OUTPUT_DIR} here~~
}
  \item{GO_OUTPUT_DIR}{ a string, which specifies the GO output directory, default: 
getwd()
%%     ~~Describe \code{GO_OUTPUT_DIR} here~~
}
  \item{GO_EXPERIMENT_NAME}{ a string, which specifies folder name generated in the GO output directory that contains all the results, default: "GO_result"
%%     ~~Describe \code{GO_EXPERIMENT_NAME} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Results will be saved two folders (default: "respect_output" and "GO_result") in the specified directory (default: current working directory).
1.In the exomPeak folder, store the information of all peak in bed and xls format
2.In the GO_result folder, store the information of reads counts, rpkm in ip bams and input bams, fold change with ip and input in xls format
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Lian Liu <liulian19860905@163.com>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# the RESPECT R-package has two main functions:
# 1. exome-based peak
# 2. compute reads counts, rpkm and fold change
# please feel free to contact liulian19860905@163.com for any questions

# load library and specify the parameters
library("RMT")
GENE_ANNO_GTF=system.file("extdata", "example.gtf", package="exomePeak")
f1=system.file("extdata", "IP1.bam", package="RMT")
f2=system.file("extdata", "IP2.bam", package="RMT")
f3=system.file("extdata", "IP3.bam", package="RMT")
f4=system.file("extdata", "IP4.bam", package="RMT")
IP_BAM=c(f1,f2,f3,f4)
f1=system.file("extdata", "Input1.bam", package="RMT")
f2=system.file("extdata", "Input2.bam", package="RMT")
f3=system.file("extdata", "Input3.bam", package="RMT")
INPUT_BAM=c(f1,f2,f3)
INPUT2IP=list(4)
INPUT2IP[[1]]=c(1,2,3)
INPUT2IP[[2]]=c(1,2,3)
INPUT2IP[[3]]=c(2,3)
INPUT2IP[[4]]=c(3)

#compute read counts and rpkm
RMT (INPUT_BAM=INPUT_BAM,IP_BAM=IP_BAM,
           INPUT2IP = INPUT2IP)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RNA epigenetics }
