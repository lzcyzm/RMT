\name{xls2Grangeslist}
\alias{xls2Grangeslist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
translate excel to Grangeslist
%%  ~~function to do ... ~~
}
\description{
read the excel, translate the information of excel to the style of Grangelist
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xls2Grangeslist(filepath)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filepath}{ the path of excel
%%     ~~Describe \code{filepath} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function outputs the information of peak in a Grangelist with gene information
}

\author{
Lian Liu <liulian19860905@163.com>
}

\examples{
library("RMT")
filepath=system.file("extdata", "peak.xls", package="exomePeak")
anno=xls2Grangeslist(filepath)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RNA epigenetics }

