% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{retrieveSeq}
\alias{retrieveSeq}
\title{RETRIEVE SEQUENCES}
\usage{
retrieveSeq(genome, chromNum, coord)
}
\arguments{
\item{genome}{specification, chromosome number, coordinate}
}
\value{
145bp sequence
}
\description{
Retrieve a 145bp sequence from BSgenome.Hsapiens.UCSC.hg19
(72bp downstream from specified base, 72b upstream from specified base)
}
