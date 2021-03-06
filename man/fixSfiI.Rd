% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{fixSfiI}
\alias{fixSfiI}
\title{FIX SfiI DIGESTION SITE}
\usage{
fixSfiI(sequence, coord1, coord2)
}
\arguments{
\item{145bp}{sequence containing site, coordinate of first base of
restriction site, coordinate of second base of restriction site.}
}
\value{
corrected sequence
}
\description{
Changes SfiI site (GGCCNNNNNGGCC -> GCGCNNNNNGGCC)
Will only first first digestion site in sequence.
}
