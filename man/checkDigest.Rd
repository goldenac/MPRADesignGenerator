% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{checkDigest}
\alias{checkDigest}
\title{CHECK FOR ALL THREE DIGESTION SITES}
\usage{
checkDigest(df_to_test, ez1, ez2, ez3)
}
\arguments{
\item{dataframe}{with sequences to test}
}
\value{
new dataframe containing sequences which had a digestion site
}
\description{
Check all sequences in a dataframe for digestion sites
(SfiI, KpnI, XbaI)
(column with sequences tested should be named short_seq)
}
