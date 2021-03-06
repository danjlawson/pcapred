% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rbed.R
\name{getinds}
\alias{getinds}
\title{Extract a chunk of individuals from a binary format bed object}
\usage{
getinds(bed, x)
}
\arguments{
\item{bed}{An "rbed" object as returned from \code{\link{readbed}} or \code{\link{mergeref}}}

\item{x}{The index of the chunk to extract (NOT an individual index)}
}
\value{
a matrix of dimension 4 by (bed$no.snp or length(bed$datkeep)) containing the SNP values for each individual
}
\description{
The binary bed format stores 4 individuals in a byte. This returns the SNP details for these 4 individuals, applying any SNP filtering for a merged dataset as well as flipping any alleles.
}
\seealso{
\code{\link{getind}}, \code{\link{getsnp}}
}
