% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binary.R
\name{writebinhap}
\alias{writebinhap}
\title{Convert binary haplotype format to hex formatted output file}
\usage{
writebinhap(dat, fileroot, direction = "snpwise")
}
\arguments{
\item{dat}{an hapmatrix object as created by \code{\link{bedashap}}.}

\item{fileroot}{the file root of the outputs.}

\item{direction}{Either snpwise or hapwise, controls whether we place all values for a snp together (snpwise) or for a haplotype}
}
\value{
NULL invisibly
}
\description{
Represent a binary haplotype in efficient hex file format. Creates:
\itemize{
\item <fileroot>.binhaps: The matrix in efficient binary form
\item <fileroot>.hapinfo: A tab delimited file with N rows and p columns, giving the individual identifier for each haplotype (of which there are are Np)
\item <fileroot>.fam: Plinks' fam format
\item <fileroot>.bim: Plinks' bim format
}
}
