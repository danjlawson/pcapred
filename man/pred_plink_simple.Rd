% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mappingfunctions.R
\name{pred_plink_simple}
\alias{pred_plink_simple}
\title{Create predictions for eigenvalues from plink-scaled data}
\usage{
pred_plink_simple(ndat, res, L)
}
\arguments{
\item{ndat}{normalised data matrix}

\item{res}{snp weighting matrix}

\item{L}{the number of SNPs}
}
\value{
A matrix of size N by K
}
\description{
Make a genome-wide prediction based on plink scalings
}
