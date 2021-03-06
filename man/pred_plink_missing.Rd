% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mappingfunctions.R
\name{pred_plink_missing}
\alias{pred_plink_missing}
\title{Create predictions for eigenvalues from flashpca-scaled data}
\usage{
pred_plink_missing(keepdat, ndat, res, Li)
}
\arguments{
\item{keepdat}{a matrix of presence indicators}

\item{ndat}{normalised data matrix}

\item{res}{snp weighting matrix}

\item{Li}{The number of non-missing inds for each locus}
}
\value{
A matrix of size N by K
}
\description{
Make a genome-wide prediction based on plink scalings
}
