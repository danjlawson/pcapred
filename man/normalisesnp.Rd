% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mappingfunctions.R
\name{normalisesnp}
\alias{normalisesnp}
\title{Normalise SNPs according to a reference frequency under a binomial model}
\usage{
normalisesnp(x, f, ploidy = 2, minfreq = 0.001, meanimpute = FALSE)
}
\arguments{
\item{x}{A vector of Individuals SNP data of length N}

\item{f}{The snp frequency for that SNP}

\item{ploidy}{ploidy of the data (2=diploid) ; maximum that rawdatmat can take}

\item{minfreq}{a thresholding for very raw SNPs in the reference, to prevent them from being over-weighted}

\item{meanimpute}{(default FALSE) whether to mean impute missing data}
}
\value{
A vector of length = length(x)
}
\description{
Normalise a SNP count vector for a single SNP by an external reference frequency distribution by mean and standard deviation predicted under a binomial model.
}
