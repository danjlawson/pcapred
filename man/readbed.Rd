% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rbed.R
\name{readbed}
\alias{readbed}
\title{Read the binary PLINK format (BED, BIM, and FAM) efficiently}
\usage{
readbed(bed, bim = NULL, fam = NULL, verbose = TRUE)
}
\arguments{
\item{bed}{the file location for the bed file, or the file root}

\item{bim}{(default NULL, meaning infer from bed) the file location for the bim file, if not using the file root}

\item{fam}{(default NULL, meaning infer from bed) the file location for the fam file, if not using the file root}

\item{verbose}{(default TRUE) Whether to report progress}
}
\value{
A list containing:
\itemize{
\item snp: The SNP matrix loaded from the BED file
\item snp.info (if !only.snp) the snp info from the bim file
\item ind.info (if !only.snp) the individual info from the fam file
}
}
\description{
Require the complete set of 3 files in the binary
PLINK format. It includes BED file, BIM file and BAM file. For more
information about the binary PLINK format, please check in the manual of PLINK.

Reads the binary file into memory but obnly stores binary format; we use accessor functions to
}
