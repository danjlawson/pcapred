% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output.R
\name{test_pca}
\alias{test_pca}
\title{Test PCA predictions with access to the truth}
\usage{
test_pca(testfile, pred, outfile = NULL, verbose = TRUE)
}
\arguments{
\item{testfile}{file conmtaining the truth}

\item{pred}{The predicted PCs}

\item{outfile}{optional; the file that \code{\link{writepred}} would use; the outut is written to <outfile>.comparison}

\item{verbose}{Whetyher output commentry on the tests}
}
\value{
a list of the intermediate values, invisibly. You want the $fullres element which is also written to file.
}
\description{
Read the true PCs from a file and compare this to our predicted PCs, for error checking.
}
