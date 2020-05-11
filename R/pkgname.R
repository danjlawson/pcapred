#' @title pcapred: Produce PC Eigenvector Predictions From Plink Input Data And Loadings
#'
#' @details
#' This is the software for "pcapred", written by Daniel Lawson (dan.lawson@bristol.ac.uk) and supporting a paper by Aliya Sarmanova, Tim Morris and Daniel Lawson (to appear).
#'
#' See \url{https://github.com/danjlawson/pcapred} for details.
#'
#' 
#'
#' @examples
#' \dontrun{
#' library("pcapred")
#' 
#' mydata=pcapred.ref::onek_genomes_tiny() # Gets the file location of
#'   # the tiny bim/bed/fam data included in pcapred's data package, pcapred.ref
#' 
#' dat=readbed(mydata) # Read "your" data
#' dat=mergeref(dat) # Merge with the reference
#'                 # (using the included standard reference of
#'                 # 18 UK Biobank Pcs by default)
#' pred=predictpcs(dat)  # Predict the first 18 UK Biobank PCs
#' writepred("projected.eigenvals",dat,pred) # Write output in plink --covar format
#' }
#' 
#'
#' @keywords internal
#' @references To follow
"_PACKAGE"
#> [1] "_PACKAGE"
