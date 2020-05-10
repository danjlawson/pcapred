#' @title Create direct SNP weights from a plink-like (scaled) reference
#' @description
#' PCA can be represented in an a number of different ways, provided that the structure
#' $$Y = U D V^T$$
#' holds, where $D$ is a diagonal matrix of the Eigenvalues. In particular, plink and FastPCA use \eqn{U' = U \sqrt{D}} using $Y=U'{V'}^T$.
#'
#' This function converts $V'$ into the canonical form, which also corrects for the factor \eqn{\sqrt{2}} from diploid genomes.
#' @param snpweights A weighting matrix of dimension L by K.
#' @param lambda A n eigenvalue vector of length K
#' @return a canonical scaled version of the snpweights matrix
#' @export
transform_plinkref=function(snpweights,lambda){
  ## Weight each SNP by 1/sqrt(lambda)
  if(class(snpweights)=="matrix"){
    if(length(lambda)!=dim(snpweights)[2]) stop("Error in transform.plinkref: must provide 1 eigenvalue per snpweight column")
    ret=snpweights %*% diag(1/sqrt(lambda))
  }else{
    ret=snpweights/sqrt(lambda)
  }
  ret/sqrt(2)
}

#' @title Create direct SNP weights from a flashpca-like (scaled) reference
#' @description
#' PCA can be represented in an a number of different ways, provided that the structure
#' $$Y = U D V^T$$
#' flashpca uses this canonical form but the weights are out by a factor \eqn{\sqrt{2}} due to diploid genomes.
#' @param snpweights A weighting matrix of dimension L by K
#' @param ... Ignored parameters (for using the same call structure as other transforms)
#' @return a canonical scaled version of the snpweights matrix
#' @export
#' @seealso transform.plinkref
transform_flashref=function(snpweights,...){
  ret=snpweights/sqrt(2)
}

#' @title Create direct SNP weights from a reference
#' @description
#' PCA can be represented in an a number of different ways, provided that the structure
#' $$Y = U D V^T$$
#' flashpca uses this canonical form but the weights are out by a factor \eqn{\sqrt{2}} due to diploid genomes. Plink-like loadings instead scale U by the square root of the eigenvalues.
#' @param snpweights A weighting matrix of dimension L by K
#' @param pcvals Eigenvalues from the original PCA
#' @param mode "plink" pr "flashpca"
#' @return a canonical scaled version of the snpweights matrix
#' @export
#' @seealso transform.plinkref
transform_pcs=function(snpweights,pcvals,mode){
    if(mode=="plink"){
        res=transform_plinkref(t(snpweights),pcvals)
    }else if(mode=="flashpca"){
        res=transform_flashref(t(snpweights))
    }
    res
}

#' @title Create predictions for eigenvalues from plink-scaled data
#' @description
#' Make a genome-wide prediction based on plink scalings
#' @param ndat normalised data matrix
#' @param res snp weighting matrix
#' @param L the number of SNPs
#' @return A matrix of size N by K
#' @export
pred_plink_simple=function(ndat,res,L){
  pred=ndat %*% res/L
}

#' @title Create predictions for eigenvalues from flashpca-scaled data
#' @description
#' Make a genome-wide prediction based on plink scalings
#' @param ndat normalised data matrix
#' @param res snp weighting matrix
#' @param L the number of SNPs
#' @return A matrix of size N by K
#' @export
pred_flash_simple=function(ndat,res,L){
  pred=ndat %*% res/sqrt(L)
  return(pred)
}

#' @title Create predictions for eigenvalues from flashpca-scaled data
#' @description
#' Make a genome-wide prediction based on plink scalings
#' @return A matrix of size N by K
#' @param keepdat a matrix of presence indicators 
#' @param ndat normalised data matrix
#' @param res snp weighting matrix
#' @param Li The number of non-missing inds for each locus
#' @export
pred_plink_missing=function(keepdat,ndat,res,Li){
  pred=t(sapply(1:dim(ndat)[1],function(i){
    ndat[i,keepdat[,i],drop=FALSE] %*% res[keepdat[,i],]/Li[i]
  }))
  return(pred)
}

#' @title Create predictions for eigenvalues from flashpca-scaled data
#' @description
#' Make a genome-wide prediction based on plink scalings
#' @param keepdat a matrix of presence indicators 
#' @param ndat normalised data matrix
#' @param res snp weighting matrix
#' @param L The number of non-missing SNPs
#' @return A matrix of size N by K
#' @export
pred_flash_missing=function(keepdat,ndat,res,L){
  pred=t(sapply(1:dim(ndat)[1],function(i){
    ndat[i,keepdat[,i],drop=FALSE] %*% res[keepdat[,i],]/sqrt(L)
  }))
  return(pred)
}

#' @title Create predictions for eigenvalues when there is no missing data
#' @description
#' Make a genome-wide prediction based on complete data, which makes the computation faster
#' @param ndat normalised data matrix
#' @param res snp weighting matrix
#' @param L the number of SNPs
#' @param mode either "flashpca" or "plink"
#' @return A matrix of size N by K
#' @export
pred_simple=function(ndat,res,L,mode){
  if(mode=="plink") return(pred_plink_simple(ndat,res,L))
  else if(mode=="flashpca") return(pred_flash_simple(ndat,res,L))
  else stop("error: invalid mode")
}

#' @title Create predictions for eigenvalues when there is missing data
#' @description
#' Make a genome-wide prediction based on incomplete data
#' @param ndat normalised data matrix
#' @param res snp weighting matrix
#' @param mode either "flashpca" or "plink"
#' @return A matrix of size N by K
#' @export
pred_missing=function(ndat,res,mode){
  keepdat=apply(ndat,1,function(x)!is.na(x))
  Li=colSums(keepdat)
  L = mean(Li)
  if(mode=="plink") return(pred_plink_missing(keepdat,ndat,res,Li))
  else if(mode=="flashpca") return(pred_flash_missing(keepdat,ndat,res,L)) ## TODO: These seem the same!
  else stop("error: invalid mode")
}

#' @title Normalise SNPs according to a reference frequency under a binomial model
#' @description
#' Normalise a SNP count matrix by an external reference frequency distribution by mean and standard deviation predicted under a binomial model.
#' @param rawdatmat A matrix of SNP counts of the alternative allele; of dimension N by L
#' @param f A snp frequency vector of length L
#' @param ploidy ploidy of the data (2=diploid) ; maximum that rawdatmat can take
#' @param minfreq a thresholding for very raw SNPs in the reference, to prevent them from being over-weighted
#' @param meanimpute (default FALSE) whether to mean impute missing data
#' @return A matrix of size N by L which is rawdatmat zero meaned and with theoretical s.d. 1
#' @export
normalise=function(rawdatmat,f,ploidy=2,minfreq=NULL,meanimpute=FALSE){
    if(all(is.null(minfreq))) minfreq=1/length(f)
    f[f>1-minfreq]=1-minfreq
    f[f<minfreq]=minfreq
    fnorm=sqrt(f*(1-f))
    mu=ploidy*f
    
    ndat=t(apply(rawdatmat,1,function(x){
        ret=(x-mu)/fnorm
        ret
    }))
    if(meanimpute) ndat[is.na(ndat)]=0
    ndat
}

#' @title Normalise SNPs according to a reference frequency under a binomial model
#' @description
#' Normalise a SNP count vector for a single SNP by an external reference frequency distribution by mean and standard deviation predicted under a binomial model.
#' @param x A vector of Individuals SNP data of length N
#' @param f The snp frequency for that SNP
#' @param ploidy ploidy of the data (2=diploid) ; maximum that rawdatmat can take
#' @param minfreq a thresholding for very raw SNPs in the reference, to prevent them from being over-weighted
#' @param meanimpute (default FALSE) whether to mean impute missing data
#' @return A vector of length = length(x)
#' @export
normalisesnp=function (x, f, ploidy = 2, minfreq = 0.001,meanimpute=FALSE)
{
    f[f > 1 - minfreq] = 1 - minfreq
    f[f < minfreq] = minfreq
    fnorm = sqrt(f * (1 - f))
    mu = ploidy * f
    ndat = (x - mu)/fnorm
    if(meanimpute) ndat[is.na(ndat)]=0
    ndat
}
