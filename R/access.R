
#' @title Get normalised data from an rbed object in an efficient format
#' @description
#' This function is for efficient bulk processing. It accesses data in the format in which it is stored, rather than by individual IDs.
#' @param bed The "rbed" object to read
#' @param i The index of data to read. 4 individuals are stored in each binary object. so i refers to the index of this set, i.e. i=2 gives individuals 5:8. You can give a single value or a vector. It returns fewer entries if you are at the end of the data.
#' @param meanimpute default FALSE. If normalise=TRUE, whether to replace missing values with their mean.
#' @param normalise default TRUE. Whether to return the raw data or normalise by frequency.
#' @param verbose (default TREE) whether to show progress for individual access
#' @return A matrix of data with the additional attribute "idx" referring to the indexes of the returned individuals.
#' @export
get_data=function(bed,i,meanimpute=FALSE,normalise=TRUE,verbose=TRUE){
    if(!is(bed,"rbed")) stop("bed must be of class rbed")
    if(length(i)>1){
        r=lapply(i,function(ii)get_data(bed,ii,meanimpute=meanimpute,
                                        normalise=normalise,verbose=verbose))
        ridx=lapply(r,function(x)attr(x,"idx"))
        rr=do.call("rbind",r)
        attr(rr,"idx")=do.call("c",ridx)
        return(rr)
    }
    localdata=(getinds(bed,i))
    if(is(bed,"mergedrbed")) { 
        if(normalise)localdata=normalise(localdata,
                            bed$snpstats$ALT_FREQS[bed$datkeep],
                            meanimpute=meanimpute)
    }else{
        if(all(is.null(bed$snpstats$ALT_FREQS)))stop("rbed object requires SNP frequencies when normalise=TRUE. Provide these with bed$snpstats = readsnpstats(afreqfile) or bed$snpstats=data.frame(ALT_FREQS=getfreqs(bed))")
        localdata=normalise(localdata,
                            bed$snpstats$ALT_FREQS,
                            meanimpute=meanimpute)
    }
    if(verbose)cat(paste0("Reading individual set ",i," of ",ceiling(bed$no.ind/4),"\n"))
    ret=localdata
    tidx=(1:4)+ (i-1)*4
    tidd=tidx<=bed$no.ind
    tidx=tidx[tidd]
    if(any(!tidd)){
        ret=ret[tidd,,drop=FALSE]
    }
    attr(ret,"idx")=tidx
    ret
}

#' @title Get Bed Data As A Matrix
#' @description
#' Gets a bed file as a regular R matrix.
#' @param bed The "rbed" object to read
#' @param meanimpute (default FALSE) Whether to replace missing values with their mean (if normalise=TRUE only)
#' @param normalise (default TRUE) Whether to normalise data by SNP frequency, or return counts.
#' @param verbose (default TRUE) whether to show progress for individual access
#' @return A matrix of data with the additional attribute "idx" referring to the indexes of the returned individuals.
#' @seealso \code{\link{bedasbigmatrix}}
#' @export
bedasmatrix=function(bed,meanimpute=FALSE,normalise=TRUE,verbose=TRUE){
    ret=matrix(NA,nrow=bed$no.ind,ncol=bed$no.snp)
    colnames(ret)=paste0("SNP",1:dim(ret)[2])
    for(i in 1:ceiling(bed$no.ind/4)){
        v=get_data(bed,i,meanimpute,normalise,verbose)
        ret[attr(v,"idx"),]=v
    }
    ret
}

#' @title Get Bed Data As A BigMatrix
#' @description
#' Gets a bed file as a BigMatrix object from package "big.memory".
#' @param bed The "rbed" object to read
#' @param filename (default NULL, meaning use in-memory representation) The location of a file to use to store the matrix
#' @param meanimpute (default FALSE) Whether to replace missing values with their mean (if normalise=TRUE only)
#' @param normalise (default TRUE) Whether to normalise data by SNP frequency, or return counts.
#' @param verbose (default TRUE) whether to show progress for individual access
#' @return A bigmatrix of
#' @seealso \code{\link{bedasbigmatrix}}
#' @export
bedasbigmatrix=function(bed,filename=NULL,meanimpute=FALSE,normalise=TRUE,verbose=TRUE){
    if (!requireNamespace("bigmemory", quietly = TRUE)){
        stop("Package \"bigmemory\" needed for this function to work. Please install it with install.packages(\"bigmemory\").",
             call. = FALSE)
    }
    if(is.null(filename)){
        ret <- bigmemory::big.matrix(nrow=bed$no.ind, ncol=bed$no.snp,
                                type="double",
                          dimnames=c(NULL,NULL))
    }else{
        unlink(paste0(filename,".bk"))
        unlink(paste0(filename,".desc"))
        ret <- bigmemory::filebacked.big.matrix(nrow=bed$no.ind, ncol=bed$no.snp,
                                type="double",
                                backingfile=paste0(filename,".bk"),
                                backingpath=".",
                                descriptorfile=paste0(filename,".desc"),
                                dimnames=c(NULL,NULL))
    }
    for(i in 1:ceiling(bed$no.ind/4)){
        v=get_data(bed,i,meanimpute,verbose)
        ret[attr(v,"idx"),]=v
    }
    ret
}
