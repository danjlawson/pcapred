

#' @title Read Snp Weights/Loadings from a file
#' @description
#' Read SNP weights from a variety of formats (plink/flashpca, either regular or transposed)
#' @param snpweightsfile The file to read
#' @param mode either flashpca or plink; default "flashpca"
#' @param transpose default FALSE, meaning treat columns as PCs and rows as SNPs. if TRUE, assume the laoding file is transposed from this.
#' @param colnames (default NULL, meaning do no naming) names for the SNPs' comes from 
#' @param verbose (default TRuE) whether to show progress for file reading
#' @return a matrix of SNP weights with K rows and L columns
#' @export
readsnpweights <- function(snpweightsfile,mode="flashpca",
                           transpose=FALSE,colnames=NULL,verbose=TRUE){
    if(mode=="plink"){
        if(transpose){
            snpweightsraw=data.table::fread(snpweightsfile,
                                header=F,skip=4,showProgress=verbose)
            snpweights=as.matrix(snpweightsraw[,-1])
            rownames(snpweights)=as.character(as.matrix(snpweightsraw[,1]))
        }else{
            snpweightsraw=data.table::fread(snpweightsfile,
                                header=T,showProgress=verbose)
            snpweights=t(as.matrix(snpweightsraw[,-(1:4)]))
        }
    }else if(mode=="flashpca"){
        if(transpose){
            snpweightsraw=data.table::fread(snpweightsfile,
                            header=F,showProgress=verbose,skip=2)
#colnames(snpweightsraw)=as.character(read.table(snpweightsfile,nrows=1,as.is=T))
            snpweights=as.matrix(snpweightsraw[,-1])
        rownames(snpweights)=as.character(as.matrix(snpweightsraw[,1]))
        }else{ ### XXX UNTESTED CODE?
            snpweightsraw=data.table::fread(snpweightsfile,
                                header=T,showProgress=verbose)
            snpweights=t(as.matrix(snpweightsraw[,-(1:2)]))
        }
    }
    if(any(!is.null(snpweights))) colnames(snpweights)=colnames
    return(snpweights)
}

#' @title Read the target dataset
#' @description
#' Read a bim/bed
#' @param pcvalsfile file to read in
#' @return A matrix of the PCs
#' @param ... extra parameters to read.table
#' @export
readpcs <-function(pcvalsfile,...){
    as.matrix(utils::read.table(pcvalsfile,header=F,...))[,1]
}

#' @title Read the snp stats
#' @description
#' Read a afreq file from plink
#' @param afreqfile file to read in
#' @param ... extra parameters to read.table
#' @return A matrix of the SNP stats
#' @export
readsnpstats<-function(afreqfile,...){
    utils::read.table(afreqfile,header=T,comment.char ="",as.is=TRUE,...)
}


#' @title Read the Reference Data for a PC projection
#' @description
#' Read in SNP weightings, eigenvalues and SNP info.
#'
#' 
#' SNP info must be in <fileroot>.afreq.gz
#' Weights are read according to mode and transpose flag.
#'
#' Plink weights are in <fileroot>.eigenvec.var.gz/<fileroot>.eigenvec.tvar.gz
#' flashpca weights are in <fileroot>.tload.gz/<fileroot>.tload.gz
#'
#' Plink (Eigen) values are in <fileroot>.eigenval
#' flashpca (Singular) values are in <fileroot>.val
#' 
#' @param fileroot fileset to read in
#' @param mode Data mode; one of flashpca (default) or plink
#' @param transpose Whether SNP weights are transposed (default FALSE)
#' @param checkonly default FALSE. If TRUE, no file reading is done but the files are checked for existing
#' @param minfreq (default 0.001) minimum allele frequency allowed in the normalisation phase, to prevent ultra-rare SNPs from being upweighted too much
#' @param verbose default TRUE. Whether to output progress.
#' @return An object of class "referencedata" which is a list containing:
#' \itemize{
#' \item snpweights (if !checkonly) The matrix of snp weights, with SNPs as rows and PCs as columns
#' \item snpstats (if !checkonly) The matrix of snp statistics, with SNPs as rows
#' \item pcvals (if !checkonly) A vector of PC Eigenvalues (plink) or SVD Singular vectors (flashpca) as appropriate.
#' \item mode As provided in the call.
#' \item transpose As provided in the call
#' \item minfreq As provided in the call
#' \item pcvalsfile The file intended to be used for pcvals
#' \item snpweightsfile The file intended to be used for snpweights
#' \item snpstatsfile The file intended to be used for snpstats
#' }
#' @export
#' @examples
#' \dontrun{
#' tfile=pcapred.ref::ukb_pcs_18()
#' ref=readreference(tfile)
#' }
readreference <- function(fileroot,mode="flashpca",
                          transpose=FALSE,checkonly=FALSE,minfreq=0.001,verbose=TRUE) {
    ## Figure out the names of everything
    afreqfile=paste0(fileroot,".afreq.gz")
    if(mode=="plink"){
        pcvalsfile=paste0(fileroot,".eigenval")
        if(transpose){
            snpweightsfile=paste0(fileroot,".eigenvec.tvar.gz")
        }else snpweightsfile=paste0(fileroot,".eigenvec.var.gz")
    }else if(mode=="flashpca"){
        pcvalsfile=paste0(fileroot,".val")
        if(transpose){
            snpweightsfile=paste0(fileroot,".tload.gz")
        }else snpweightsfile=paste0(fileroot,".load.gz")
    }else{
        stop(paste("Unsupported mode", mode))
    }
    ## Check
    testfiles(c(snpweightsfile,pcvalsfile,snpweightsfile),verbose)
    if(checkonly) return(invisible(list(
                      mode=mode,
                      transpose=transpose,
                      minfreq=minfreq,
                      pcvalsfile=pcvalsfile,
                      snpweightsfile=snpweightsfile,
                      afreqfile=afreqfile)))
    
    if(verbose) cat(paste("Reading Refernce data from root ",fileroot,"\n"))

    ## Read the eigenvalues
    if(verbose) cat(paste("... Reading Eigenvalues from ",pcvalsfile,"\n"))
    pcvals=readpcs(pcvalsfile)

    ## Read the SNP stats
    if(verbose) cat(paste("... Reading SNP info from ",afreqfile,"\n"))
    snpstats=readsnpstats(afreqfile)

    ## Read the SNP weights
    if(verbose) cat(paste("... Reading SNP weights from ",snpweightsfile,"\n"))
    snpweights=readsnpweights(snpweightsfile,
                              mode,transpose,
                              colnames=as.character(snpstats[,"ID"]),verbose)

    ## Sanity checking
    if(dim(snpstats)[1] != dim(snpweights)[2]) {
        stop("Error: SNP stats don't describe the SNP weights file.")
    }
    if(length(pcvals) != dim(snpweights)[1]) {
        stop("Error: SNP stats don't describe the SNP weights file.")
    }

    if(verbose){
        cat(paste("... Found",dim(snpstats)[1],"SNPs in the reference\n"))
        cat(paste("... Found",length(pcvals),"eigenvalues for processing.\n"))
    }
    snpstats=readsnpstats(afreqfile)
    ref=list(snpweights=snpweights,
             snpstats=snpstats,
             pcvals=pcvals,
             mode=mode,
             minfreq=minfreq,
             transpose=transpose,
             pcvalsfile=pcvalsfile,
             snpweightsfile=snpweightsfile,
             afreqfile=afreqfile
             )
    class(ref)="referencedata"
    ref
}

#' @title Read the target dataset
#' @description
#' Read a bim/bed
#' @param inputfile file root to read
#' @param raw if TRUE, treat the input as a text file as output by plink --convert 12. Otherwise (default) read the .bim/.bed/.fam
#' @param verbose whether to report file reading progress
#' @return See readbed
#' @export
readinput <- function(inputfile,raw=FALSE,verbose=TRUE){
    if(raw){
        rawdat=data.table::fread(inputfile,showProgress=verbose)
        rawdatmat=as.matrix(rawdat[,7:dim(rawdat)[2]])
        dat=list(snp=rawdatmat,
                 indinfo=matrix(colnames(rawdat),ncol=1),
                 snpinfo=rawdat[,1:7])
    }else{
        if(verbose) cat(paste("Reading bed file:",inputfile,".\n"))
        dat=readbed(inputfile)
        dat$snp=2-dat$snp
        colnames(dat$snp)=paste0(as.character(dat$snpinfo[,"ID"]),
                                   "_",
                                   as.character(dat$snpinfo[,"allele1"]))
    }
    if(verbose) cat(paste("Found",dim(dat$snp)[1],"individuals and",dim(dat$snp)[2],"SNPS.\n"))
    class(dat)="rawdata"
    dat
}

testfiles <- function(testfiles,verbose=FALSE){
    if(verbose) cat("Testing files exist ...\n")
    for(f in testfiles){
        if(file.exists(f)){
            if(verbose) cat(paste("...",f,"exists.\n"))
        }else{
            stop(paste("Required input file",f,"not found!"))
        }
    }
    return(invisible(TRUE))
}

## Read the raw data
