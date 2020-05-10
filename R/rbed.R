#' @title Read the binary PLINK format (BED, BIM, and FAM) efficiently
#'
#' @description Require the complete set of 3 files in the binary
#' PLINK format. It includes BED file, BIM file and BAM file. For more
#' information about the binary PLINK format, please check in the manual of PLINK.
#'
#' Reads the binary file into memory but obnly stores binary format; we use accessor functions to 
#' 
#' @param bed the file location for the bed file, or the file root
#' @param bim (default NULL, meaning infer from bed) the file location for the bim file, if not using the file root
#' @param fam (default NULL, meaning infer from bed) the file location for the fam file, if not using the file root
#' @param verbose (default TRUE) Whether to report progress
#'
#' @return A list containing:
#' \itemize{
#' \item snp: The SNP matrix loaded from the BED file
#' \item snp.info (if !only.snp) the snp info from the bim file
#' \item ind.info (if !only.snp) the individual info from the fam file
#' }
#' @export
readbed<-function(bed, bim=NULL, fam=NULL,verbose=TRUE) {
    if(verbose)cat(paste0("Opening file as new dataset\n"))
    
    ret = NA
    if(!endsWith(bed, ".bed")){
        if(is.null(bim)) bim=paste0(bed,".bim")
        if(is.null(fam)) fam=paste0(bed,".fam")
        bed=paste0(bed,".bed")
    }
    
    if (!file.exists(bed)){
        cat(paste0("Error: BED file doesn't exist: ",bed,"\n"))
        return(ret)
    }
    if (!file.exists(bim)){
        bim=paste0(bim,".gz")
        if (!file.exists(bim)){
            cat(paste0("Error: BIM file doesn't exist: ",bim,"\n"))
            return(ret)
        }
    }
    if (!file.exists(fam)){
        fam=paste0(fam,".gz")
        if (!file.exists(fam)){
            cat(paste0("Error: FAM file doesn't exist: ",fam,"\n"))
            return(ret)
        }
    }
                                        #Read BIM file
    if(verbose)cat(paste0("... Reading SNP info from ",bim,"\n"))
    snpinfo <- utils::read.table(bim,header=FALSE,as.is=TRUE)
    colnames(snpinfo) = c("chr","ID","GD","position","allele1","allele2")
    
                                        #Read FAM file
    if(verbose)cat(paste0("... Reading Individual info from ",fam,"\n"))
    indinfo = utils::read.table(fam,header=FALSE,as.is=TRUE)
    colnames(indinfo) = c("FamID","IndID","PatID","MatID","sex","phenotype")
    
    no.ind = dim(indinfo)[1]
    no.snp = dim(snpinfo)[1]
    
    if(verbose)cat(paste0("... Found ",no.ind," individuals and ",no.snp," SNPs\n"))
    ##Read BED file
    fh = file(bed, 'rb')
    ##Read the first three bytes to check file format
    buff = readBin(fh, what="raw",n=3)
    if (sum(buff[1:2] == c('6c','1b')) != 2){
        cat(paste0("Error: BED file is not in a correct format: ",bed,"\n"))
        return(ret)
    }
    no.byte.to.read = NA
    no.loop = NA
    if (buff[3] == '01'){
        no.byte.to.read = ceiling(no.ind/4.0)
        no.loop = no.snp
    }else{
        no.byte.to.read = ceiling(no.snp/4.0)
        no.loop = no.ind
    }
    if(verbose)cat(paste0("... Reading BED file in units of ",no.byte.to.read,"\n"))
    bindata = readBin(fh, what="raw",n=(no.byte.to.read*no.loop))
    if(verbose)cat(paste0("... Loaded raw data\n"))

    close(fh)
    bed=list(
        snpinfo=snpinfo,
        indinfo=indinfo,
        bim=bim,
        fam=fam,
        bed=bed,
        no.ind=no.ind,
        no.snp=no.snp,
        no.byte.to.read=no.byte.to.read,
        no.loop=no.loop,
        datkeep=numeric(),
        bindata=bindata
    )
    class(bed)="rbed"
    return(bed)
}

#' @title Extract a SNP from a binary format bed object
#' @description Convert the raw BED data into an integer string, apply any allele flipping inferred from the reference
#' 
#' @param bed An "rbed" object as returned from \code{\link{readbed}} or \code{\link{mergeref}}
#' @param x The index of the SNP to extract
#'
#' @return a vector of length bed$no.ind containing the SNP values for each individual
#' @seealso \code{\link{getind}}, \code{\link{getinds}}
#' @export
getsnp=function(bed,x){
    idx1=(x-1)*bed$no.byte.to.read+1
    idx2=(x)*(bed$no.byte.to.read)
    snpbits=matrix(as.integer(rawToBits(bed$bindata[idx1:idx2])),
                   ncol=2,byrow=T)[1:bed$no.ind,]    
    ret=rowSums(snpbits)
    ret[which((snpbits[,2]-snpbits[,1])<0)]=NA
    if(length(bed$flip)>0 && bed$flip[x]) ret=2-ret
    ret
}

#' @title Compute the allele frequencies for bed file
#' @description Adds the SNP frequency calculations to an "rbed" object.
#'
#' @param bed An "rbed" object as returned from \code{\link{readbed}} or \code{\link{mergeref}}
#' @param verbose (default TRUE) whether to report on activities
#'
#' @return The bed object provided with the $freq
#' @export
getfreqs=function(bed,verbose=TRUE){
    if(verbose)cat("... Calculating SNP frequencies\n")
    ret=sapply(1:bed$no.snp,function(x){
        snp=getsnp(bed,x)
        sum(snp[!is.na(snp)])/sum(!is.na(snp))/2
    })
    if(length(bed$datkeep)>0) ret=ret[bed$datkeep]
    if(length(bed$flip)>0) ret[bed$flip]=1-ret[bed$flip]
    ret
}

#' @title Extract a chunk of individuals from a binary format bed object
#' @description The binary bed format stores 4 individuals in a byte. This returns the SNP details for these 4 individuals, applying any SNP filtering for a merged dataset as well as flipping any alleles.
#' 
#' @param bed An "rbed" object as returned from \code{\link{readbed}} or \code{\link{mergeref}}
#' @param x The index of the chunk to extract (NOT an individual index)
#'
#' @return a matrix of dimension 4 by (bed$no.snp or length(bed$datkeep)) containing the SNP values for each individual
#' @seealso \code{\link{getind}}, \code{\link{getsnp}}
#' @export
getinds=function(bed,x){
    myidx=x + bed$no.byte.to.read*((0:(bed$no.loop-1)))
    snpbits=matrix(as.integer(rawToBits(bed$bindata[myidx])),
                   ncol=2,byrow=T) 
    snpinds=rowSums(snpbits)
    snpinds[which((snpbits[,2]-snpbits[,1])<0)]=NA
    mymat=matrix(snpinds,nrow=4,byrow=F)
    if(length(bed$datkeep)>0) mymat=mymat[,bed$datkeep,drop=FALSE]
    if(length(bed$flip)>0) mymat[,bed$flip]=2-mymat[,bed$flip]
    mymat
}

#' @title Extract a individual from a binary format bed object
#' @description The binary bed format stores 4 individuals in a byte. This returns the SNP details for these one individual only. It is more efficient to get the chunks with \code{\link{getinds}}.
#' 
#' @param bed An "rbed" object as returned from \code{\link{readbed}} or \code{\link{mergeref}}
#' @param x The index of the individual to access
#'
#' @return a vector of length (bed$no.snp or length(bed$datkeep)) containing the SNP values for each individual
#' @seealso \code{\link{getinds}}, \code{\link{getsnp}}
#' @export
getind=function(bed,x){
    xx=(x%%4)
    if(xx==0) xx=4
    getinds(bed,ceiling(x/4))[xx,]
}

#' @title Merge your dataset into a reference using the efficient binary "rbed" format
#' @description
#'
#' Harmonise your dataset with a reference, performing sanity checks and standardization of REF/ALT alleles.
#' 
#' @param dat An object of class "rawdata" as returned by \code{\link{readinput}}.
#' @param ref (default: the provided \code{\link[pcapred.ukbpcs18]{ukb_pcs_18}}) The reference dataset of class "referencedata" with which you want to merge your data. This is as returned by \code{link{readreference}}.
#' @param mergeon (default ID) how to merge the SNPs. "ID" is always available but assumes you have rsids for the data. Alternative is "CHRPOS" which uses the extended freq information we provide in the reference dataset.
#' @param verbose whether to say which algorithm is being used
#' @return An object of class "mergeddata" which is a list containing:
#' \itemize{
#' \item snpstats (as in ref, see \code{link{readreference}})
#' \item snpweights (as in ref, see \code{link{readreference}})
#' \item pcvals (as in ref, see \code{link{readreference}})
#' \item mode  (as in ref, see \code{link{readreference}})
#' \item transpose (as in ref, see \code{link{readreference}})
#' \item minfreq (as applied, taken from parameters or ref)
#' \item refkeep A vector of indices for the retained SNPs in the reference
#' \item datkeep A vector of indices for the retained SNPs in the target dataset
#' \item flip A boolean vector of length = length(datkeep) for whether the alleles of that SNP need flipping
#' }
#' @export
mergeref=function(dat,ref=readreference(pcapred.ukbpcs18::ukb_pcs_18()),mergeon="ID",
                     verbose=TRUE){
    dat$snpstats=ref$snpstats
    dat$snpweights=ref$snpweights
    dat$pcvals=ref$pcvals
    dat$mode=ref$mode
    dat$transpose=ref$transpose
        dat$minfreq=ref$minfreq

    if(mergeon=="ID"){
        if(verbose) cat("Trying to merge on SNP IDs.\n")
        rawsnps=as.character(dat$snpinfo[,"ID"])
        panelsnps=as.character(dat$snpstats[,"ID"])
        if(all(rawsnps%in%panelsnps==FALSE)){
            rawsnps=sapply(rawsnps,function(x)strsplit(x,"_")[[1]][1])
        }
        if(all(rawsnps%in%panelsnps==FALSE)){
            panelsnps=sapply(panelsnps,function(x)strsplit(x,"_")[[1]][1])
        }
    }else if(mergeon=="CHRPOS"){
        if(verbose) cat("Trying to merge on chromosome and position.\n")
        rawsnps=paste0(as.character(dat$snpinfo[,"chr"]),"_",as.character(dat$snpinfo[,"position"]))
        panelsnps=paste0(as.character(dat$snpstats[,1]),"_",as.character(dat$snpstats[,"position"]))
    }else stop("Error: Unreckognised \"mergeon\" option; choices are \"ID\" or \"CHRPOS\".")
    
    if(all(rawsnps%in%panelsnps==FALSE)){
        stop("ERROR: Unable to find your SNPs in the reference data! Try a different \"mergeon\" option. Do they have the same labels? Are they in the same build? The column names should be in the ID column of the .afreq.gz file.")
    }
    foundtab=sapply(c(TRUE,FALSE),function(x)sum(rawsnps%in%panelsnps==x))
    foundtab2=sapply(c(TRUE,FALSE),function(x)sum(panelsnps%in%rawsnps==x))
    if(verbose) {
        cat(paste("Compatability check:
... Data contains",foundtab[1],"SNPs in common with the reference,
... Data contains",foundtab[2],"SNPs not in the reference.
... Reference contains",foundtab2[2],"SNPs not in the data.\n"))
    }
    dat$refkeep=which(panelsnps%in%rawsnps)
    dat$datkeep=which(rawsnps%in%panelsnps)
    rawsnpskeep=rawsnps[dat$datkeep]
    
    if(verbose)cat(paste("... Retaining",length(dat$datkeep),"SNPs in the reference\n"))
    
    ## Standardize SNPs against the reference
    if(verbose) cat(paste("Standardizing SNPs against the reference\n"))
    rawrefalleles=dat$snpinfo[dat$datkeep,"allele1"]

    tref=rawrefalleles==dat$snpstats[dat$refkeep,"REF"]
    tref2=rawrefalleles==dat$snpstats[dat$refkeep,"ALT"]
    if(verbose) {
        cat(paste("...",sum(tref),"Alleles match the REF allele\n"))
        cat(paste("...",sum(tref2),"Alleles match the ALT allele\n"))
    }

    dat$flip=tref2[tref | tref2]
    
    ## Detect problems
    tref3=!(tref | tref2)
    if(sum(tref3)>0){
        if(verbose) cat(paste("...",sum(tref3),"Alleles match no allele and are omitted\n"))
        dat$refkeep=dat$refkeep[!tref3]
        dat$datkeep=dat$datkeep[!tref3]
        if(mean(tref3)>0.1) cat(paste("WARNING! MANY INCONSISTENCIES BETWEEN PANEL AND DATA!\n"))
    }else if(verbose) {
        cat(paste("...",sum(tref3),"Alleles match no allele so none are omitted\n"))
    }
    class(dat)=c("rbed","mergedrbed")
    dat
}

#' @title Create predictions for eigenvalues from binary stored data
#' @description
#' Make a genome-wide prediction based on incomplete data
#' @param dat An object of class "mergedrbed" as returned by \code{\link{mergeref}}
#' @param minfreq The minimum reference frequency to allow, to prevent rare SNPs having a disproportionate effect
#' @param verbose whether to say which algorithm is being used
#' @return A matrix of size N by K of predicted Eigenvectors
#' @export
predictpcs=function(dat,minfreq=NULL,verbose=TRUE){
    if(verbose) cat(paste("Predicting eigenvectors\n"))
    if(is.null(minfreq)) minfreq=dat$minfreq
    if(all(is.null(dat$weights))) {
        ## Scale the weights into canonical form
        if(verbose) cat(paste("... Creating normalised SNP weights\n"))
        dat$weights=transform_pcs(dat$snpweights[,dat$refkeep],dat$pcvals,mode=dat$mode)
    }

    nloop=ceiling(dat$no.ind/4)
    reported=c(FALSE,FALSE)
    pred=matrix(NA,nrow=nloop*4,ncol=dim(dat$weights)[2])
    colnames(pred)=rownames(dat$snpweights)
    for(i in 1:nloop){
        ndat=normalise(getinds(dat,i),dat$snpstats[dat$refkeep,"ALT_FREQS"],minfreq=minfreq)
        natab=sapply(c(TRUE,FALSE),function(x)sum(is.na(as.numeric(ndat))==x))
        if(natab[1]==0){
            if(verbose && (!reported[1])) {
                cat("Sections without missing data detected, using fast algorithm where possible.\n")
                reported[1]<-TRUE
            }
            ret = pred_simple(ndat,dat$weights,L=dim(ndat)[2],dat$mode)
        }else{
            if(verbose&& (!reported[2])) {
                cat("Sections with missing data detected, using missing-data aware algorithm where necessary.\n")
                reported[2]<-TRUE
            }
            ret = pred_missing(ndat,dat$weights,dat$mode)
        }
        pred[(1:4+(i-1)*4),]=ret
    }
    return(pred=pred[1:dat$no.ind,,drop=FALSE])
}
