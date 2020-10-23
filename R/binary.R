#' @title Convert Rbed object to binary haplotype format
#'
#' @description Convert an rbed object from plink's native file format into a haplotype-based binary matrix. This involves pseudo-phasing, where alleles are assigned at random to either of the two haplotypes.
#' @param dat an rbed object as output by \code{\link{readbed}}
#' @param ploidy the ploidy of the data. The only plausible value for plink format is 2.
#' @return a "hapmatrix" object, which is a list containing:
#' \itemize{
#' \item data: The SNP matrix in haplotype format with pN rows and L columns (p=ploidy)
#' \item snpinfo: as provided from the rbed object
#' \item indinfo: as provided from the rbed object
#' \item hapinfo: an N by p matrix giving individual identifiers for the haplotypes
#' }
#' @export
#' 
bedashap=function(dat,ploidy=2){
    ## Convert genotype data into haplotype data
    if(!is(dat,"rbed")) stop("Must provide rbed object as produced by readbed!")
    mat=bedasmatrix(dat,normalise=FALSE)
    ## We create "unphased" haplotypes by randomly assigning SNPs to haplotypes
    ret=matrix(0,dim(mat)[1]*ploidy,dim(mat)[2])
    for(i in 1:dim(mat)[1]){
        v=mat[i,]
        pr=stats::runif(dim(mat)[2],0,1)
        for(p in 1:ploidy){
            ret[(i-1)*ploidy + p,v==2]=1
            tp=(pr>((p-1)/ploidy)) & (pr<((p)/ploidy))
            ret[(i-1)*ploidy + p,(v==1)& tp]=1
        }
    }
    r=list(data=ret,
           snpinfo=dat$snpinfo,
           hapinfo=matrix(rep(dat$indinfo$IndID,each=ploidy),
                          ncol=ploidy,byrow=TRUE),
           indinfo=dat$indinfo
           )
    class(r)="hapmatrix"
    r
}

#' @title Convert binary haplotype format to hex format
#' @description Convert a matrix of 0s and 1s into a matrix of 0-f hex values, each describing 4 snps
#' @param hmat A binary matrix of size N by L
#' @param direction Either snpwise or hapwise, controls whether we place all values for a snp together (snpwise) or for a haplotype
#' @return A character matrix of size N by ceiling(L/4) of hex values
#' @export
haps2hex=function(hmat,direction="snpwise"){
    requireNamespace("BMS")
    if(direction=="hapwise"){
        if(dim(hmat)[2]/4 != floor(dim(hmat)[2]/4)){
            hmat=cbind(hmat,matrix(0,nrow=dim(hmat)[1],dim(hmat)[2]%%4))
        }
        lhex=dim(hmat)[2]/4
        r=sapply(1:lhex,function(h){
            strsplit(BMS::bin2hex(t(hmat[,4*(h-1)+1:4])),"")[[1]]
        })
    }else if(direction=="snpwise"){
        if(dim(hmat)[1]/4 != floor(dim(hmat)[1]/4)){
            hmat=rbind(hmat,matrix(0,nrow=dim(hmat)[1]%%4,dim(hmat)[2]))
        }
        lhex=dim(hmat)[1]/4
        r=sapply(1:lhex,function(h){
            strsplit(BMS::bin2hex(hmat[4*(h-1)+1:4,]),"")[[1]]
        })
    }else{
        stop("Unrecognised direction in haps2hex. Options are snpwise or hapwise.")
    }
    r
}
    
#' @title Convert binary haplotype format to hex formatted output file
#'
#' @description Represent a binary haplotype in efficient hex file format. Creates:
#' \itemize{
#' \item <fileroot>.binhaps: The matrix in efficient binary form
#' \item <fileroot>.hapinfo: A tab delimited file with N rows and p columns, giving the individual identifier for each haplotype (of which there are are Np)
#' \item <fileroot>.fam: Plinks' fam format
#' \item <fileroot>.bim: Plinks' bim format
#' }
#' 
#' @param dat an hapmatrix object as created by \code{\link{bedashap}}.
#' @param fileroot the file root of the outputs.
#' @param direction Either snpwise or hapwise, controls whether we place all values for a snp together (snpwise) or for a haplotype
#' @return NULL invisibly
#' @export
#' 
writebinhap=function(dat,fileroot,direction="snpwise"){
    ## hapwise == row major, means give each haplotype as a contiguous string
    ## snpwise == column major, means give each snp as a contiguous string
    ## writes a haplotype object
    if(!is(dat,"hapmatrix")) stop("dat must be a hapmatrix object as returned by bedashap")
    if(direction=="hapwise") {
        fn=paste0(fileroot,".binhaps")
    }else if(direction=="snpwise"){
        fn=paste0(fileroot,".binsnps")
    }else{
        stop("Invalid direction in writebinhap! options are snpwise or hapwise.")
    }
    d=paste(haps2hex(dat$data,direction),collapse="")
    writeLines(d,fn)
    
    utils::write.table(dat$indinfo,paste0(fileroot,".fam"),
                quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
    utils::write.table(dat$snpinfo,paste0(fileroot,".bim"),
                quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
    utils::write.table(dat$hapinfo,paste0(fileroot,".hapinfo"),
                quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
    return(invisible(NULL))
}

