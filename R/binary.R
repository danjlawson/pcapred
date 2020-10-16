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
    if(!is(mat,"rbed")) stop("Must provide rbed object as produced by readbed!")
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

haps2hex=function(hmat){
    requireNamespace("BMS")
    if(dim(hmat)[2]/4 != floor(dim(hmat)[2]/4)){
        hmat=cbind(hmat,matrix(0,nrow=dim(hmat)[1],dim(hmat)[2]%%4))
    }
    lhex=dim(hmat)[2]/4
    r=sapply(1:lhex,function(h){
        strsplit(BMS::bin2hex(t(hmat[,4*(h-1)+1:4])),"")[[1]]
    })
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
#' @param major (default="hap") whether to flatten the matrix by reporting each haplotype in sequence, or (if major="snp") each set of 4 snps in sequence.
#' @return NULL invisibly
#' @export
#' 
writebinhap=function(dat,fileroot,major="hap"){
    ## hap == row major, means give each haplotype as a contiguous string
    ## snp == column major, means give each set of 4-snps as a contiguous string
    ## writes a haplotype object
    if(!is(dat,"hapmatrix")) stop("dat must be a hapmatrix object as returned by bedashap")
    if(major=="hap"){
        d=paste(t(haps2hex(dat$data)),collapse="")
        writeLines(d,paste0(fileroot,".binhaps"))
    }else if(major=="snp"){
        d=paste(haps2hex(dat$data),collapse="")
        writeLines(d,paste0(fileroot,".tbinhaps"))
    }else{
        stop("Invalid major; valid options are \"hap\" or \"snp\"")
    }
    utils::write.table(dat$indinfo,paste0(fileroot,".fam"),
                quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
    utils::write.table(dat$snpinfo,paste0(fileroot,".bim"),
                quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
    utils::write.table(dat$hapinfo,paste0(fileroot,".hapinfo"),
                quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
    return(invisible(NULL))
}

