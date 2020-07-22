#' @title Write Predicted PCs In Standard Format
#' @description
#' Output the predicted PCs into a file, in the plink covariate format.
#' @param outfile file to write to
#' @param dat object of class "mergeddata" as returned by \code{\link{mergeref}}, or a dataframe with "FamID" amd "IndID" columns
#' @param pred The predicted PCs
#' @return NULL, invisibly
#' @export
writepred=function(outfile,dat,pred){
    if(is(dat,"mergeddata")) {
        indinfo=dat$indinfo
    }else if(is(dat,"mergedrbed")){
        indinfo=dat$indinfo
    }else{
        indinfo=dat
    }
    fid=as.character(indinfo[,"FamID"])
    iid=as.character(indinfo[,"IndID"])
    outdat=data.frame("#FID"=fid,
                      "IID"=iid,
                      pred)
    cat("#FID\tIID\t",file=outfile)
    cat(paste(colnames(pred),paste=("\t")),file=outfile,append=T)
    cat("\n",file=outfile,append=T)
    utils::write.table(outdat,file=outfile,quote=FALSE,
                       sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)
    invisible(NULL)
}


#' @title Test PCA predictions with access to the truth
#' @description
#' Read the true PCs from a file and compare this to our predicted PCs, for error checking.
#' @param testfile file conmtaining the truth
#' @param pred The predicted PCs
#' @param outfile optional; the file that \code{\link{writepred}} would use; the outut is written to <outfile>.comparison
#' @param verbose Whetyher output commentry on the tests
#' @return a list of the intermediate values, invisibly. You want the $fullres element which is also written to file.
#' @export
test_pca=function(testfile,pred,outfile=NULL,verbose=TRUE){
    trueres=as.matrix(utils::read.table(testfile,header=T,row.names=1,comment.char ="")[,-1])

    pred=pred[,1:min(dim(pred)[2],dim(trueres)[2])]
    trueres=trueres[,1:min(dim(pred)[2],dim(trueres)[2])]
    
    allvals=data.frame(true=as.numeric(trueres),pred=as.numeric(pred))
    tests=function(data){
        tmp=summary(stats::lm(true~pred,data=data))$coefficients[,1]
        tmp-c(0,1)
    }
    testbypc=t(sapply(1:dim(trueres)[2],function(i){
        tvals=data.frame(true=as.numeric(trueres[,i]),pred=as.numeric(pred[,i]))
        tests(tvals)
    }))
    rownames(testbypc)=colnames(trueres)
    fullres=rbind(full=tests(allvals),
                  testbypc)
    maxdev=max(abs(fullres))
    if(verbose) cat(paste("Maximum deviation is ",maxdev,"\n"))
    if(maxdev<0.01) cat(paste("TESTING CONCLUSION: Probably OK\n"))
    if(maxdev>=0.01) {
        cat(paste("TESTING CONCLUSION: PROBLEM!\n"))
        cat("TRUTH head:\n")
        print(trueres[1:5,1:5])
        cat("PREDICTED head:\n")
        print(pred[1:5,1:5])
    }
    if(!is.null(outfile)) {
        outfile2=paste0(outfile,".comparison")
        if(verbose) cat(paste("Writing regression deviations to",outfile2,"\n"))
        utils::write.table(fullres,outfile2,quote=FALSE,col.names=FALSE)
    }
    return(invisible(list(truth=trueres,
                          pred=pred,
                          maxdev=maxdev,
                          testbypc=testbypc,
                          tests=tests,
                          fullres=fullres)))
}
