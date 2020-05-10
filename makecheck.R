
library("devtools")
library(roxygen2)
setwd("/Users/madjl/code")

setwd("pcapred")
document()
check()


## Make distributable packages
setwd("/Users/madjl/code")
system("rm -f pcapred.tar.gz")
system("tar --exclude=.git -czvf pcapred.tar.gz pcapred")
install.packages("pcapred.tar.gz",repos = NULL, type="source")

library("pcapred")
