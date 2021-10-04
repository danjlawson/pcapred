# R package *pcapred*

This R package allows you to predict PCA values in the UK Biobank reference data (or your own source of data). It includes key reference data derived from the UK Biobank.

If you are not an R user, or want to use a command line interface, one is provided at [pcapred-script](https://github.com/danjlawson/pcapred-script).

You may want [pcapred.largedata with 200 PCs](https://github.com/danjlawson/pcapred.largedata) in the same format as the 18 provided with this package (via the included [pcapred.ref](https://github.com/danjlawson/pcapred.largedata) package).

## Getting started

Install the package using:

```{r}
remotes::install_github("danjlawson/pcapred")
```

## A complete example

We have provided a sample dataset from 1000 Genomes in bim/bed/fam format.

```{r}
library("pcapred")
mydata=pcapred.ref::onek_genomes_tiny() # Gets the file location of the tiny bim/bed/fam data included in pcapred's helper data package, pcapred.ref
dat=readbed(mydata) # Read "your" data
dat=mergeref(dat)     # Merge with the reference (using the included standard reference of 18 UK Biobank Pcs by default)
pred=predictpcs(dat)  # Predict the first 18 UK Biobank PCs
writepred("projected.eigenvals",dat,pred) # Write output in plink --covar format
```

## Modification for your data:

You just change the file location to point to your own data.

```{r}
library("pcapred")
dat=readbed("path/to/plink/bimbedfamroot") # omit the .bed file ending
dat=mergeref(dat) # Merge your dataset with the provided reference
pred=predictpcs(dat) # Predict the first 18 UK Biobank PCs
writepred("/path/to/projected.eigenvals",dat,pred) # Write the PCs in the correct format for plink's --covar command.
```

## A real example

If you want to get the 1000 Genomes data, you can analyse it similarly. This is how you get the data and **prune it to only have SNPs that are in the reference**. The [cog-genomics](https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3) site has up-to-date links, though the  [PlinkQC package vignette](https://cran.r-project.org/web/packages/plinkQC/vignettes/Genomes1000.pdf) is also helpful.

To get the 1000 Genomes SNPs in the correct format:
```{r}
download.file("https://www.dropbox.com/s/afvvf1e15gqzsqo/all_phase3.pgen.zst?dl=1",
              destfile = "all_phase3.pgen.zst")
download.file("https://www.dropbox.com/s/op9osq6luy3pjg8/all_phase3.pvar.zst?dl=1",
              destfile = "all_phase3.pvar.zst")
download.file("https://www.dropbox.com/s/yozrzsdrwqej63q/phase3_corrected.psam?dl=1",
destfile = "all_phase3.psam")
system(paste0("zless ", pcapred.ref::ukb_pcs_18(),".load.gz | cut -f1-6 > ref.bim"))
system("plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen")
system("plink2 --pfile all_phase3 vzs --extract ref.bim --make-bed --out 1000G_forPCA")
unlink(c("all_phase3.pgen.zst","all_phase3.psam","all_phase3.pvar.zst","all_phase3.pgen","ref.bim"))
```

Now you have nice small datasets in the file "1000G_forPCA".

```{r}
library("pcapred")
t1=as.numeric(Sys.time())
mydata="1000G_forPCA"
dat=readbed(mydata)
t2=as. numeric(Sys.time())
dat=mergeref(dat)     # Merge with the reference (using the included standard reference of 18 UK Biobank Pcs by default)
t3=as.numeric(Sys.time())
pred=predictpcs(dat)  # Predict the first 18 UK Biobank PCs
t4=as.numeric(Sys.time())
writepred("1000G_forPCA.projected.eigenvals",dat,pred) # Write output in plink --covar format
print(paste("Total time: ",round(t4-t1)))
print(paste("Data reading time: ",round(t2-t1)))
print(paste("Merging time: ",round(t3-t2)))
print(paste("Prediction time: ",round(t4-t3)))
```

On a good compute node, this takes 57 seconds total. There are 2504 individuals. The code is not parallelized but uses memory efficiently and therefore scales to UK Biobank. The cost is linear in the number of SNPs provided in the data and the number of individuals.

The missing-data aware algorithm is automatically used where necessary. It is a little slower by a factor of around 2.

## Limitations

There is no parallelisation or memory management. It takes around 6-7 hours to process all 500K UK Biobank participants. You will need about 32Gb of memory.

Merging is currently done best on IDs.

## Additional information

The UK Biobank SNPs chosen for PC analysis may not all have been available in your dataset. We scale the remaining SNPs assuming that they are missing at random to generate predictions that should be similar if you do not have all of the SNPs. This should work fine if you use 1000 Genomes imputation.

## Licence Information

Author: Daniel Lawson (dan.lawson@bristol.ac.uk)
Date: 2020
Licence: This code is released under the GPLv3.
