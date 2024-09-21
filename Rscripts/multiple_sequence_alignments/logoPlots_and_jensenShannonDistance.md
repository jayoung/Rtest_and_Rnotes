logoPlots_and_jensenShannonDistance
================
Janet Young

2024-09-20

Goals:  
- show how to make a logo plot  
- show how to calculate Jensen-Shannon divergence between two (aligned)
alignments at each position

JSD is basically working - xxx tidy it up a lot and use tidyverse code

xxx add logo plot

# Load libraries

``` r
knitr::opts_chunk$set(echo = TRUE)
# library(ggtree)
# library(ape)
# library(here)
library(Biostrings)
library(entropy)
```

# Define some functions

We define a function that matches Mike Doud’s function to quantify
measure of AA preference difference using the [Jensen-Shannon
distance](https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence)

``` r
# myEnt is a function to get Shannon entropy. 
# x is a vector of numbers that sums to 1 representing amino acid frequencies at a given position in a single alignment. In my case it has length 21 (includes gap character, -)
myEnt <- function (x) {  entropy.empirical(x, "log2") }

# myJSD is a function to get JSD comparing two frequency profiles (prefs1 and prefs2)
myJSD <- function( prefs1, prefs2 ) {
    firstTerm <- myEnt( (prefs1 + prefs2)/2 )
    secondTerm <- (myEnt(prefs1) + myEnt(prefs2))/2
    myResult <- sqrt( firstTerm - secondTerm )
    return(myResult)
}

# myJSD_3classes is a function to get JSD comparing three frequency profiles (prefs1, prefs2 and prefs3). The upper bound of the three-class case is sqrt(log2(3)) = 1.258953
#  (that is correct according the wikipedia page)
myJSD_3classes <- function( prefs1, prefs2, prefs3 ) {
    firstTerm <- myEnt( (prefs1 + prefs2 + prefs3)/3 )
    secondTerm <- (myEnt(prefs1) + myEnt(prefs2) + myEnt(prefs3))/3
    myResult <- sqrt( firstTerm - secondTerm )
    return(myResult)
}

# myJSD_4classes is a function to get JSD comparing three frequency profiles (prefs1, prefs2, prefs3 and prefs4). The upper bound of the four-class case is sqrt(log2(4)) = 1.414214
myJSD_4classes <- function( prefs1, prefs2, prefs3, prefs4 ) {
    firstTerm <- myEnt( (prefs1 + prefs2 + prefs3 + prefs4)/4 )
    secondTerm <- (myEnt(prefs1) + myEnt(prefs2) + myEnt(prefs3) + myEnt(prefs4))/4
    myResult <- sqrt( firstTerm - secondTerm )
    return(myResult)
}
```

# read in alignments

# Read in example data: short H2A alignment (Antoine)

In this case I have a single alignment of the histone fold domain of
various H2A family members. It contains canonical H2A, H2A.B, H2A.L,
H2A.P (8 species each) and the marsupial-specific H2A.R (5 sequences)
(total of 37 sequences).

I read in the alignment, figure out which H2A variant each sequence is
from, and split the single master alignment into 5 individual
alignments, one for each H2A variant.

``` r
aln_file <- "exampleProtAln_shortH2As_histoneFoldDomain.fa"

masterAln <- readAAStringSet(aln_file)
# simplify the sequence names by removing the description
names(masterAln) <- sapply(strsplit(names(masterAln), " "), "[[", 1)
```

Then I figure out which variant each sequence is from and split the
alignment into one smaller alignment per variant.

``` r
# figure out which variant each sequence is from
masterAlnSeqTypes <- sapply(strsplit(names(masterAln), "_"), "[[", 1)
masterAlnSeqTypes <- gsub("R[12]$","R", masterAlnSeqTypes)
# table(masterAlnSeqTypes)

# split the alignment into those categories. We get a list object containing all 5 alignments
masterAlnSplit <- split(masterAln, masterAlnSeqTypes)
```

xxxx below here can be tidied

``` r
### get count matrices in a way that includes every amino acid and the gap character:
myAAtoTabulate <- c("-", AA_STANDARD)
masterAlnSplit_counts <- lapply( masterAlnSplit, function(x) {
    y <- as.data.frame(consensusMatrix(x))
    missingRows <- setdiff (myAAtoTabulate, rownames(y))
    for (thisRow in missingRows) {
        y[thisRow,] <- rep(0, dim(y)[2])
    }
    y <- y[ myAAtoTabulate, ]
    y
} )

masterAlnSplit_frequencies <- lapply(masterAlnSplit_counts, function(x) {
    y <- x / colSums(x)
    y
})
```

``` r
myJSD (masterAlnSplit_frequencies[["H2A"]][,7], masterAlnSplit_frequencies[["H2A.B"]][,7] )
```

    ## [1] 0.8467096

``` r
###
myDistances <- data.frame(row.names= paste(1:82, strsplit( as.character(masterAlnSplit[["H2A"]][[1]] ), "")[[1]], sep=""), pos=1:82 )
myDistances[,"BvsH2A"] <- sapply (1:dim(masterAlnSplit_frequencies[["H2A"]])[2], function (x) {
    myJSD(masterAlnSplit_frequencies[["H2A"]][,x], masterAlnSplit_frequencies[["H2A.B"]][,x] )
} )
myDistances[,"LvsH2A"] <- sapply (1:dim(masterAlnSplit_frequencies[["H2A"]])[2], function (x) {
    myJSD(masterAlnSplit_frequencies[["H2A"]][,x], masterAlnSplit_frequencies[["H2A.L"]][,x] )
} )
myDistances[,"PvsH2A"] <- sapply (1:dim(masterAlnSplit_frequencies[["H2A"]])[2], function (x) {
    myJSD(masterAlnSplit_frequencies[["H2A"]][,x], masterAlnSplit_frequencies[["H2A.P"]][,x] )
} )

myDistances[,"RvsH2A"] <- sapply (1:dim(masterAlnSplit_frequencies[["H2A"]])[2], function (x) {
    myJSD(masterAlnSplit_frequencies[["H2A"]][,x], masterAlnSplit_frequencies[["H2A.R"]][,x] )
} )
```

``` r
plotXaxisFunction <- function() {
    axis(1, at=(myDistances[,"pos"]-0.5), labels=rownames(myDistances), las=2, cex.axis=0.5, tick=FALSE, hadj=0.4 )
}
plotYaxisFunction <- function(thisYmin=0) {
    myYmgp <- c(1.4, 0.8, 0)
    axis(2, at=c(thisYmin,1), labels=c(thisYmin,1), las=2, mgp=myYmgp)
    title(ylab="JSD", mgp=myYmgp)
}
myBarplot <- function(thisColName, thisTitle, thisYmin=0, ...) {
    barplot( myDistances[,thisColName], main=thisTitle, space=0, xaxs="i", yaxt="n", ...)
    plotXaxisFunction()
    plotYaxisFunction(thisYmin)
}
myBarplot( thisColName="BvsH2A", thisTitle="H2A.B vs H2A")
```

![](logoPlots_and_jensenShannonDistance_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

# Finished

show R version used, and package versions

``` r
sessionInfo()
```

    ## R version 4.4.0 (2024-04-24)
    ## Platform: x86_64-apple-darwin20
    ## Running under: macOS Ventura 13.6.9
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] entropy_1.3.1       Biostrings_2.72.1   GenomeInfoDb_1.40.1
    ## [4] XVector_0.44.0      IRanges_2.38.0      S4Vectors_0.42.0   
    ## [7] BiocGenerics_0.50.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] digest_0.6.36           R6_2.5.1                zlibbioc_1.50.0        
    ##  [4] fastmap_1.2.0           xfun_0.45               GenomeInfoDbData_1.2.12
    ##  [7] knitr_1.47              UCSC.utils_1.0.0        htmltools_0.5.8.1      
    ## [10] rmarkdown_2.27          cli_3.6.3               compiler_4.4.0         
    ## [13] highr_0.11              httr_1.4.7              rstudioapi_0.16.0      
    ## [16] tools_4.4.0             evaluate_0.24.0         yaml_2.3.8             
    ## [19] crayon_1.5.3            jsonlite_1.8.8          rlang_1.1.4
