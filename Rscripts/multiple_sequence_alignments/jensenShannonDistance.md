Jensen-Shannon distance calculations
================
Janet Young

2024-09-23

Goal - show how to calculate Jensen-Shannon divergence between two
(aligned) alignments at each position

# Load libraries

``` r
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(Biostrings)
library(entropy)
```

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
# class(masterAln)
# [1] "AAStringSet"
# attr(,"package")
# [1] "Biostrings"

masterAln
```

    ## AAStringSet object of length 37:
    ##      width seq                                              names               
    ##  [1]    82 SRSSRAGLQFPVGRVHRLLRKGN...TRIIPRHLQLAIRNDEELNKLL H2A_panda
    ##  [2]    82 TRSSRAGLQFPVGRVHRLLRKGN...TRIIPRHLQLAIRNDEELNKLL H2A_armadillo
    ##  [3]    82 TRSSRAGLQFPVGRVHRLLRKGN...TRIIPRHLQLAIRNDEELNKLL H2A_chineseHamster
    ##  [4]    82 SRSSRAGLQFPVGRVHRLLRKGN...TRIIPRHLQLAIRNDEELNKLL H2A_leopard
    ##  [5]    82 TRSSRAGLQFPVGRVHRLLRKGN...TRIIPRHLQLAIRNDEELNKLL H2A_mouse
    ##  ...   ... ...
    ## [33]    82 THLTTTEPQVPVSFVDHLLQEDQ...MQMTPQDVERAVDSNAEPHRQV H2A.P_pig
    ## [34]    82 SHLIRSELQCPLSYVDRLLLEDQ...MHTVPQD-DRGVGSNGQRPQNL H2A.P_leopard
    ## [35]    82 AHLITTELQVPVSYVDRLLQENQ...MPTAPQDVERAVDSSGEPYHRS H2A.P_panda
    ## [36]    82 ACLPTAELQFPVSYLDRLLQKDE...SSSVAQDVEGGVNNNREPQRQV H2A.P_rhino
    ## [37]    82 SLSARTEMEFSPSGLERLLQEDR...SHIAPLDVERGVRNNRLLRHLL H2A.P_armadillo

Then for each sequence, I figure out which variant it’s from, and I
split the alignment into one for each variant.

``` r
## figure out which variant each sequence is from
masterAlnSeqTypes <- sapply(strsplit(names(masterAln), "_"), "[[", 1)
masterAlnSeqTypes <- gsub("R[12]$","R", masterAlnSeqTypes)
# table(masterAlnSeqTypes)

## split the alignment into those categories. We get a list object containing all 5 alignments
masterAlnSplit <- split(masterAln, masterAlnSeqTypes)

# names(masterAlnSplit)
# [1] "H2A"   "H2A.B" "H2A.L" "H2A.P" "H2A.R"

## I add a sixth alignment, which is the combined B, L and P alignment
masterAlnSplit[["BandLandP"]] <- c(masterAlnSplit[["H2A.B"]], masterAlnSplit[["H2A.L"]], masterAlnSplit[["H2A.P"]])
```

# Define some useful functions

We define a function called getJSD to get [Jensen-Shannon
distance](https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence),
based on info I got from Mike Doud (Bloom lab), now published
[here](https://academic.oup.com/mbe/article/32/11/2944/982113#74376742)
(see equation 2).

Jensen-Shannon distance quantifies amino acid preference difference
between two frequency profiles. Note that JS **distance** is the square
root of JS **divergence**.

``` r
# myEnt is a function to get Shannon entropy from a single vector of frequencies.
# It's a simple wrapper function for the entropy.empirical function from the entropy package
# freqs is a vector of numbers that sums to 1 representing amino acid frequencies at a given position in a single alignment. In my case it has length 21 (includes gap character, -)
myEnt <- function (freqs) {  
    # I don't know why, but I had to put in this rounding stage before checking freqs sum to 1, otherwise some positions gave me an error
    mySum <- round(sum(freqs), digits=8)
    if(mySum != 1) {
        # cat("sum ",mySum," freqs ", freqs, "\n\n") # for troubleshooting
        stop("\n\nERROR in myEnt function - frequencies should sum to 1, but they don't\n\n")
    }
    entropy.empirical(freqs, "log2") 
}

# getJSD is a function to get JSD comparing two frequency profiles (freqs1 and freqs2). 
getJSD <- function( freqs1, freqs2 ) {
    firstTerm <- myEnt( (freqs1 + freqs2)/2 )
    secondTerm <- (myEnt(freqs1) + myEnt(freqs2))/2
    myResult <- sqrt( firstTerm - secondTerm )
    return(myResult)
}

# The upper bound of two-class JSD is sqrt(log2(2)) = 1
# The upper bound of three-class JSD is sqrt(log2(3)) = 1.258953
#  (that is correct according the wikipedia page)
# The upper bound of the four-class JSD is sqrt(log2(4)) = 1.414214
# I noticed that Wikipedia definition's of two-sample JSD does NOT divide the first term by 2, but Mike Doud's equation 2 does. I'm going to trust Mike Doud's for now.

## getJSD_multiClass is a function that works on any number of frequency profiles
getJSD_multiClass <- function(freqs_list) {
    firstTerm <- Reduce("+", freqs_list) / length(freqs_list)
    firstTerm <- myEnt(firstTerm)
    secondTerm <- sapply(freqs_list, myEnt)
    secondTerm <- sum(secondTerm) / length(freqs_list)
    myResult <- sqrt( firstTerm - secondTerm )
    return(myResult)
}

#### tests
## getJSD and getJSD_multiClass give same output on the two-class case
# getJSD(masterAlnSplit_freqs[["H2A"]][,7], 
#        masterAlnSplit_freqs[["H2A.B"]][,7]) 
# getJSD_multiClass(list(masterAlnSplit_freqs[["H2A"]][,7], 
#                        masterAlnSplit_freqs[["H2A.B"]][,7]) )
## getJSD_3classes (now removed) and getJSD_multiClass give same output on the three-class case
# getJSD_multiClass(list(masterAlnSplit_freqs[["H2A"]][,7], 
#                        masterAlnSplit_freqs[["H2A.B"]][,7],
#                        masterAlnSplit_freqs[["H2A.P"]][,7]) )
# getJSD_3classes(masterAlnSplit_freqs[["H2A"]][,7], 
#                 masterAlnSplit_freqs[["H2A.B"]][,7],
#                 masterAlnSplit_freqs[["H2A.P"]][,7]) 
```

Get amino acid (plus gap) count and frequency matrices for each
alignment

First we define a function called `getAlnCounts()` that takes an
alignment and returns a count or frequency matrix

``` r
# a tiny function that makes sure all seqs in an alignment are the same length as each other
checkAlnLengths <- function(aln) {
    if(length(unique(width(aln))) != 1) {
        stop("\n\nERROR - you supplied a ragged alignment (seqs not all the same length)\n\n")
    } else {
        return(TRUE)
    } 
}


## define the letters we want to count
# AA_STANDARD is defined in the Biostrings package and includes the usual 20 amino acids. I want to add the gap character ("-")
myAAtoTabulate <- c("-", AA_STANDARD)

## define the function
getAlnCounts <- function(aln, letters=myAAtoTabulate, as.prob=FALSE) {
    # check for ragged alns (seqs not all the same length)
    checkAlnLengths(aln)
    
    # get counts
    countsEachSeq <- lapply(1:length(aln), function(i) {
        letterFrequencyInSlidingView(aln[[i]], view.width = 1, letters=letters)
    })
    
    # if there were letters in the alignment that are not accounted for in the letters argument, the totals won't be correct.
    expectedTotals <- width(aln)[1]
    totalCountsEachSeq <- sapply(countsEachSeq, sum)
    if ( sum(totalCountsEachSeq != expectedTotals) > 0) {
        stop("\n\nERROR - the total counts didn't add up correctly. Are there letters in the alignment that are not present in the letters argument you supplied?\n\n")
    }
    
    # get total counts by position - the Reduce function takes a list object and uses the specified function on all the elements
    countTotals <- Reduce("+", countsEachSeq)
    
    # transpose so columns are positions and rows are each letter type
    countTotals <- t(countTotals)
    
    # perhaps get frequencies not counts
    if(as.prob) {
        freqs <- countTotals / colSums(countTotals)
        return(freqs)
    } else {
        return(countTotals)
    }
}
```

Now use that function on each alignment:

``` r
masterAlnSplit_counts <- lapply( masterAlnSplit, getAlnCounts, letters=myAAtoTabulate)
masterAlnSplit_freqs <- lapply( masterAlnSplit, getAlnCounts, letters=myAAtoTabulate, as.prob=TRUE)
```

Now test out the `getJSD()` function that gets Jensen-Shannon distance.

Here are frequencies for the 7th position in the H2A alignment:

``` r
masterAlnSplit_freqs[["H2A"]][,7]
```

    ## - A R N D C Q E G H I L K M F P S T W Y V 
    ## 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0

and in the H2A.B alignment

``` r
masterAlnSplit_freqs[["H2A.B"]][,7]
```

    ##     -     A     R     N     D     C     Q     E     G     H     I     L     K 
    ## 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.875 0.125 0.000 0.000 0.000 0.000 
    ##     M     F     P     S     T     W     Y     V 
    ## 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000

here’s the corresponding JSD:

``` r
getJSD (masterAlnSplit_freqs[["H2A"]][,7], masterAlnSplit_freqs[["H2A.B"]][,7] )
```

    ## [1] 0.8467096

Same thing for the 2nd position, where H2A and H2A.B are identical -
freqs in the H2A alignment:

``` r
masterAlnSplit_freqs[["H2A"]][,2]
```

    ## - A R N D C Q E G H I L K M F P S T W Y V 
    ## 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

and in the H2A.B alignment

``` r
masterAlnSplit_freqs[["H2A.B"]][,2]
```

    ## - A R N D C Q E G H I L K M F P S T W Y V 
    ## 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

here’s the corresponding JSD:

``` r
getJSD (masterAlnSplit_freqs[["H2A"]][,2], masterAlnSplit_freqs[["H2A.B"]][,2] )
```

    ## [1] 0

Now we make a wrapper function that gets JSD at every position between
two alignments

``` r
# I use the getJSD_multiClass in case I want to do it with >2 alignments
getJSDeachPos <- function(aln_list, letters=myAAtoTabulate) {
    # check none of the alignments are ragged
    null <- lapply(aln_list, checkAlnLengths)
    # check all the alns are the same length
    aln_widths <- sapply(aln_list, function(aln) { width(aln)[1] })
    if(length(unique(aln_widths))>1) {
        stop("\n\nThe alignments are not all the same length as each other\n\n")
    }
    # get freq matrices
    freqs_list <- lapply(aln_list, getAlnCounts, letters=letters, as.prob=TRUE)
    # for each column, get the JSDs
    allJSDs <- sapply(1:aln_widths[1], function(i) {
        freqs_this_position <- lapply(freqs_list, function(eachFreqMatrix) {
            eachFreqMatrix[,i]
        })
        jsd <- getJSD_multiClass(freqs_this_position)
        return(jsd)
    })
    return(allJSDs)
}

## test
# temp <- getJSDeachPos(list(masterAlnSplit[["H2A"]], 
#                            masterAlnSplit[["H2A.B"]],
#                            masterAlnSplit[["H2A.P"]]))
```

``` r
## we'll keep out JSD results in a tibble for later usage
jsd_tbl <- tibble(pos = 1:width(masterAln)[1],
                  aa_in_H2A = strsplit(as.character(masterAln[["H2A_mouse"]]), split="")[[1]]) %>% 
    unite("xaxis_label", c(pos, aa_in_H2A), sep="", remove=FALSE)

jsd_tbl$H2AvsB <- getJSDeachPos(list(masterAlnSplit[["H2A"]], masterAlnSplit[["H2A.B"]]))
jsd_tbl$H2AvsLBP <- getJSDeachPos(list(masterAlnSplit[["H2A"]], masterAlnSplit[["BandLandP"]]))
jsd_tbl$H2AvsLvsBvsP <- getJSDeachPos(list(masterAlnSplit[["H2A"]], 
                                           masterAlnSplit[["H2A.B"]], 
                                           masterAlnSplit[["H2A.L"]], 
                                           masterAlnSplit[["H2A.P"]]))
```

Make a plot

``` r
jsd_tbl %>% 
    ggplot(aes(x=pos, y=H2AvsLBP)) +
    geom_col(fill="lightgrey", col="darkgray", width=1) +
    theme_classic() +
    labs(x="alignment position, and amino acid in mouse H2A",
         y="JSD, H2A vs H2A.B+L+P") +
    theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=0.5)) + 
    scale_x_continuous(breaks=jsd_tbl$pos, labels=jsd_tbl$xaxis_label)
```

![](jensenShannonDistance_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

# Finished

show R version used, and package versions

``` r
sessionInfo()
```

    ## R version 4.4.0 (2024-04-24)
    ## Platform: x86_64-apple-darwin20
    ## Running under: macOS Ventura 13.7
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
    ##  [1] entropy_1.3.1       Biostrings_2.72.1   GenomeInfoDb_1.40.1
    ##  [4] XVector_0.44.0      IRanges_2.38.0      S4Vectors_0.42.0   
    ##  [7] BiocGenerics_0.50.0 lubridate_1.9.3     forcats_1.0.0      
    ## [10] stringr_1.5.1       dplyr_1.1.4         purrr_1.0.2        
    ## [13] readr_2.1.5         tidyr_1.3.1         tibble_3.2.1       
    ## [16] ggplot2_3.5.1       tidyverse_2.0.0    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] utf8_1.2.4              generics_0.1.3          stringi_1.8.4          
    ##  [4] hms_1.1.3               digest_0.6.36           magrittr_2.0.3         
    ##  [7] evaluate_0.24.0         grid_4.4.0              timechange_0.3.0       
    ## [10] fastmap_1.2.0           jsonlite_1.8.8          httr_1.4.7             
    ## [13] fansi_1.0.6             UCSC.utils_1.0.0        scales_1.3.0           
    ## [16] cli_3.6.3               crayon_1.5.3            rlang_1.1.4            
    ## [19] munsell_0.5.1           withr_3.0.0             yaml_2.3.8             
    ## [22] tools_4.4.0             tzdb_0.4.0              colorspace_2.1-0       
    ## [25] GenomeInfoDbData_1.2.12 vctrs_0.6.5             R6_2.5.1               
    ## [28] lifecycle_1.0.4         zlibbioc_1.50.0         pkgconfig_2.0.3        
    ## [31] pillar_1.9.0            gtable_0.3.5            glue_1.7.0             
    ## [34] highr_0.11              xfun_0.45               tidyselect_1.2.1       
    ## [37] rstudioapi_0.16.0       knitr_1.47              farver_2.1.2           
    ## [40] htmltools_0.5.8.1       labeling_0.4.3          rmarkdown_2.27         
    ## [43] compiler_4.4.0