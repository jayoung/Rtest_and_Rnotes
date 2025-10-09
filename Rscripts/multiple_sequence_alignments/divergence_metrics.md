divergence_metrics
================
Janet Young

2025-10-08

# Goal

Play with code to look at divergence metrics in a multiple sequence
alignment

``` r
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(Biostrings)
library(pwalign)
library(here)
```

# Define some utility functions

Set up a couple of utility functions to degap alignments

``` r
#### degapNucAln - function to remove columns that are entirely gap (or can use fractionOfSeqsWithGap to relax requirements of how many seqs have gap)
degapNucAln <- function(myAln, fractionOfSeqsWithGap=1) {
    maskedAln <- myAln %>% 
        DNAMultipleAlignment() %>% 
        maskGaps(min.fraction=fractionOfSeqsWithGap, 
                 min.block.width=1) %>% 
        DNAStringSet()
    return(maskedAln)
}


#### degapAAaln - function to remove columns that are entirely gap (or can use fractionOfSeqsWithGap to relax requirements of how many seqs have gap)
### this is how I can see maskGaps uses >= for the fractionOfSeqsWithGap threshold:
# selectMethod("maskGaps", "AAMultipleAlignment")
# newmask <- (m["-", ]/colSums(m)) >= min.fraction
degapAAaln <- function(myAln, fractionOfSeqsWithGap=1) {
    maskedAln <- myAln %>% 
        AAMultipleAlignment() %>% 
        maskGaps(min.fraction=fractionOfSeqsWithGap, 
                 min.block.width=1) %>% 
        AAStringSet()
    return(maskedAln)
}
```

# Read example alignments (DNA and protein)

Histon amino acid alignment

``` r
aa_aln_histone <- here("Rscripts/multiple_sequence_alignments/example_alignment_files/exampleProtAln_shortH2As_histoneFoldDomain.fa") %>% 
    readAAStringSet()

names(aa_aln_histone) <- sapply(strsplit(names(aa_aln_histone), " "), "[[", 1)

aa_aln_histone
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

Get a smaller alignment, for testing

``` r
aa_aln_h2a_l <- aa_aln_histone[ grep("H2A.L", names(aa_aln_histone)) ] %>% 
    degapAAaln()

aa_aln_h2a_l_firstBit <- aa_aln_h2a_l %>% 
    narrow(start=1, end=12)
aa_aln_h2a_l_firstBit <- aa_aln_h2a_l_firstBit[1:4]

aa_aln_h2a_l_firstBit
```

    ## AAStringSet object of length 4:
    ##     width seq                                               names               
    ## [1]    12 SRITRGQLQFSL                                      H2A.L_chineseHamster
    ## [2]    12 SRSRRAELQFPV                                      H2A.L_rat
    ## [3]    12 TRSQRGEL--PL                                      H2A.L_mouse
    ## [4]    12 SCSSRAELQFPM                                      H2A.L_pig

## stringDist function

`stringDist()` (was in Biostrings, moved to pwalign package)

It counts number of changes between all pairs of sequences, using
various methods - “levenshtein”, “hamming”, “quality”, or
“substitutionMatrix”.

(I checked - hamming and levenshtein give the same result for this
particular alignment.)

Show Hamming distances for the tiny alignment:

``` r
aa_aln_h2a_l_firstBit %>% 
    stringDist(diag=TRUE, upper=TRUE, method="hamming")
```

    ##                      H2A.L_chineseHamster H2A.L_rat H2A.L_mouse H2A.L_pig
    ## H2A.L_chineseHamster                    0         6           7         7
    ## H2A.L_rat                               6         0           6         3
    ## H2A.L_mouse                             7         6           0         7
    ## H2A.L_pig                               7         3           7         0

Further below I explore the other distance metrics.

## Understand how gaps are counted

To understand how gaps are counted, I make a very tiny alignment. I can
see that a gap-to-AA mismatch is counted as distance 1. I can also see
that gap-to-gap positions count as 0 distance.

``` r
aa_aln_h2a_l_firstBit_tiny <- aa_aln_h2a_l_firstBit[2:3] %>% 
    narrow(start=8, end=11)

temp <- aa_aln_h2a_l_firstBit_tiny[2]
names(temp) <- "H2A.L_mouse_again"
aa_aln_h2a_l_firstBit_tiny <- c(aa_aln_h2a_l_firstBit_tiny,
                                temp)

aa_aln_h2a_l_firstBit_tiny
```

    ## AAStringSet object of length 3:
    ##     width seq                                               names               
    ## [1]     4 LQFP                                              H2A.L_rat
    ## [2]     4 L--P                                              H2A.L_mouse
    ## [3]     4 L--P                                              H2A.L_mouse_again

### gap handling by stringDist

``` r
aa_aln_h2a_l_firstBit_tiny %>% 
    stringDist(diag=TRUE, upper=TRUE, method="hamming")
```

    ##                   H2A.L_rat H2A.L_mouse H2A.L_mouse_again
    ## H2A.L_rat                 0           2                 2
    ## H2A.L_mouse               2           0                 0
    ## H2A.L_mouse_again         2           0                 0

# Hamming versus Levenshtein distances

This [wikipedia
page](https://en.wikipedia.org/wiki/Levenshtein_distance) explains the
difference between hamming and Levenshtein distances. Levenshtein
distance essentially allows indels.

Example - lawn versus flaw:

``` r
y <- c("lawn","flaw")
x <- BStringSet(y)
names(x) <- y
x
```

    ## BStringSet object of length 2:
    ##     width seq                                               names               
    ## [1]     4 lawn                                              lawn
    ## [2]     4 flaw                                              flaw

``` r
x %>% stringDist(diag=TRUE, upper=TRUE, method="hamming")
```

    ##      lawn flaw
    ## lawn    0    4
    ## flaw    4    0

``` r
x %>% stringDist(diag=TRUE, upper=TRUE, method="levenshtein")
```

    ##      lawn flaw
    ## lawn    0    2
    ## flaw    2    0

Show levenshtein distances for the real tiny alignment:

``` r
aa_aln_h2a_l_firstBit %>% 
    stringDist(diag=TRUE, upper=TRUE, method="levenshtein")
```

    ##                      H2A.L_chineseHamster H2A.L_rat H2A.L_mouse H2A.L_pig
    ## H2A.L_chineseHamster                    0         6           7         7
    ## H2A.L_rat                               6         0           6         3
    ## H2A.L_mouse                             7         6           0         7
    ## H2A.L_pig                               7         3           7         0

Haven’t run substitutionMatrix distances, because it requires providing
a matrix using substitutionMatrix argument, and I haven’t tracked down
an amino acid matrix that works. This code failes:

``` r
# data(BLOSUM62)
# aa_aln_h2a_l_firstBit %>% 
#     stringDist(diag=TRUE, upper=TRUE, 
#                method="substitutionMatrix",
#                substitutionMatrix=BLOSUM62)
# Error in .Call2("XStringSet_align_distance", x, type, typeCode, gapOpening,  : 
#   key 45 not in lookup table
```

Show quality distances for the tiny alignment. I have no idea what this
is - I can only imagine it’s something intended for phred scores and
it’s inappropriate to use here.

``` r
aa_aln_h2a_l_firstBit %>% 
    stringDist(diag=TRUE, upper=TRUE, method="quality")
```

    ##                      H2A.L_chineseHamster H2A.L_rat H2A.L_mouse H2A.L_pig
    ## H2A.L_chineseHamster              0.00000  13.82201     7.51834   7.51834
    ## H2A.L_rat                        13.82201   0.00000    13.82201  32.73302
    ## H2A.L_mouse                       7.51834  13.82201     0.00000   7.51834
    ## H2A.L_pig                         7.51834  32.73302     7.51834   0.00000

# Play with `pwalign` package

See
[documentation](https://bioconductor.org/packages/release/bioc/vignettes/pwalign/inst/doc/PairwiseAlignments.pdf)

`PairwiseAlignmentsSingleSubject()` :

``` r
temp <- aa_aln_h2a_l_firstBit[1:2]
temp2 <- temp %>% PairwiseAlignmentsSingleSubject()
nmatch(temp2)
```

    ## [1] 6

``` r
nmismatch(temp2)
```

    ## [1] 6

``` r
# mismatchTable(temp2)
deletion(temp2)
```

    ## IRangesList object of length 1:
    ## [[1]]
    ## IRanges object with 0 ranges and 0 metadata columns:
    ##        start       end     width
    ##    <integer> <integer> <integer>

``` r
insertion(temp2)
```

    ## IRangesList object of length 1:
    ## [[1]]
    ## IRanges object with 0 ranges and 0 metadata columns:
    ##        start       end     width
    ##    <integer> <integer> <integer>

``` r
indel(temp2)
```

    ## An object of class "InDel"
    ## Slot "insertion":
    ## IRangesList object of length 1:
    ## [[1]]
    ## IRanges object with 0 ranges and 0 metadata columns:
    ##        start       end     width
    ##    <integer> <integer> <integer>
    ## 
    ## 
    ## Slot "deletion":
    ## IRangesList object of length 1:
    ## [[1]]
    ## IRanges object with 0 ranges and 0 metadata columns:
    ##        start       end     width
    ##    <integer> <integer> <integer>

``` r
nchar(temp2)
```

    ## [1] 12

``` r
pid(temp2)
```

    ## [1] 50

# Finished

``` r
sessionInfo()
```

    ## R version 4.5.1 (2025-06-13)
    ## Platform: x86_64-apple-darwin20
    ## Running under: macOS Sequoia 15.7.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.5-x86_64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
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
    ##  [1] here_1.0.1          pwalign_1.4.0       Biostrings_2.76.0  
    ##  [4] GenomeInfoDb_1.44.0 XVector_0.48.0      IRanges_2.42.0     
    ##  [7] S4Vectors_0.46.0    BiocGenerics_0.54.0 generics_0.1.4     
    ## [10] lubridate_1.9.4     forcats_1.0.0       stringr_1.5.1      
    ## [13] dplyr_1.1.4         purrr_1.0.4         readr_2.1.5        
    ## [16] tidyr_1.3.1         tibble_3.3.0        ggplot2_3.5.2      
    ## [19] tidyverse_2.0.0    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] stringi_1.8.7           hms_1.1.3               digest_0.6.37          
    ##  [4] magrittr_2.0.3          evaluate_1.0.4          grid_4.5.1             
    ##  [7] timechange_0.3.0        RColorBrewer_1.1-3      fastmap_1.2.0          
    ## [10] rprojroot_2.0.4         jsonlite_2.0.0          httr_1.4.7             
    ## [13] UCSC.utils_1.4.0        scales_1.4.0            cli_3.6.5              
    ## [16] crayon_1.5.3            rlang_1.1.6             withr_3.0.2            
    ## [19] yaml_2.3.10             tools_4.5.1             tzdb_0.5.0             
    ## [22] GenomeInfoDbData_1.2.14 vctrs_0.6.5             R6_2.6.1               
    ## [25] lifecycle_1.0.4         pkgconfig_2.0.3         pillar_1.10.2          
    ## [28] gtable_0.3.6            glue_1.8.0              xfun_0.52              
    ## [31] tidyselect_1.2.1        rstudioapi_0.17.1       knitr_1.50             
    ## [34] dichromat_2.0-0.1       farver_2.1.2            htmltools_0.5.8.1      
    ## [37] rmarkdown_2.29          compiler_4.5.1
