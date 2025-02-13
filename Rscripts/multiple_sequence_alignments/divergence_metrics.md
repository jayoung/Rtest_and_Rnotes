divergence_metrics
================
Janet Young

2025-02-13

Goal - play with code to look at divergence metrics in a multiple
sequence alignment

Load libraries

``` r
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(Biostrings)
library(pwalign)
library(here)
```

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

Read example alignments (DNA and protein)

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

``` r
## get a smaller alignment
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

`stringDist()` counts number of changes between all pairs of sequences.
(I checked - hamming and levenshtein give the same result for this
particular alignment.)

``` r
aa_aln_h2a_l_firstBit %>% 
    stringDist(diag=TRUE, upper=TRUE, method="hamming")
```

    ##                      H2A.L_chineseHamster H2A.L_rat H2A.L_mouse H2A.L_pig
    ## H2A.L_chineseHamster                    0         6           7         7
    ## H2A.L_rat                               6         0           6         3
    ## H2A.L_mouse                             7         6           0         7
    ## H2A.L_pig                               7         3           7         0

To understand how gaps are counted, I make a very tiny alignment. I can
see that a gap-to-AA mismatch is counted as distance 1.

``` r
aa_aln_h2a_l_firstBit_tiny <- aa_aln_h2a_l_firstBit[2:3] %>% 
    narrow(start=8, end=11)
aa_aln_h2a_l_firstBit_tiny
```

    ## AAStringSet object of length 2:
    ##     width seq                                               names               
    ## [1]     4 LQFP                                              H2A.L_rat
    ## [2]     4 L--P                                              H2A.L_mouse

``` r
aa_aln_h2a_l_firstBit_tiny %>% 
    stringDist(diag=TRUE, upper=TRUE, method="hamming")
```

    ##             H2A.L_rat H2A.L_mouse
    ## H2A.L_rat           0           2
    ## H2A.L_mouse         2           0

Just double-checking - gap-to-gap is counted as distance 0:

``` r
aa_aln_h2a_l_firstBit_tiny_fakeAln <- aa_aln_h2a_l_firstBit_tiny[c(2,2)]
aa_aln_h2a_l_firstBit_tiny_fakeAln
```

    ## AAStringSet object of length 2:
    ##     width seq                                               names               
    ## [1]     4 L--P                                              H2A.L_mouse
    ## [2]     4 L--P                                              H2A.L_mouse

``` r
aa_aln_h2a_l_firstBit_tiny_fakeAln %>% 
    stringDist(diag=TRUE, upper=TRUE, method="hamming")
```

    ##             H2A.L_mouse H2A.L_mouse
    ## H2A.L_mouse           0           0
    ## H2A.L_mouse           0           0

``` r
# dna_aln_cenH3 <- here("Rscripts/multiple_sequence_alignments/example_alignment_files/cenH3_aln8.nt.fa") %>% 
#     readDNAStringSet()
# dna_aln_cenH3
```

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

xxxx

Play with more pwalign package functions - see
[documentation](https://bioconductor.org/packages/release/bioc/vignettes/pwalign/inst/doc/PairwiseAlignments.pdf)

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

    ## R version 4.4.0 (2024-04-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 18.04.6 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: FlexiBLAS OPENBLAS;  LAPACK version 3.11.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] here_1.0.1          pwalign_1.0.0       Biostrings_2.72.0  
    ##  [4] GenomeInfoDb_1.40.0 XVector_0.44.0      IRanges_2.38.0     
    ##  [7] S4Vectors_0.42.0    BiocGenerics_0.50.0 lubridate_1.9.4    
    ## [10] forcats_1.0.0       stringr_1.5.1       dplyr_1.1.4        
    ## [13] purrr_1.0.4         readr_2.1.5         tidyr_1.3.1        
    ## [16] tibble_3.2.1        ggplot2_3.5.1       tidyverse_2.0.0    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] generics_0.1.3          stringi_1.8.4           hms_1.1.3              
    ##  [4] digest_0.6.35           magrittr_2.0.3          evaluate_0.23          
    ##  [7] grid_4.4.0              timechange_0.3.0        fastmap_1.2.0          
    ## [10] rprojroot_2.0.4         jsonlite_1.8.9          httr_1.4.7             
    ## [13] UCSC.utils_1.0.0        scales_1.3.0            cli_3.6.3              
    ## [16] crayon_1.5.2            rlang_1.1.5             munsell_0.5.1          
    ## [19] withr_3.0.2             yaml_2.3.8              tools_4.4.0            
    ## [22] tzdb_0.4.0              colorspace_2.1-0        GenomeInfoDbData_1.2.12
    ## [25] vctrs_0.6.5             R6_2.5.1                lifecycle_1.0.4        
    ## [28] zlibbioc_1.50.0         pkgconfig_2.0.3         pillar_1.10.1          
    ## [31] gtable_0.3.5            glue_1.8.0              xfun_0.44              
    ## [34] tidyselect_1.2.1        rstudioapi_0.16.0       knitr_1.46             
    ## [37] htmltools_0.5.8.1       rmarkdown_2.26          compiler_4.4.0
