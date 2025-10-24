multiple_sequence_alignments_functions_testCode.Rmd
================
Janet Young

2025-10-22

# Goal

Make a test alignment and some code to work on how to treat alignments
with gaps and ambiguities

Set up test DNA alignment, with in-frame gap, frameshifting gap and
codons containing Ns

``` r
testAln <- DNAStringSet(c(seq1="ATG-ACAAT---TGG",
                          seq2="ATNGANGGT---ATG",
                          seq3="ATG-ACGGTTGGATG"))
testAln
```

    ## DNAStringSet object of length 3:
    ##     width seq                                               names               
    ## [1]    15 ATG-ACAAT---TGG                                   seq1
    ## [2]    15 ATNGANGGT---ATG                                   seq2
    ## [3]    15 ATG-ACGGTTGGATG                                   seq3

`getCodons` is a small function that splits each sequence into
constituent codons:

``` r
tempCodons <- getCodons(testAln)
tempCodons
```

    ## $seq1
    ## [1] "ATG" "-AC" "AAT" "---" "TGG"
    ## 
    ## $seq2
    ## [1] "ATN" "GAN" "GGT" "---" "ATG"
    ## 
    ## $seq3
    ## [1] "ATG" "-AC" "GGT" "TGG" "ATG"

Translate the alignment using my `translateGappedAln` function. The
defaults are:

- `frameshiftTranslatesTo="X"` (e.g. for codons containing gap and DNA
  base),
- `unknownCodonTranslatesTo="-"` (e.g. for codons containing
  ambiguities)
- `quiet=FALSE` : will emit warnings

``` r
translateGappedAln(testAln)
```

    ## Warning in FUN(X[[i]], ...): 
    ## 
    ## WARNING - seq seq2 contains ambiguous codons:1_ATN,2_GAN

    ## 
    ## Warning in seq seq1 - there were codons I could not translate. Using this character: -

    ## The codons in question were: -AC

    ## 
    ## Warning in seq seq2 - there were codons I could not translate. Using this character: -

    ## The codons in question were: ATN,GAN

    ## 
    ## Warning in seq seq3 - there were codons I could not translate. Using this character: -

    ## The codons in question were: -AC

    ## AAStringSet object of length 3:
    ##     width seq                                               names               
    ## [1]     5 MXN-W                                             seq1
    ## [2]     5 --G-M                                             seq2
    ## [3]     5 MXGWM                                             seq3

Can simply turn off the warnings:

``` r
translateGappedAln(testAln, quiet=TRUE)
```

    ## AAStringSet object of length 3:
    ##     width seq                                               names               
    ## [1]     5 MXN-W                                             seq1
    ## [2]     5 --G-M                                             seq2
    ## [3]     5 MXGWM                                             seq3

Can choose other translations for the ambiguous codons and frameshifts:

``` r
translateGappedAln(testAln,
                   frameshiftTranslatesTo="b", 
                   unknownCodonTranslatesTo="z", 
                   quiet=TRUE)
```

    ## AAStringSet object of length 3:
    ##     width seq                                               names               
    ## [1]     5 MBN-W                                             seq1
    ## [2]     5 ZZG-M                                             seq2
    ## [3]     5 MBGWM                                             seq3

``` r
translateGappedAln(testAln,
                   frameshiftTranslatesTo="-", 
                   unknownCodonTranslatesTo="-", 
                   quiet=TRUE)
```

    ## AAStringSet object of length 3:
    ##     width seq                                               names               
    ## [1]     5 M-N-W                                             seq1
    ## [2]     5 --G-M                                             seq2
    ## [3]     5 M-GWM                                             seq3

# Degapping alignments

Reminder of what the alignment looks like:

``` r
testAln
```

    ## DNAStringSet object of length 3:
    ##     width seq                                               names               
    ## [1]    15 ATG-ACAAT---TGG                                   seq1
    ## [2]    15 ATNGANGGT---ATG                                   seq2
    ## [3]    15 ATG-ACGGTTGGATG                                   seq3

``` r
translateGappedAln(testAln, quiet=TRUE)
```

    ## AAStringSet object of length 3:
    ##     width seq                                               names               
    ## [1]     5 MXN-W                                             seq1
    ## [2]     5 --G-M                                             seq2
    ## [3]     5 MXGWM                                             seq3

``` r
degapAln(translateGappedAln(testAln, quiet=TRUE),
         fractionOfSeqsWithGap = 0.5)
```

    ## AAStringSet object of length 3:
    ##     width seq                                               names               
    ## [1]     4 MXNW                                              seq1
    ## [2]     4 --GM                                              seq2
    ## [3]     4 MXGM                                              seq3

``` r
degapAln(testAln,
         fractionOfSeqsWithGap = 0.5)
```

    ## DNAStringSet object of length 3:
    ##     width seq                                               names               
    ## [1]    11 ATGACAATTGG                                       seq1
    ## [2]    11 ATNANGGTATG                                       seq2
    ## [3]    11 ATGACGGTATG                                       seq3

# Finished

``` r
sessionInfo()
```

    ## R version 4.5.1 (2025-06-13)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Sequoia 15.6.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
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
    ## [1] Biostrings_2.76.0   GenomeInfoDb_1.44.1 XVector_0.48.0     
    ## [4] IRanges_2.42.0      S4Vectors_0.46.0    BiocGenerics_0.54.0
    ## [7] generics_0.1.4      magrittr_2.0.4      here_1.0.2         
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] crayon_1.5.3            httr_1.4.7              cli_3.6.5              
    ##  [4] knitr_1.50              rlang_1.1.6             xfun_0.53              
    ##  [7] UCSC.utils_1.4.0        jsonlite_2.0.0          rprojroot_2.1.1        
    ## [10] htmltools_0.5.8.1       rmarkdown_2.29          evaluate_1.0.5         
    ## [13] fastmap_1.2.0           yaml_2.3.10             compiler_4.5.1         
    ## [16] rstudioapi_0.17.1       digest_0.6.37           R6_2.6.1               
    ## [19] GenomeInfoDbData_1.2.14 tools_4.5.1
