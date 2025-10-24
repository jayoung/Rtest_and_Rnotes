stringDist_questions_2025_10_22
================
Janet Young

2025-10-24

# Issues to look at, 2025/10/22

Using pwalign_1.4.0 and Biostrings_2.76.0 (these are both the current
release version)

[Issue 1](https://github.com/Bioconductor/pwalign/issues/15): can
StringDist handle gaps?

[Issue 2](https://github.com/Bioconductor/pwalign/issues/14): weird
scores on the diagonal

The two plain R docs stringDist_questions_2025_10_22.issue1.R and
stringDist_questions_2025_10_22.issue2.R are useful to make reprex.
After loading `library(reprex)`, you can highlight the code you want as
a reprex, copy it into the clipboard, and type `reprex()` on the
console.

## Setup

``` r
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(Biostrings)
library(pwalign)
library(here)
data(BLOSUM62)
```

Make small test alignments. `x` has no gaps, `y` has a gap.

``` r
aln_no_gaps <-  AAStringSet(c(seq1="VF",
                              seq2="VF",
                              seq3="LF"))

aln_with_gap <-  AAStringSet(c(seq1="VF",
                               seq2="V-",
                               seq3="LF"))
```

## Issue 1: can StringDist handle gaps?

stringDist works fine on the gapped alignment with default method
(levenshtein)

``` r
stringDist(aln_with_gap,
           diag=TRUE, upper=TRUE)
```

    ##      seq1 seq2 seq3
    ## seq1    0    1    1
    ## seq2    1    0    2
    ## seq3    1    2    0

substitutionMatrix method does work without gaps:

``` r
stringDist(aln_no_gaps, 
           diag=TRUE, upper=TRUE,
           method="substitutionMatrix",
           substitutionMatrix=BLOSUM62)
```

    ##      seq1 seq2 seq3
    ## seq1    0  -10   -7
    ## seq2  -10    0   -7
    ## seq3   -7   -7    0

But substitutionMatrix method gives an error when there’s a gap. Maybe
the gap penalties are not being handled:

``` r
stringDist(aln_with_gap, 
           diag=TRUE, upper=TRUE,
           method="substitutionMatrix",
           substitutionMatrix=BLOSUM62)
```

    ## Error in .Call2("XStringSet_align_distance", x, type, typeCode, gapOpening, : key 45 not in lookup table

## Issue 2: inconsistent scores on the diagonal

seq1 and seq2 are identical to each other. Therefore the scores for
seq1-versus-seq1 and seq1-versus-seq2 should be identical to each other.
But they’re not. Self-matches always get a score of 0. This makes sense
for many distance metrics but maybe not for substitutionMatrix methods.

``` r
stringDist(aln_no_gaps, 
           diag=TRUE, upper=TRUE,
           method="substitutionMatrix",
           substitutionMatrix=BLOSUM62)
```

    ##      seq1 seq2 seq3
    ## seq1    0  -10   -7
    ## seq2  -10    0   -7
    ## seq3   -7   -7    0

Just to make sure I understand how the off-diagonal scores come about,
here are the relevant bits of the BLOSUM62 matrix:

``` r
BLOSUM62[c("F","L","V"),c("F","L","V")]
```

    ##    F L  V
    ## F  6 0 -1
    ## L  0 4  1
    ## V -1 1  4

Understanding the scores:

- VF aligned to VF scores -10, because BLOSUM62 scores F-F as 6 and V-V
  as 4
- VF aligned to LF score -7, because BLOSUM62 scores F-F as 6 and V-L as
  1

# Herve’s workaround

<https://github.com/Bioconductor/pwalign/issues/15>

I put his new functions in
`useful_functions/multiple_sequence_alignments_functions.R`

``` r
source(here("useful_functions/multiple_sequence_alignments_functions.R"))
```

``` r
x <- c("ARND", "AR-D", "CRNA")
x1 <- c(x, x[1:2])
x1 %>% 
    AAStringSet()
```

    ## AAStringSet object of length 5:
    ##     width seq
    ## [1]     4 ARND
    ## [2]     4 AR-D
    ## [3]     4 CRNA
    ## [4]     4 ARND
    ## [5]     4 AR-D

``` r
alignedStringDistMatrix(x1)
```

    ##      [,1] [,2] [,3] [,4] [,5]
    ## [1,]    0    1    2    0    1
    ## [2,]    1    0    3    1    0
    ## [3,]    2    3    0    2    3
    ## [4,]    0    1    2    0    1
    ## [5,]    1    0    3    1    0

``` r
alignedStringDistMatrix(x1, indel.weight=10)
```

    ##      [,1] [,2] [,3] [,4] [,5]
    ## [1,]    0   10    2    0   10
    ## [2,]   10    0   12   10    0
    ## [3,]    2   12    0    2   12
    ## [4,]    0   10    2    0   10
    ## [5,]   10    0   12   10    0

``` r
alignedStringDistMatrix(x1, indel.weight=0)
```

    ##      [,1] [,2] [,3] [,4] [,5]
    ## [1,]    0    0    2    0    0
    ## [2,]    0    0    2    0    0
    ## [3,]    2    2    0    2    2
    ## [4,]    0    0    2    0    0
    ## [5,]    0    0    2    0    0

``` r
alignedStringDistMatrix(x1, weightmat=-BLOSUM62)
```

    ##      [,1] [,2] [,3] [,4] [,5]
    ## [1,]  -21  -14   -9  -21  -14
    ## [2,]  -14  -15   -2  -14  -15
    ## [3,]   -9   -2  -24   -9   -2
    ## [4,]  -21  -14   -9  -21  -14
    ## [5,]  -14  -15   -2  -14  -15

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
    ##  [1] here_1.0.2          pwalign_1.4.0       Biostrings_2.76.0  
    ##  [4] GenomeInfoDb_1.44.1 XVector_0.48.0      IRanges_2.42.0     
    ##  [7] S4Vectors_0.46.0    BiocGenerics_0.54.0 generics_0.1.4     
    ## [10] lubridate_1.9.4     forcats_1.0.0       stringr_1.5.2      
    ## [13] dplyr_1.1.4         purrr_1.1.0         readr_2.1.5        
    ## [16] tidyr_1.3.1         tibble_3.3.0        ggplot2_3.5.2      
    ## [19] tidyverse_2.0.0    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] stringi_1.8.7           hms_1.1.3               digest_0.6.37          
    ##  [4] magrittr_2.0.4          evaluate_1.0.5          grid_4.5.1             
    ##  [7] timechange_0.3.0        RColorBrewer_1.1-3      fastmap_1.2.0          
    ## [10] rprojroot_2.1.1         jsonlite_2.0.0          httr_1.4.7             
    ## [13] UCSC.utils_1.4.0        scales_1.4.0            cli_3.6.5              
    ## [16] rlang_1.1.6             crayon_1.5.3            withr_3.0.2            
    ## [19] yaml_2.3.10             tools_4.5.1             tzdb_0.5.0             
    ## [22] GenomeInfoDbData_1.2.14 vctrs_0.6.5             R6_2.6.1               
    ## [25] lifecycle_1.0.4         pkgconfig_2.0.3         pillar_1.11.1          
    ## [28] gtable_0.3.6            glue_1.8.0              xfun_0.53              
    ## [31] tidyselect_1.2.1        rstudioapi_0.17.1       knitr_1.50             
    ## [34] farver_2.1.2            htmltools_0.5.8.1       rmarkdown_2.29         
    ## [37] compiler_4.5.1
