bioconductor_demo
================
Janet Young

2026-03-16

``` r
## the above is a good chunk header for chunks that load libraries
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(here)

library(GenomicRanges)
library(GenomicAlignments)

# example data
library(pasillaBamSubset)
```

# Coverage demo

Explore code from [GenomicRanges
vignette](https://bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesHOWTOs.pdf)

We use example data from the pasillaBamSubset - a bam file containing
reads mapped to the Drosophila genome.

First, get the bam file name (`un1`), read it in as a `GAlignments`
object (`reads1`), and get coverage at each base (`cvg1`) using
`GenomicAlignments::coverage()` function.

The result is an RleList object (one item per chromosome). Can get a
similar coverage object by importing a bigwig file.

``` r
un1 <- untreated1_chr4() 
reads1 <- readGAlignments(un1)
cvg1 <- coverage(reads1)
```

This RleList stores coverage for each chromosome, each as a simple Rle
object (see `?Rle`):

``` r
cvg1[["chr4"]]
```

    ## integer-Rle of length 1351857 with 122061 runs
    ##   Lengths:  891   27    5   12   13   45    5 ...    3  106   75 1600   75 1659
    ##   Values :    0    1    2    3    4    5    4 ...    6    0    1    0    1    0

We can get genome info (chromosome lengths) and turn it into a GRanges
object representing each chromosome’s full length: `genome_gr`

``` r
genome_gr <- reads1 |> seqinfo() |> GRanges()
genome_gr
```

    ## GRanges object with 8 ranges and 0 metadata columns:
    ##           seqnames     ranges strand
    ##              <Rle>  <IRanges>  <Rle>
    ##     chr2L    chr2L 1-23011544      *
    ##     chr2R    chr2R 1-21146708      *
    ##     chr3L    chr3L 1-24543557      *
    ##     chr3R    chr3R 1-27905053      *
    ##      chr4     chr4  1-1351857      *
    ##      chrM     chrM    1-19517      *
    ##      chrX     chrX 1-22422827      *
    ##   chrYHet  chrYHet   1-347038      *
    ##   -------
    ##   seqinfo: 8 sequences from an unspecified genome

Use `genome_gr` to get sliding windows for the whole genome (1Mb,
sliding by 0.5Mb):

``` r
genome_1mb_windows <- slidingWindows(genome_gr, width=10^6, step=0.5*10^6) |> 
    unlist()
# The names for each window can be annoying later so I'll strip them off
names(genome_1mb_windows) <- NULL
```

Now we can use the `Views` function to consider coverage for each 1mb
window separately

``` r
Views(cvg1, genome_1mb_windows)
```

    ## RleViewsList object of length 8:
    ## $chr2L
    ## Views on a 23011544-length Rle subject
    ## 
    ## views:
    ##         start      end   width
    ##  [1]        1  1000000   1e+06 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...]
    ##  [2]   500001  1500000   1e+06 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...]
    ##  [3]  1000001  2000000   1e+06 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...]
    ##  [4]  1500001  2500000   1e+06 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...]
    ##  [5]  2000001  3000000   1e+06 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...]
    ##  [6]  2500001  3500000   1e+06 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...]
    ##  [7]  3000001  4000000   1e+06 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...]
    ##  [8]  3500001  4500000   1e+06 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...]
    ##  [9]  4000001  5000000   1e+06 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...]
    ##  ...      ...      ...     ... ...
    ## [38] 18500001 19500000   1e+06 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...]
    ## [39] 19000001 20000000   1e+06 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...]
    ## [40] 19500001 20500000   1e+06 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...]
    ## [41] 20000001 21000000   1e+06 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...]
    ## [42] 20500001 21500000   1e+06 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...]
    ## [43] 21000001 22000000   1e+06 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...]
    ## [44] 21500001 22500000   1e+06 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...]
    ## [45] 22000001 23000000   1e+06 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...]
    ## [46] 22500001 23011544  511544 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...]
    ## 
    ## ...
    ## <7 more elements>

And we can feed the result into the `viewMeans()` function to get mean
for each window (it returns a list, one per chromosome).

There’s also `viewSums()`, `viewMins()`, `viewMaxs()`, and then
`viewApply` to do more complex functions for each window.

``` r
Views(cvg1, genome_1mb_windows) |> 
    viewMeans()
```

    ## NumericList of length 8
    ## [["chr2L"]] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ... 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    ## [["chr2R"]] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ... 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    ## [["chr3L"]] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ... 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    ## [["chr3R"]] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ... 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    ## [["chr4"]] 12.9399 8.31981776284048
    ## [["chrM"]] 0
    ## [["chrX"]] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ... 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    ## [["chrYHet"]] 0

I like to store the output (mean coverage) as a column on the original
GRanges object - the `unlist` function is useful here:

``` r
genome_1mb_windows$cov_mean <- Views(cvg1, genome_1mb_windows) |> 
    viewMeans() |> 
    unlist(use.names = FALSE)
```

Show results for just the chr4 windows:

``` r
genome_1mb_windows[which(seqnames(genome_1mb_windows)=="chr4")]
```

    ## GRanges object with 2 ranges and 1 metadata column:
    ##       seqnames         ranges strand |  cov_mean
    ##          <Rle>      <IRanges>  <Rle> | <numeric>
    ##   [1]     chr4      1-1000000      * |  12.93990
    ##   [2]     chr4 500001-1351857      * |   8.31982
    ##   -------
    ##   seqinfo: 8 sequences from an unspecified genome

# Finished

``` r
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Tahoe 26.3.1
    ## 
    ## Matrix products: default
    ## BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
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
    ##  [1] pasillaBamSubset_0.48.0     GenomicAlignments_1.46.0   
    ##  [3] Rsamtools_2.26.0            Biostrings_2.78.0          
    ##  [5] XVector_0.50.0              SummarizedExperiment_1.40.0
    ##  [7] Biobase_2.70.0              MatrixGenerics_1.22.0      
    ##  [9] matrixStats_1.5.0           GenomicRanges_1.62.1       
    ## [11] Seqinfo_1.0.0               IRanges_2.44.0             
    ## [13] S4Vectors_0.48.0            BiocGenerics_0.56.0        
    ## [15] generics_0.1.4              here_1.0.2                 
    ## [17] lubridate_1.9.5             forcats_1.0.1              
    ## [19] stringr_1.6.0               dplyr_1.2.0                
    ## [21] purrr_1.2.1                 readr_2.1.6                
    ## [23] tidyr_1.3.2                 tibble_3.3.1               
    ## [25] ggplot2_4.0.2               tidyverse_2.0.0            
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6        xfun_0.56           lattice_0.22-7     
    ##  [4] tzdb_0.5.0          vctrs_0.7.1         tools_4.5.2        
    ##  [7] bitops_1.0-9        parallel_4.5.2      pkgconfig_2.0.3    
    ## [10] Matrix_1.7-4        RColorBrewer_1.1-3  S7_0.2.1           
    ## [13] cigarillo_1.0.0     lifecycle_1.0.5     compiler_4.5.2     
    ## [16] farver_2.1.2        codetools_0.2-20    htmltools_0.5.9    
    ## [19] yaml_2.3.12         pillar_1.11.1       crayon_1.5.3       
    ## [22] BiocParallel_1.44.0 DelayedArray_0.36.0 abind_1.4-8        
    ## [25] tidyselect_1.2.1    digest_0.6.39       stringi_1.8.7      
    ## [28] rprojroot_2.1.1     fastmap_1.2.0       grid_4.5.2         
    ## [31] cli_3.6.5           SparseArray_1.10.6  magrittr_2.0.4     
    ## [34] S4Arrays_1.10.1     withr_3.0.2         scales_1.4.0       
    ## [37] timechange_0.4.0    rmarkdown_2.30      otel_0.2.0         
    ## [40] hms_1.1.4           evaluate_1.0.5      knitr_1.51         
    ## [43] rlang_1.1.7         glue_1.8.0          rstudioapi_0.18.0  
    ## [46] R6_2.6.1
