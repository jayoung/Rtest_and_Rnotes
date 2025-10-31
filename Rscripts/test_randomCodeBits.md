Untitled
================
2025-06-12

# Goal

Temporary script to test things I’m trying to figure out

``` r
# gr <- GRanges(seqnames="NC_001133.9",
#               ranges=IRanges(start=c(10,20,30,40)),
#               score=c(5,2,-4, 10))
# gr
```

``` r
# y <- coverage(gr, weight = "score")
# y[which(y==0)] <- NA
# y
```

``` r
# export(y, con=here("temp/test.bw"))
```

``` r
# import(con=here("temp/test.bw"))
```

``` r
# import(con="/Volumes/malik_h/user/jayoung/forOtherPeople/forGrantKing/SATAY_stuff/SATAY/240103_SATAY_Pilot_4.0/data/bedfiles/fitness_raw/A_cerulenin_fitness0cent_raw.bw")
```

``` r
# import(con="/Volumes/malik_h/user/jayoung/forOtherPeople/forGrantKing/SATAY_stuff/SATAY/240103_SATAY_Pilot_4.0/data/bedfiles/fitness_raw/A_cerulenin_fitness0cent_raw.bw")
```

``` r
gr <- GRanges(seqnames=rep("chr1",5),
        ranges=IRanges(start=c(5,10,15,20,25),width=1), 
        score=c(0.72869552845629478, 1.12030172459286, 1.4891584, 1.0056807, 0.9010696) )
gr
```

    ## GRanges object with 5 ranges and 1 metadata column:
    ##       seqnames    ranges strand |     score
    ##          <Rle> <IRanges>  <Rle> | <numeric>
    ##   [1]     chr1         5      * |  0.728696
    ##   [2]     chr1        10      * |  1.120302
    ##   [3]     chr1        15      * |  1.489158
    ##   [4]     chr1        20      * |  1.005681
    ##   [5]     chr1        25      * |  0.901070
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome; no seqlengths

``` r
cov <- coverage(gr, weight="score") 

as.numeric(cov[["chr1"]])
```

    ##  [1] 0.0000000 0.0000000 0.0000000 0.0000000 0.7286955 0.0000000 0.0000000
    ##  [8] 0.0000000 0.0000000 1.1203017 0.0000000 0.0000000 0.0000000 0.0000000
    ## [15] 1.4891584 0.0000000 0.0000000 0.0000000 0.0000000 1.0056807 0.0000000
    ## [22] 0.0000000 0.0000000 0.0000000 0.9010696

``` r
sessionInfo()
```

    ## R version 4.4.1 (2024-06-14)
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
    ## [1] here_1.0.1           rtracklayer_1.64.0   GenomicRanges_1.56.2
    ## [4] GenomeInfoDb_1.40.0  IRanges_2.38.1       S4Vectors_0.42.0    
    ## [7] BiocGenerics_0.50.0 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Matrix_1.7-0                jsonlite_1.8.9             
    ##  [3] compiler_4.4.1              rjson_0.2.21               
    ##  [5] crayon_1.5.3                SummarizedExperiment_1.34.0
    ##  [7] Biobase_2.64.0              Rsamtools_2.20.0           
    ##  [9] bitops_1.0-7                Biostrings_2.72.0          
    ## [11] GenomicAlignments_1.40.0    parallel_4.4.1             
    ## [13] BiocParallel_1.38.0         yaml_2.3.8                 
    ## [15] fastmap_1.2.0               lattice_0.22-6             
    ## [17] R6_2.6.1                    XVector_0.44.0             
    ## [19] S4Arrays_1.4.0              curl_7.0.0                 
    ## [21] knitr_1.50                  XML_3.99-0.16.1            
    ## [23] DelayedArray_0.30.1         MatrixGenerics_1.16.0      
    ## [25] rprojroot_2.0.4             GenomeInfoDbData_1.2.12    
    ## [27] rlang_1.1.5                 xfun_0.53                  
    ## [29] SparseArray_1.4.3           cli_3.6.3                  
    ## [31] zlibbioc_1.50.0             grid_4.4.1                 
    ## [33] digest_0.6.36               evaluate_0.24.0            
    ## [35] codetools_0.2-20            abind_1.4-5                
    ## [37] RCurl_1.98-1.14             restfulr_0.0.15            
    ## [39] rmarkdown_2.27              httr_1.4.7                 
    ## [41] matrixStats_1.3.0           tools_4.4.1                
    ## [43] BiocIO_1.14.0               htmltools_0.5.8.1          
    ## [45] UCSC.utils_1.0.0
