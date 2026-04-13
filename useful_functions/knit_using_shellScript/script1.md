script1
================
Janet Young

2026-04-13

``` r
mtcars_efficient <- mtcars[which(mtcars$mpg >= 30),]
```

This script filters the `mtcars` dataset (contains 32 rows) to keep only
cars where mpg is \>= 30. There are 4 rows after this filtering.

Save `mtcars_efficient` for use in script2.

``` r
# it is saved in the same location as the script
save(mtcars_efficient, file="mtcars_efficient.Rdata")
```

Show .libPaths

``` r
.libPaths()
```

    ## [1] "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library"

Show R and package version information

``` r
sessionInfo()
```

    ## R version 4.5.3 (2026-03-11)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Tahoe 26.4.1
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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_4.5.3    fastmap_1.2.0     cli_3.6.5         tools_4.5.3      
    ##  [5] htmltools_0.5.9   otel_0.2.0        rstudioapi_0.18.0 yaml_2.3.12      
    ##  [9] rmarkdown_2.31    knitr_1.51        xfun_0.57         digest_0.6.39    
    ## [13] rlang_1.1.7       evaluate_1.0.5
