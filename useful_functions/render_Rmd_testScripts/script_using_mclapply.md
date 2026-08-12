script_using_mclapply
================
Janet Young

2026-08-11

``` r
num_interations <- 10^5
num_threads <- 4
```

# Do something slow using lapply

`?system.time` gives me an idea for a slow function: `mad(runif(1000)`

``` r
system.time(
    result <- lapply(1:num_interations, function(i) {
        mad(runif(1000))
    })
)
```

    ##    user  system elapsed 
    ##  11.232   0.144  11.374

# Do something slow using mclapply

Running using Rstudio and 4 cores, I do see a successful speedup

``` r
system.time(
    result <- mclapply(1:num_interations, function(i) {
        mad(runif(1000))
    }, mc.cores=num_threads)
)
```

    ##    user  system elapsed 
    ##  11.578   0.328   3.388

# Finished

``` r
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Etc/UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] digest_0.6.37     fastmap_1.2.0     xfun_0.60         knitr_1.50       
    ##  [5] htmltools_0.5.8.1 rmarkdown_2.30    lifecycle_1.0.5   cli_3.6.5        
    ##  [9] textshaping_1.0.4 systemfonts_1.3.1 compiler_4.5.2    rstudioapi_0.17.1
    ## [13] tools_4.5.2       ragg_1.5.0        evaluate_1.0.5    yaml_2.3.10      
    ## [17] rlang_1.3.0
