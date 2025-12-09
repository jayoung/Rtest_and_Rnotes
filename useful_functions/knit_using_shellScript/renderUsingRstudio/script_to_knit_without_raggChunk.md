script_to_knit
================
Janet Young

2025-12-09

If I “knit” from Rstudio, plots are saved as png. I like that (files are
smaller so better to host on github)

If I “render” from the command line (on the cluster) using Rscript,
plots would be saved as svg (sometimes too big for github).

The fix (thanks to Dan Tenenbaum) is to include a chunk of code in in my
Rmd script - that code now lives in `~/.ragg_png_functions_from_dan.R`
and is sourced when R starts up, via `~/.Rprofile`. For projects where I
use `renv`, I would also want to source it in the LOCAL `.Rprofile` file
for that project.

The shell script to render looks like this:

    #!/bin/bash
    source /app/lmod/lmod/init/profile
    module purge
    module load fhR/4.4.1-foss-2023b-R-4.4.1
    module load Pandoc/2.13

    Rscript -e 'rmarkdown::render("script_to_knit.Rmd", output_format="github_document", clean=TRUE)' > script_to_knit.Rrender.Rout.txt 2> script_to_knit.Rrender.Rerr.txt

``` r
plot(mtcars$mpg, mtcars$disp)
```

![](script_to_knit_without_raggChunk_files/figure-gfm/mtcars%20plot-1.png)<!-- -->

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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_4.4.0    fastmap_1.2.0     cli_3.6.3         ragg_1.3.1       
    ##  [5] htmltools_0.5.8.1 tools_4.4.0       rstudioapi_0.16.0 yaml_2.3.8       
    ##  [9] rmarkdown_2.26    knitr_1.50        digest_0.6.35     xfun_0.53        
    ## [13] textshaping_0.3.7 lifecycle_1.0.4   systemfonts_1.1.0 rlang_1.1.5      
    ## [17] evaluate_0.23
