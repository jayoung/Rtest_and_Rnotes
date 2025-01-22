alignment_coordinates
================
Janet Young

2025-01-21

Goal - given coordinates of one or more features in a reference
sequence, figure out what the coordinates are in an alignment that
contains the reference sequence

# Load libraries

``` r
knitr::opts_chunk$set(echo = TRUE)
library(Biostrings)
library(tidyverse)
library(here)
```

# Read example alignment:

``` r
aln <- readDNAStringSet("example_alignment_files/cenH3_aln8.nt.fa")
# get rid of descriptions
names(aln) <- sapply( strsplit(names(aln), " "), "[[", 1)
```

Define which sequence is our reference - the sequence to which
coordinates refer

``` r
ref_name <- grep("human_", names(aln), value=TRUE)
ref_name
```

    ## [1] "human_CENPA_ORF"

Make up some example feature coordinates in the reference sequence

``` r
human_features <- tibble(
    name=c("region1","region2","region3"),
    start=c(1, 101, 331),
    end=c(90,140,400)
) 
human_features
```

    ## # A tibble: 3 × 3
    ##   name    start   end
    ##   <chr>   <dbl> <dbl>
    ## 1 region1     1    90
    ## 2 region2   101   140
    ## 3 region3   331   400

Set up a couple of functions

`getUngappedPosOneSeq` takes a single gapped sequence (i.e. a sequence
in the alignment), and gets a tibble of position in each sequence versus
position in alignment. We use cumulative sum of non-gap bases.

`makeLookupTibble` is a bigger function that applies the
`getUngappedPosOneSeq` to every sequence in an alignment and returns a
tibble

``` r
## in getUngappedPosOneSeq, myGappedSeq is a DNAString (or AAString/BString, etc)
getUngappedPosOneSeq <- function(myGappedSeq) {
    mySeq <- strsplit(as.character(myGappedSeq),"")[[1]]
    myCounts <- cumsum(mySeq != "-")
    myCounts[which(mySeq=="-")] <- NA
    return(myCounts)
}

makeLookupTibble <- function(alignment) {
    output <- tibble(aln_pos=1:width(alignment)[1])
    each_seq_lookup <- sapply( names(alignment), function(each_seq_name) {
        getUngappedPosOneSeq( alignment[[each_seq_name]] )
    } ) 
    output <- bind_cols(output, each_seq_lookup)
    return(output)
}
```

Make a position lookup tibble for our alignment, considering just the
reference sequence that our features are defined in. Show the last few
rows of the tibble.

``` r
alnPos_lookup_table <- makeLookupTibble( aln[ref_name] )

alnPos_lookup_table %>% 
    tail()
```

    ## # A tibble: 6 × 2
    ##   aln_pos human_CENPA_ORF
    ##     <int>           <int>
    ## 1     598             418
    ## 2     599             419
    ## 3     600             420
    ## 4     601             421
    ## 5     602             422
    ## 6     603             423

Set up another function called `addAlnCoords` that will use the lookup
table to convert start and end coordinates in the reference sequence to
start and end coordinates in the alignment

``` r
addAlnCoords <- function(feature_tbl,
                         lookup_tbl, 
                         refseq_name) {
    ## check the ref sequence is present in the lookup tbl
    if(! refseq_name %in% colnames(lookup_tbl)) {
        stop("\n\nERROR - could not find refseq called ",
             refseq_name, "in the lookup_tbl colnames\n\n")
    }
    
    ## get lookup tbl in a useful format
    lookup_tbl <- lookup_tbl %>% 
        select(aln_pos, ref_pos=matches(refseq_name))
    
    # look up start aln_pos
    feature_tbl <- left_join(feature_tbl,
                             lookup_tbl,
                             by=c("start"="ref_pos") ) %>% 
        rename(start_aln=aln_pos)
    # look up end aln_pos
    feature_tbl <- left_join(feature_tbl,
                             lookup_tbl,
                             by=c("end"="ref_pos") ) %>% 
        rename(end_aln=aln_pos)
    return(feature_tbl)
}
```

Use `addAlnCoords` with our lookup table to convert coordinates in the
example features tibble (`human_features`)

``` r
addAlnCoords( feature_tbl=human_features, 
              lookup_tbl=alnPos_lookup_table, 
              refseq_name=ref_name )
```

    ## # A tibble: 3 × 5
    ##   name    start   end start_aln end_aln
    ##   <chr>   <dbl> <dbl>     <int>   <int>
    ## 1 region1     1    90       106     267
    ## 2 region2   101   140       278     320
    ## 3 region3   331   400       511     580

# Finished

``` r
sessionInfo()
```

    ## R version 4.4.2 (2024-10-31)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Sequoia 15.2
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
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
    ##  [1] here_1.0.1          lubridate_1.9.4     forcats_1.0.0      
    ##  [4] stringr_1.5.1       dplyr_1.1.4         purrr_1.0.2        
    ##  [7] readr_2.1.5         tidyr_1.3.1         tibble_3.2.1       
    ## [10] ggplot2_3.5.1       tidyverse_2.0.0     Biostrings_2.74.1  
    ## [13] GenomeInfoDb_1.42.1 XVector_0.46.0      IRanges_2.40.1     
    ## [16] S4Vectors_0.44.0    BiocGenerics_0.52.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] utf8_1.2.4              generics_0.1.3          stringi_1.8.4          
    ##  [4] hms_1.1.3               digest_0.6.37           magrittr_2.0.3         
    ##  [7] evaluate_1.0.3          grid_4.4.2              timechange_0.3.0       
    ## [10] fastmap_1.2.0           rprojroot_2.0.4         jsonlite_1.8.9         
    ## [13] httr_1.4.7              UCSC.utils_1.2.0        scales_1.3.0           
    ## [16] cli_3.6.3               rlang_1.1.4             crayon_1.5.3           
    ## [19] munsell_0.5.1           withr_3.0.2             yaml_2.3.10            
    ## [22] tools_4.4.2             tzdb_0.4.0              colorspace_2.1-1       
    ## [25] GenomeInfoDbData_1.2.13 vctrs_0.6.5             R6_2.5.1               
    ## [28] lifecycle_1.0.4         zlibbioc_1.52.0         pkgconfig_2.0.3        
    ## [31] pillar_1.10.1           gtable_0.3.6            glue_1.8.0             
    ## [34] xfun_0.50               tidyselect_1.2.1        rstudioapi_0.17.1      
    ## [37] knitr_1.49              htmltools_0.5.8.1       rmarkdown_2.29         
    ## [40] compiler_4.4.2
