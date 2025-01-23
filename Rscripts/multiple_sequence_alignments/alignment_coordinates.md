alignment_coordinates
================
Janet Young

2025-01-22

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
    end=c(90,140,400),
    strand="+"
) 
human_features
```

    ## # A tibble: 3 × 4
    ##   name    start   end strand
    ##   <chr>   <dbl> <dbl> <chr> 
    ## 1 region1     1    90 +     
    ## 2 region2   101   140 +     
    ## 3 region3   331   400 +

Set up a couple of functions

`getUngappedPosOneSeq` takes a single gapped sequence (i.e. a sequence
in the alignment), and gets a tibble of position in each sequence versus
position in alignment. We use cumulative sum of non-gap bases.

`makeLookupTibble` is a bigger function that applies the
`getUngappedPosOneSeq` to every sequence in an alignment and returns a
tibble

``` r
myGappedSeq <- aln[[ref_name]]

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
human_features <- addAlnCoords( feature_tbl=human_features, 
                                lookup_tbl=alnPos_lookup_table, 
                                refseq_name=ref_name )
```

Get those three regions from the alignment

``` r
human_feature_alns <- lapply(1:nrow(human_features), 
                             function(i) {
    narrow(aln, 
           start=human_features$start_aln[i],
           end=human_features$end_aln[i] )
})
```

We can use the `xscat` function to join those alignment pieces together.

Here’s the simplest way to use `xscat`.

``` r
xscat(human_feature_alns[[1]],
      human_feature_alns[[2]],
      human_feature_alns[[3]])
```

    ## DNAStringSet object of length 35:
    ##      width seq
    ##  [1]   275 ATG---------------------GGCCCGCGT...CCAAAGACATTCAGTTGACCAGGAGAATCCGAG
    ##  [2]   275 ATGGTCGGGCGCCGCAAGCCAGGGACCCCGAGG...CCAAAGATGTGCAGCTAGCCAGGAGGATCCGAG
    ##  [3]   275 ATG---------------------GGCCCGCGC...CCAAGGACATACAGCTCACCAGGAGGATCCGAG
    ##  [4]   275 ATG---------------------GGCCCGCGC...CGAAGGATGTGCAGCTGGCCAGGAGGATCCGAG
    ##  [5]   275 ATG---------------------GGCCCGCGC...CGAAGGATGTGCAGTTGGCCAGGAGGATCCGAG
    ##  ...   ... ...
    ## [31]   275 ATG---------------------GGCCCGCGC...CAAAGGATGTGCAACTGGCCCGGAGGATCCGGG
    ## [32]   275 ATG---------------------GGCCCGCGC...CAAAGGATGTGCAACTGGCCCGGAGGATCCGGG
    ## [33]   275 ATG---------------------GGCCCGCGC...CAAAGGATGTGCAACTGGCCCGGAGGATCCGGG
    ## [34]   275 ATG---------------------GGCCCGCGC...CAAAGGATGTGCAACTGGCCCGGAGGATCCGGG
    ## [35]   275 ATG---------------------GGCCCGCGC...CAAAGGATGTGCAACTGGCCCGGAGGATCCGGG

However, using `xscat` that way is not very generalizable, e.g. if we
have a large list, or if we don’t know head of time how long our list
is:

A better way, perhaps is to wrap `xscat` inside the `do.call` function -
it is a function that I think of as “unpacking” a list and running a
function on all the elements together:

``` r
do.call("xscat", human_feature_alns)
```

    ## DNAStringSet object of length 35:
    ##      width seq
    ##  [1]   275 ATG---------------------GGCCCGCGT...CCAAAGACATTCAGTTGACCAGGAGAATCCGAG
    ##  [2]   275 ATGGTCGGGCGCCGCAAGCCAGGGACCCCGAGG...CCAAAGATGTGCAGCTAGCCAGGAGGATCCGAG
    ##  [3]   275 ATG---------------------GGCCCGCGC...CCAAGGACATACAGCTCACCAGGAGGATCCGAG
    ##  [4]   275 ATG---------------------GGCCCGCGC...CGAAGGATGTGCAGCTGGCCAGGAGGATCCGAG
    ##  [5]   275 ATG---------------------GGCCCGCGC...CGAAGGATGTGCAGTTGGCCAGGAGGATCCGAG
    ##  ...   ... ...
    ## [31]   275 ATG---------------------GGCCCGCGC...CAAAGGATGTGCAACTGGCCCGGAGGATCCGGG
    ## [32]   275 ATG---------------------GGCCCGCGC...CAAAGGATGTGCAACTGGCCCGGAGGATCCGGG
    ## [33]   275 ATG---------------------GGCCCGCGC...CAAAGGATGTGCAACTGGCCCGGAGGATCCGGG
    ## [34]   275 ATG---------------------GGCCCGCGC...CAAAGGATGTGCAACTGGCCCGGAGGATCCGGG
    ## [35]   275 ATG---------------------GGCCCGCGC...CAAAGGATGTGCAACTGGCCCGGAGGATCCGGG

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
    ##  [1] here_1.0.1          lubridate_1.9.3     forcats_1.0.0      
    ##  [4] stringr_1.5.1       dplyr_1.1.4         purrr_1.0.2        
    ##  [7] readr_2.1.5         tidyr_1.3.1         tibble_3.2.1       
    ## [10] ggplot2_3.5.1       tidyverse_2.0.0     Biostrings_2.72.0  
    ## [13] GenomeInfoDb_1.40.0 XVector_0.44.0      IRanges_2.38.0     
    ## [16] S4Vectors_0.42.0    BiocGenerics_0.50.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] utf8_1.2.4              generics_0.1.3          stringi_1.8.4          
    ##  [4] hms_1.1.3               digest_0.6.35           magrittr_2.0.3         
    ##  [7] timechange_0.3.0        evaluate_0.23           grid_4.4.0             
    ## [10] fastmap_1.2.0           rprojroot_2.0.4         jsonlite_1.8.8         
    ## [13] httr_1.4.7              fansi_1.0.6             UCSC.utils_1.0.0       
    ## [16] scales_1.3.0            cli_3.6.2               rlang_1.1.4            
    ## [19] crayon_1.5.2            munsell_0.5.1           withr_3.0.0            
    ## [22] yaml_2.3.8              tools_4.4.0             tzdb_0.4.0             
    ## [25] colorspace_2.1-0        GenomeInfoDbData_1.2.12 vctrs_0.6.5            
    ## [28] R6_2.5.1                lifecycle_1.0.4         zlibbioc_1.50.0        
    ## [31] pkgconfig_2.0.3         pillar_1.9.0            gtable_0.3.5           
    ## [34] glue_1.7.0              xfun_0.44               tidyselect_1.2.1       
    ## [37] rstudioapi_0.16.0       knitr_1.46              htmltools_0.5.8.1      
    ## [40] rmarkdown_2.26          compiler_4.4.0
