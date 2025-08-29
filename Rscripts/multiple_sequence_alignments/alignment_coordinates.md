alignment_coordinates
================
Janet Young

2025-08-29

# Goal

given coordinates of one or more features in a reference sequence,
figure out what the coordinates are in an alignment that contains the
reference sequence

# Load libraries

``` r
knitr::opts_chunk$set(echo = TRUE)

library(tidyverse)
library(here)

## Biostrings is a bioconductor package, not a regular R package. This page has instructions on how to install it:
## https://bioconductor.org/packages/release/bioc/html/Biostrings.html
library(Biostrings)
```

# Read example alignment

``` r
## readDNAStringSet is a Biostrings function that reads a multi-sequence fasta file
## readAAStringSet is the amino acid equivalent
aln <- readDNAStringSet("example_alignment_files/cenH3_aln9.nt.fa")

# get rid of descriptions in the sequence names
names(aln) <- sapply( strsplit(names(aln), " "), "[[", 1)
```

What does it look like?

34 sequences in the alignment (rows) and 492 aligned positions (columns)

``` r
aln
```

    ## DNAStringSet object of length 34:
    ##      width seq                                              names               
    ##  [1]   492 ATG--------------------...CTTCGAGGGCGGACTCCCCTAA mouse_Cenpa_ORF
    ##  [2]   492 ATGGTCGGGCGCCGCAAGCCAGG...CATTGAGGGAGGACTCGGCTGA rat_Cenpa_ORF
    ##  [3]   492 ATG--------------------...CATTGAGGGAGGACTTGGCTGA Chinese_hamster_C...
    ##  [4]   492 ATG--------------------...CATTCAGGAAGGGCTTGGCTGA leopard_CENPA_ORF
    ##  [5]   492 ATG--------------------...CATTCAGGAAGGGCTTGGCTGA panda_cenH3_ORF_A...
    ##  ...   ... ...
    ## [30]   492 ATG--------------------...CATTGAGGCGGGACTCGGCTGA crab-eating_macaq...
    ## [31]   492 ATG--------------------...CATTGAGGCGGGACTCGGCTGA Rhesus_monkey_CEN...
    ## [32]   492 ATG--------------------...CATTGAGGCGGGACTCGGCTGA gelada_CENPA_ORF
    ## [33]   492 ATG--------------------...CATTGAGGCGGGACTCGGCTGA olive_baboon_CENP...
    ## [34]   492 ATG--------------------...CATTGAGGCGGGACTCGGCTGA green_monkey_CENP...

Figure out the name of our reference sequence (human) - the sequence
which we want coordinates to refer to:

``` r
ref_name <- grep("human_", names(aln), value=TRUE)
ref_name
```

    ## [1] "human_CENPA_ORF"

Set up a couple of functions

`getUngappedPosOneSeq` takes a single gapped sequence (i.e. a sequence
in the alignment), and gets a tibble of position in each sequence versus
position in alignment. We use cumulative sum of non-gap bases.

`makeLookupTibble` is a bigger function that applies the
`getUngappedPosOneSeq` to every sequence in an alignment and returns a
tibble

``` r
###### getUngappedPosOneSeq is a function to take a single aligned sequence, and to return a number for each position in the alignment that represents the position in the original sequence
## myGappedSeq is a DNAString (or AAString/BString, etc)
getUngappedPosOneSeq <- function(myGappedSeq) {
    # split the aligned sequence into individual characters
    mySeq <- strsplit(as.character(myGappedSeq),"")[[1]]
    # as we go along the aligned sequence from start to finish, count how many non-gap characters we have encountered
    myCounts <- cumsum(mySeq != "-")
    # gaps in the aligned sequence should have no coordinate - change those counts to NA
    myCounts[which(mySeq=="-")] <- NA
    return(myCounts)
}

##### makeLookupTibble - a function to take an alignment, and for every single sequence, get a position lookup table, and join the together
makeLookupTibble <- function(alignment) {
    output <- tibble(aln_pos=1:width(alignment)[1])
    each_seq_lookup <- sapply( names(alignment), function(each_seq_name) {
        getUngappedPosOneSeq( alignment[[each_seq_name]] )
    } ) 
    output <- bind_cols(output, each_seq_lookup)
    return(output)
}
```

Make a position lookup tibble for our alignment, just the reference
sequence. Show the last few rows of the tibble.

``` r
alnPos_lookup_table <- makeLookupTibble( aln[ref_name] )

alnPos_lookup_table %>% 
    tail()
```

    ## # A tibble: 6 × 2
    ##   aln_pos human_CENPA_ORF
    ##     <int>           <int>
    ## 1     487             418
    ## 2     488             419
    ## 3     489             420
    ## 4     490             421
    ## 5     491             422
    ## 6     492             423

# Tabulate every position in the alignment

The `consensusMatrix()` function (Biostrings) can tabulate what’s
present at each alignment site (column):

It gives us a matrix object with a row for each nucleotide (class) it
counted, and a column for each alignment position. Here are counts for
the first 10 positions

``` r
## as a reminder, our alignment has 34 seqs and 492 positions
aln_counts <- consensusMatrix(aln, baseOnly=TRUE)
aln_counts[,1:10]
```

    ##       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
    ## A       34    0    0    0    0    0    0    0    0     0
    ## C        0    0    0    0    0    1    0    0    0     1
    ## G        0    0   34    1    0    0    1    1    1     0
    ## T        0   34    0    0    1    0    0    0    0     0
    ## other    0    0    0   33   33   33   33   33   33    33

Maybe we want to use the lookup table we made earlier to add human
positions on to that alignment

First, let’s turn that matrix 90 degrees and make it into a tibble,
which will be easier to work with

``` r
aln_counts <- aln_counts %>% 
    t() %>% 
    as.data.frame() %>% 
    as_tibble(rownames="aln_pos") %>% 
    mutate(aln_pos=as.integer(aln_pos))

aln_counts %>% 
    head()
```

    ## # A tibble: 6 × 6
    ##   aln_pos     A     C     G     T other
    ##     <int> <int> <int> <int> <int> <int>
    ## 1       1    34     0     0     0     0
    ## 2       2     0     0     0    34     0
    ## 3       3     0     0    34     0     0
    ## 4       4     0     0     1     0    33
    ## 5       5     0     0     0     1    33
    ## 6       6     0     1     0     0    33

Use `left_join()` to add the human positions

``` r
aln_counts <- left_join(aln_counts,
                        alnPos_lookup_table, 
                        by="aln_pos") %>% 
    relocate(human_CENPA_ORF, .after=aln_pos) %>% 
    dplyr::rename(human_pos=human_CENPA_ORF)
aln_counts
```

    ## # A tibble: 492 × 7
    ##    aln_pos human_pos     A     C     G     T other
    ##      <int>     <int> <int> <int> <int> <int> <int>
    ##  1       1         1    34     0     0     0     0
    ##  2       2         2     0     0     0    34     0
    ##  3       3         3     0     0    34     0     0
    ##  4       4        NA     0     0     1     0    33
    ##  5       5        NA     0     0     0     1    33
    ##  6       6        NA     0     1     0     0    33
    ##  7       7        NA     0     0     1     0    33
    ##  8       8        NA     0     0     1     0    33
    ##  9       9        NA     0     0     1     0    33
    ## 10      10        NA     0     1     0     0    33
    ## # ℹ 482 more rows

Add the nucleotide in human

``` r
human_nucs <- aln[[ref_name]] %>% 
    as.character() %>% 
    str_split("")
human_nucs <- human_nucs[[1]]

aln_counts$human_nuc <- human_nucs

aln_counts <- aln_counts %>% 
    relocate(human_nuc, .after=human_pos)
aln_counts
```

    ## # A tibble: 492 × 8
    ##    aln_pos human_pos human_nuc     A     C     G     T other
    ##      <int>     <int> <chr>     <int> <int> <int> <int> <int>
    ##  1       1         1 A            34     0     0     0     0
    ##  2       2         2 T             0     0     0    34     0
    ##  3       3         3 G             0     0    34     0     0
    ##  4       4        NA -             0     0     1     0    33
    ##  5       5        NA -             0     0     0     1    33
    ##  6       6        NA -             0     1     0     0    33
    ##  7       7        NA -             0     0     1     0    33
    ##  8       8        NA -             0     0     1     0    33
    ##  9       9        NA -             0     0     1     0    33
    ## 10      10        NA -             0     1     0     0    33
    ## # ℹ 482 more rows

Perhaps we drop the positions where human has a gap - we might not care
about those

``` r
aln_counts <- aln_counts %>% 
    filter(human_nuc != "-")
aln_counts
```

    ## # A tibble: 423 × 8
    ##    aln_pos human_pos human_nuc     A     C     G     T other
    ##      <int>     <int> <chr>     <int> <int> <int> <int> <int>
    ##  1       1         1 A            34     0     0     0     0
    ##  2       2         2 T             0     0     0    34     0
    ##  3       3         3 G             0     0    34     0     0
    ##  4      25         4 G             1     0    33     0     0
    ##  5      26         5 G             0     2    32     0     0
    ##  6      27         6 C             1    33     0     0     0
    ##  7      28         7 C             0    34     0     0     0
    ##  8      29         8 C             0    34     0     0     0
    ##  9      30         9 G             0     0    34     0     0
    ## 10      31        10 C             1    33     0     0     0
    ## # ℹ 413 more rows

# Look at regions of interest

Generate example features, with coordinates in the reference sequence

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

Set up another function called `addAlnCoords` that will use the lookup
table, and the feature table, and will convert coordinates in the
reference sequence into coordinates in the alignment

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
        ## the rename column is weird and misbehaves when certain Bioconductor packages are loaded, unless I specify that I want to use the rename function from the dplyr package, like this:
        dplyr::rename(start_aln=aln_pos)
    # look up end aln_pos
    feature_tbl <- left_join(feature_tbl,
                             lookup_tbl,
                             by=c("end"="ref_pos") ) %>% 
        dplyr::rename(end_aln=aln_pos)
    return(feature_tbl)
}
```

Use `addAlnCoords` with our lookup table to convert coordinates in the
example features tibble (`human_features`)

``` r
human_features <- addAlnCoords( feature_tbl=human_features, 
                                lookup_tbl=alnPos_lookup_table, 
                                refseq_name=ref_name )
human_features
```

    ## # A tibble: 3 × 6
    ##   name    start   end strand start_aln end_aln
    ##   <chr>   <dbl> <dbl> <chr>      <int>   <int>
    ## 1 region1     1    90 +              1     156
    ## 2 region2   101   140 +            167     209
    ## 3 region3   331   400 +            400     469

## Extract aligned regions

Get those three regions from the alignment - we use the `narrow`
function

``` r
human_feature_alns <- lapply(1:nrow(human_features), 
                             function(i) {
    narrow(aln, 
           start=human_features$start_aln[i],
           end=human_features$end_aln[i] )
})
human_feature_alns
```

    ## [[1]]
    ## DNAStringSet object of length 34:
    ##      width seq                                              names               
    ##  [1]   156 ATG--------------------...----GGACCCTCGCGACAGAGC mouse_Cenpa_ORF
    ##  [2]   156 ATGGTCGGGCGCCGCAAGCCAGG...----GGACCCTCGCGACGGAGC rat_Cenpa_ORF
    ##  [3]   156 ATG--------------------...----GGACCCTCGCGACGCAGC Chinese_hamster_C...
    ##  [4]   156 ATG--------------------...----AGCCCCCCCCGGCGGGGC leopard_CENPA_ORF
    ##  [5]   156 ATG--------------------...----AGCCCTCCCCGGCGGGGC panda_cenH3_ORF_A...
    ##  ...   ... ...
    ## [30]   156 ATG--------------------...----GGCCCCTCCCGGCGGGGC crab-eating_macaq...
    ## [31]   156 ATG--------------------...----GGCCCCTCCCGGCGGGGC Rhesus_monkey_CEN...
    ## [32]   156 ATG--------------------...----GGCCCCTCCCGGCGGGGC gelada_CENPA_ORF
    ## [33]   156 ATG--------------------...----GGCCCCTCCCGGCGGGGC olive_baboon_CENP...
    ## [34]   156 ATG--------------------...----GGCCCCTCCCGGCGGGGC green_monkey_CENP...
    ## 
    ## [[2]]
    ## DNAStringSet object of length 34:
    ##      width seq                                              names               
    ##  [1]    43 GCTCTCAGACACTGCGCAGAAGA---CAGAAATTC---ATGTG      mouse_Cenpa_ORF
    ##  [2]    43 GCCCTCAGGCACTACACAGAAGA---CGGAGATTC---CTGTG      rat_Cenpa_ORF
    ##  [3]    43 GTAAGAGG------------------CGGAAATTT---CTGTG      Chinese_hamster_C...
    ##  [4]    43 GCGCTTCCTCCCGTCAGCGTGGTCCCCGGAGAAGT---CGGGT      leopard_CENPA_ORF
    ##  [5]    43 GCGCTCCCTCCCGTCAGCGTCGTCCCCGGAGACAT---CGGGT      panda_cenH3_ORF_A...
    ##  ...   ... ...
    ## [30]    43 GCCCTTCCTCCCGTCAACATGGTCGGCGGAGACAA---GCTTG      crab-eating_macaq...
    ## [31]    43 GCCCTTCCTCCCGTCAACATGGTCGGCGGAGACAA---GCTTG      Rhesus_monkey_CEN...
    ## [32]    43 GCCCTTCCTCCCGTCAACATGGTCGGCGGAGACAA---GCTTG      gelada_CENPA_ORF
    ## [33]    43 GCCCTTCCTCCCGTCAACATGGTCGGCGGAGACAA---GCTTG      olive_baboon_CENP...
    ## [34]    43 GCCCTTCCTCCCATCAACATGGTCGGCGGAGACAA---GCTTG      green_monkey_CENP...
    ## 
    ## [[3]]
    ## DNAStringSet object of length 34:
    ##      width seq                                              names               
    ##  [1]    70 CTCCTCTCCTTACATGCTGGTCG...CAGTTGACCAGGAGAATCCGAG mouse_Cenpa_ORF
    ##  [2]    70 CTCCTCTCCTTACATGCTGGCCG...CAGCTAGCCAGGAGGATCCGAG rat_Cenpa_ORF
    ##  [3]    70 CTCCTCACCTTACATGCCGGCAG...CAGCTCACCAGGAGGATCCGAG Chinese_hamster_C...
    ##  [4]    70 CTTCTCTCCTTACATGCCGGCCG...CAGCTGGCCAGGAGGATCCGAG leopard_CENPA_ORF
    ##  [5]    70 CTCCTGTCCTTACATGCCGGCCG...CAGTTGGCCAGGAGGATCCGAG panda_cenH3_ORF_A...
    ##  ...   ... ...
    ## [30]    70 CTCCTCGCCTTACATGCCGGCCG...CAACTGGCCCGGAGGATCCGGG crab-eating_macaq...
    ## [31]    70 CTCCTCGCCTTACATGCCGGCCG...CAACTGGCCCGGAGGATCCGGG Rhesus_monkey_CEN...
    ## [32]    70 CTCCTCACCTTACATGCCGGCCG...CAACTGGCCCGGAGGATCCGGG gelada_CENPA_ORF
    ## [33]    70 CTCCTCACCTTACATGCCGGCCG...CAACTGGCCCGGAGGATCCGGG olive_baboon_CENP...
    ## [34]    70 CTCCTCACCTTACATGCCGGCCG...CAACTGGCCCGGAGGATCCGGG green_monkey_CENP...

We can use the `xscat` function to join those alignment pieces together.

Here’s the simplest way to use `xscat`.

``` r
xscat(human_feature_alns[[1]],
      human_feature_alns[[2]],
      human_feature_alns[[3]])
```

    ## DNAStringSet object of length 34:
    ##      width seq
    ##  [1]   269 ATG---------------------GGCCCGCGT...CCAAAGACATTCAGTTGACCAGGAGAATCCGAG
    ##  [2]   269 ATGGTCGGGCGCCGCAAGCCAGGGACCCCGAGG...CCAAAGATGTGCAGCTAGCCAGGAGGATCCGAG
    ##  [3]   269 ATG---------------------GGCCCGCGC...CCAAGGACATACAGCTCACCAGGAGGATCCGAG
    ##  [4]   269 ATG---------------------GGCCCGCGC...CGAAGGATGTGCAGCTGGCCAGGAGGATCCGAG
    ##  [5]   269 ATG---------------------GGCCCGCGC...CGAAGGATGTGCAGTTGGCCAGGAGGATCCGAG
    ##  ...   ... ...
    ## [30]   269 ATG---------------------GGCCCGCGC...CAAAGGATGTGCAACTGGCCCGGAGGATCCGGG
    ## [31]   269 ATG---------------------GGCCCGCGC...CAAAGGATGTGCAACTGGCCCGGAGGATCCGGG
    ## [32]   269 ATG---------------------GGCCCGCGC...CAAAGGATGTGCAACTGGCCCGGAGGATCCGGG
    ## [33]   269 ATG---------------------GGCCCGCGC...CAAAGGATGTGCAACTGGCCCGGAGGATCCGGG
    ## [34]   269 ATG---------------------GGCCCGCGC...CAAAGGATGTGCAACTGGCCCGGAGGATCCGGG

However, using `xscat` that way is not very generalizable, e.g. if we
have a large list, or if we don’t know ahead of time how long our list
is:

A better way, perhaps is to wrap `xscat` inside the `do.call` function -
it is a function that I think of as “unpacking” a list and running a
function on all the elements together:

``` r
do.call("xscat", human_feature_alns)
```

    ## DNAStringSet object of length 34:
    ##      width seq
    ##  [1]   269 ATG---------------------GGCCCGCGT...CCAAAGACATTCAGTTGACCAGGAGAATCCGAG
    ##  [2]   269 ATGGTCGGGCGCCGCAAGCCAGGGACCCCGAGG...CCAAAGATGTGCAGCTAGCCAGGAGGATCCGAG
    ##  [3]   269 ATG---------------------GGCCCGCGC...CCAAGGACATACAGCTCACCAGGAGGATCCGAG
    ##  [4]   269 ATG---------------------GGCCCGCGC...CGAAGGATGTGCAGCTGGCCAGGAGGATCCGAG
    ##  [5]   269 ATG---------------------GGCCCGCGC...CGAAGGATGTGCAGTTGGCCAGGAGGATCCGAG
    ##  ...   ... ...
    ## [30]   269 ATG---------------------GGCCCGCGC...CAAAGGATGTGCAACTGGCCCGGAGGATCCGGG
    ## [31]   269 ATG---------------------GGCCCGCGC...CAAAGGATGTGCAACTGGCCCGGAGGATCCGGG
    ## [32]   269 ATG---------------------GGCCCGCGC...CAAAGGATGTGCAACTGGCCCGGAGGATCCGGG
    ## [33]   269 ATG---------------------GGCCCGCGC...CAAAGGATGTGCAACTGGCCCGGAGGATCCGGG
    ## [34]   269 ATG---------------------GGCCCGCGC...CAAAGGATGTGCAACTGGCCCGGAGGATCCGGG

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
    ##  [1] Biostrings_2.76.0   GenomeInfoDb_1.44.1 XVector_0.48.0     
    ##  [4] IRanges_2.42.0      S4Vectors_0.46.0    BiocGenerics_0.54.0
    ##  [7] generics_0.1.4      here_1.0.1          lubridate_1.9.4    
    ## [10] forcats_1.0.0       stringr_1.5.1       dplyr_1.1.4        
    ## [13] purrr_1.1.0         readr_2.1.5         tidyr_1.3.1        
    ## [16] tibble_3.3.0        ggplot2_3.5.2       tidyverse_2.0.0    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] utf8_1.2.6              stringi_1.8.7           hms_1.1.3              
    ##  [4] digest_0.6.37           magrittr_2.0.3          evaluate_1.0.4         
    ##  [7] grid_4.5.1              timechange_0.3.0        RColorBrewer_1.1-3     
    ## [10] fastmap_1.2.0           rprojroot_2.1.0         jsonlite_2.0.0         
    ## [13] httr_1.4.7              UCSC.utils_1.4.0        scales_1.4.0           
    ## [16] cli_3.6.5               crayon_1.5.3            rlang_1.1.6            
    ## [19] withr_3.0.2             yaml_2.3.10             tools_4.5.1            
    ## [22] tzdb_0.5.0              GenomeInfoDbData_1.2.14 vctrs_0.6.5            
    ## [25] R6_2.6.1                lifecycle_1.0.4         pkgconfig_2.0.3        
    ## [28] pillar_1.11.0           gtable_0.3.6            glue_1.8.0             
    ## [31] xfun_0.52               tidyselect_1.2.1        rstudioapi_0.17.1      
    ## [34] knitr_1.50              farver_2.1.2            htmltools_0.5.8.1      
    ## [37] rmarkdown_2.29          compiler_4.5.1
