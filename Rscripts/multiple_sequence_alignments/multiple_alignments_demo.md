multiple_alignments_demo
================
Janet Young

2026-09-03

# Goal

Demonstrate various code tricks for multiple alignments

1.  Coordinates. Given coordinates of one or more features in a
    reference sequence, figure out what the coordinates are in an
    alignment that contains the reference sequence

2.  Getting alignment slices and concatenating them

# Load libraries

``` r
knitr::opts_chunk$set(echo = TRUE)

library(tidyverse)
library(here)

## Biostrings and GenomicRanges are Bioconductor packages, not regular R packages. This page has instructions on how to install them:
## https://bioconductor.org/packages/release/bioc/html/Biostrings.html
library(Biostrings)
library(GenomicRanges)
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

`getAlnPosLookupTable` is a bigger function that applies the
`getUngappedPosOneSeq` to every sequence in an alignment and returns a
tibble.

A combined version of those functions is also called
`getAlnPosLookupTable` and can be found in the file
`useful_functions/multiple_sequence_alignments_functions.R`, which I
often read into other scripts so I don’t have to redefine the function
each time.

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

##### getAlnPosLookupTable - a function to take an alignment, and for every single sequence, get a position lookup table, and join the together
getAlnPosLookupTable <- function(alignment) {
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
alnPos_lookup_table <- getAlnPosLookupTable( aln[ref_name] )

alnPos_lookup_table |> 
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
aln_counts <- aln_counts |> 
    t() |> 
    as.data.frame() |> 
    as_tibble(rownames="aln_pos") |> 
    mutate(aln_pos=as.integer(aln_pos))

aln_counts |> 
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
                        by="aln_pos") |> 
    relocate(human_CENPA_ORF, .after=aln_pos) |> 
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
human_nucs <- aln[[ref_name]] |> 
    as.character() |> 
    str_split("")
human_nucs <- human_nucs[[1]]

aln_counts$human_nuc <- human_nucs

aln_counts <- aln_counts |> 
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
aln_counts <- aln_counts |> 
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
    lookup_tbl <- lookup_tbl |> 
        select(aln_pos, ref_pos=matches(refseq_name))
    
    # look up start aln_pos
    feature_tbl <- left_join(feature_tbl,
                             lookup_tbl,
                             by=c("start"="ref_pos") ) |> 
        ## the rename column is weird and misbehaves when certain Bioconductor packages are loaded, unless I specify that I want to use the rename function from the dplyr package, like this:
        dplyr::rename(start_aln=aln_pos)
    # look up end aln_pos
    feature_tbl <- left_join(feature_tbl,
                             lookup_tbl,
                             by=c("end"="ref_pos") ) |> 
        dplyr::rename(end_aln=aln_pos)
    return(feature_tbl)
}
```

Use `addAlnCoords` with our lookup table to convert coordinates in the
example features tibble (`human_features`), adding them as extra
columns.

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

More often I store the features as a GRanges object, in which case I can
use the `convertCoordsOneSeq` function defined in the
`useful_functions/multiple_sequence_alignments_functions.R` file, which
I often read in to otehr scripts.

Here’s how we convert the original human_features object to GRanges -
get a version BEFORE I added alignment coordinates, so I can demo how we
convert them using that other function.

``` r
human_features_gr <- human_features |> 
    select(start, end, strand, "name") |> 
    mutate(seqnames="human_CENPA_ORF") |> 
    GRanges()
human_features_gr
```

    ## GRanges object with 3 ranges and 1 metadata column:
    ##              seqnames    ranges strand |        name
    ##                 <Rle> <IRanges>  <Rle> | <character>
    ##   [1] human_CENPA_ORF      1-90      + |     region1
    ##   [2] human_CENPA_ORF   101-140      + |     region2
    ##   [3] human_CENPA_ORF   331-400      + |     region3
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome; no seqlengths

Now let’s convert coordinates using `convertCoordsOneSeq`:

``` r
source( here("useful_functions/multiple_sequence_alignments_functions.R") )

### by default the function returns get position in the alignment, but you could also ask for positions in other sequences in the alignment if those are in the lookup table
human_features_gr <- convertCoordsOneSeq( human_features_gr, 
                                          lookup_tbl=alnPos_lookup_table )
human_features_gr
```

    ## GRanges object with 3 ranges and 5 metadata columns:
    ##       seqnames    ranges strand |   orig_seqnames orig_start  orig_end
    ##          <Rle> <IRanges>  <Rle> |        <factor>  <integer> <integer>
    ##   [1]  aln_pos     1-156      + | human_CENPA_ORF          1        90
    ##   [2]  aln_pos   167-209      + | human_CENPA_ORF        101       140
    ##   [3]  aln_pos   400-469      + | human_CENPA_ORF        331       400
    ##       orig_width        name
    ##        <numeric> <character>
    ##   [1]         90     region1
    ##   [2]         40     region2
    ##   [3]         70     region3
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome; no seqlengths

## Extract aligned regions and concatenate

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

# Concatenate alignment slices

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

    ## R version 4.6.1 (2026-06-24)
    ## Platform: aarch64-apple-darwin23
    ## Running under: macOS Tahoe 26.6.2
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
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
    ##  [1] GenomicRanges_1.64.0 Biostrings_2.80.1    Seqinfo_1.2.0       
    ##  [4] XVector_0.52.0       IRanges_2.46.0       S4Vectors_0.50.2    
    ##  [7] BiocGenerics_0.58.1  generics_0.1.4       here_1.0.2          
    ## [10] lubridate_1.9.5      forcats_1.0.1        stringr_1.6.0       
    ## [13] dplyr_1.2.1          purrr_1.2.2          readr_2.2.0         
    ## [16] tidyr_1.3.2          tibble_3.3.1         ggplot2_4.0.3       
    ## [19] tidyverse_2.0.0     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] utf8_1.2.6         stringi_1.8.9      hms_1.1.4          digest_0.6.39     
    ##  [5] magrittr_2.0.5     evaluate_1.0.5     grid_4.6.1         timechange_0.4.0  
    ##  [9] RColorBrewer_1.1-3 fastmap_1.2.0      rprojroot_2.1.1    scales_1.4.0      
    ## [13] textshaping_1.0.5  cli_3.6.6          crayon_1.5.3       rlang_1.3.0       
    ## [17] withr_3.0.3        yaml_2.3.12        otel_0.2.0         tools_4.6.1       
    ## [21] tzdb_0.5.0         vctrs_0.7.3        R6_2.6.1           lifecycle_1.0.5   
    ## [25] ragg_1.5.2         pkgconfig_2.0.3    pillar_1.11.1      gtable_0.3.6      
    ## [29] glue_1.8.1         systemfonts_1.3.2  xfun_0.60          tidyselect_1.2.1  
    ## [33] rstudioapi_0.19.0  knitr_1.51         farver_2.1.2       htmltools_0.5.9   
    ## [37] rmarkdown_2.32     compiler_4.6.1     S7_0.2.2
