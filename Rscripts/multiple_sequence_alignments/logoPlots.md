sequence logo plots
================
Janet Young

2025-08-29

Goal - show how to make logo plots

# Load libraries

``` r
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(Biostrings)
library(ggseqlogo)
library(DiffLogo)
library(ggmsa)
```

# Logo plots

There are many packages that can make logo plots. Some are:  
\* [ggseqlogo](https://omarwagih.github.io/ggseqlogo/) - see demo
below.  
\*
[DiffLogo](https://bioconductor.org/packages/release/bioc/html/DiffLogo.html) -
see demo below  
\* [ggmsa](http://yulab-smu.top/ggmsa/) - see demo below  
\* [RWebLogo](https://github.com/WebLogo/weblogo) - see
[manual](http://weblogo.threeplusone.com/manual.html). Was not available
for R version 3.4.0 - I haven’t checked more recent versions. Seems more
a python thing. I think it does regular logos but not difference logos  
\*
[motifStack](https://bioconductor.org/packages/release/bioc/vignettes/motifStack/inst/doc/motifStack_HTML.html) -
I don’t think it does difference logos (although there is something
called affinity logos, for TF affinities)  
\*
[seqLogo](https://bioconductor.org/packages/release/bioc/vignettes/seqLogo/inst/doc/seqLogo.html) -
I think this only works with DNAseqs

## Read in example data: short H2A alignment (Antoine)

In this example I have a single alignment of the histone fold domain of
various H2A family members.

I read in the alignment

``` r
## same code is found in jensenShannonDistance.Rmd
aln_file <- "example_alignment_files/exampleProtAln_shortH2As_histoneFoldDomain.fa"

shortH2Aaln <- readAAStringSet(aln_file)
# simplify the sequence names by removing the description
names(shortH2Aaln) <- sapply(strsplit(names(shortH2Aaln), " "), "[[", 1)
```

## ggseqlogo package - quick demo

See [ggseqlogo documentation](https://omarwagih.github.io/ggseqlogo/).

Nice.

Looks like you can make custom color scheme (see `?make_col_scheme`) but
I haven’t tried it.

You can make ‘custom height’ logos that allow negative values. This may
provide a way to make difference logos, but I think we’d have to do the
calculations ourselves about letter heights.

I think you can also supply frequency or count matrices (but I haven’t
tried it)

``` r
# ?ggseqlogo
# ?geom_logo
# ?list_col_schemes
# ?make_col_scheme
ggseqlogo(as.character(shortH2Aaln),
          col_scheme="chemistry2") +
    # turn the x axis labels 90 degrees and change font size
    theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=0.5, size=7)) +
    # suppress the color scheme legend:
    guides(fill = "none")  
```

![](logoPlots_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

## ggmsa quick demo

Makes a basic logo plot above the corresponding alignment.

See [ggmsa documentation](http://yulab-smu.top/ggmsa/).

Not sure if there’s a way to get rid of the alignment. I don’t think
there’s a way to make the height of each letter stack reflect the
entropy. This is probably too basic for our use.

``` r
ggmsa(shortH2Aaln,
      # font = NULL,
      color = "Chemistry_AA") + 
    geom_seqlogo(adaptive=FALSE)  # adaptive=FALSE makes the logo plot taller, but whether T or F the overall heights of each stack are the same
```

![](logoPlots_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Comparing several alignments

To get example data, I’ll first split the short H2A alignment into
groups.

It contains canonical H2A, H2A.B, H2A.L, H2A.P (8 species each) and the
marsupial-specific H2A.R (5 sequences) (total of 37 sequences).

Here I figure out which H2A variant each sequence is from and split the
single master alignment into 5 individual alignments accordingly.

``` r
## figure out which variant each sequence is from
shortH2AalnSeqTypes <- sapply(strsplit(names(shortH2Aaln), "_"), "[[", 1)
shortH2AalnSeqTypes <- gsub("R[12]$","R", shortH2AalnSeqTypes)

## split the alignment into those categories. We get a list object containing all 5 alignments
shortH2AalnSplit <- split(shortH2Aaln, shortH2AalnSeqTypes)

## I add a sixth alignment, which is the combined B, L and P alignment
shortH2AalnSplit[["BandLandP"]] <- c(shortH2AalnSplit[["H2A.B"]], shortH2AalnSplit[["H2A.L"]], shortH2AalnSplit[["H2A.P"]])
```

To show several ggseqlogo plots above each other:

``` r
shortH2AalnSplit_chars <- lapply(shortH2AalnSplit, as.character)

ggseqlogo(shortH2AalnSplit_chars, ncol=1) +
    theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=0.5, size=7)) +
    guides(fill = "none")  
```

![](logoPlots_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## DiffLogo demo, default color scheme

Difflogo is helpful in looking at frequency differences between groups,
but it wants the data in frequency matrix format, rather than as
sequence alignments.

So first we define a function called `getAlnCounts()` that takes an
alignment and returns a count or frequency matrix

``` r
# checkAlnLengths is a tiny function that makes sure all seqs in an alignment are the same length as each other
checkAlnLengths <- function(aln) {
    if(length(unique(width(aln))) != 1) {
        stop("\n\nERROR - you supplied a ragged alignment (seqs not all the same length)\n\n")
    } else {
        return(TRUE)
    } 
}

## define the letters we want to count
# AA_STANDARD is defined in the Biostrings package and includes the usual 20 amino acids. I want to add the gap character ("-")
myAAtoTabulate <- c("-", AA_STANDARD)

## define the function
getAlnCounts <- function(aln, letters=myAAtoTabulate, as.prob=FALSE) {
    # check for ragged alns (seqs not all the same length)
    checkAlnLengths(aln)
    
    # get counts
    countsEachSeq <- lapply(1:length(aln), function(i) {
        letterFrequencyInSlidingView(aln[[i]], view.width = 1, letters=letters)
    })
    
    # if there were letters in the alignment that are not accounted for in the letters argument, the totals won't be correct.
    expectedTotals <- width(aln)[1]
    totalCountsEachSeq <- sapply(countsEachSeq, sum)
    if ( sum(totalCountsEachSeq != expectedTotals) > 0) {
        stop("\n\nERROR - the total counts didn't add up correctly. Are there letters in the alignment that are not present in the letters argument you supplied?\n\n")
    }
    
    # get total counts by position - the Reduce function takes a list object and uses the specified function on all the elements
    countTotals <- Reduce("+", countsEachSeq)
    
    # transpose so columns are positions and rows are each letter type
    countTotals <- t(countTotals)
    
    # perhaps get frequencies not counts
    if(as.prob) {
        freqs <- countTotals / colSums(countTotals)
        return(freqs)
    } else {
        return(countTotals)
    }
}
```

Now use the `getAlnCounts()` function on each of our alignments to get
frequency matrices:

``` r
shortH2AalnSplit_freqs <- lapply( shortH2AalnSplit, 
                                  getAlnCounts, 
                                  letters=myAAtoTabulate, 
                                  as.prob=TRUE)
```

See [DiffLogo
documentation](https://bioconductor.org/packages/release/bioc/html/DiffLogo.html).

I used DiffLogo for Rossana’s data, because I was very motivated to make
**difference** logos, not just single-alignment logo plots. In Mar 21,
2022 I think some of the code given in the vignette may have been broken

(it always gives an annoying bonus message, “pwm must be of class matrix
or data.frame. Trying to convert”, but if there are no other errors we
can ignore it)

First try a simple seqLogo plot for the H2A.B alignment

    ## [1] "pwm must be of class matrix or data.frame. Trying to convert"

![](logoPlots_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

## DiffLogo plot, controlling the color scheme

### first, set up a custom color scheme

this color scheme mimics the “chemistry” color scheme from motifStack
(defined
[here](https://github.com/jianhong/motif%20Stack/blob/5aa80388b44bc8f93738315310cfb56c3495130c/R/publicUtilities.R))

``` r
# changeColors is a function I wrote that works on objects of class "Alphabet" to change the colors for some amino acids (or nucleotides)
# myAlphabet is the Alphabet object (the color scheme we're working on)
# myAA is a vector of the AA (or nuc) we want to change
# myCol is the new color (e.g. "black")
changeColors <- function(myAlphabet, myAA, myCol) {
    # get indices of things we want to change
    whichChar <- which(myAlphabet$chars %in% myAA)
    # get existing color scheme
    tempColors <- myAlphabet$cols
    # change the ones we selected
    tempColors[whichChar] <- myCol
    # add the altered color scheme back to the Alphabet object
    myAlphabet$cols <- tempColors
    # done
    return(myAlphabet)
}

### ASN is a DiffLogo object (see ?ASN) that defines their default color scheme
### make custom version of the ASN Alphabet object called FULL_ALPHABET_JYchemistryColors
ASN_JYchemistryColors <- ASN

ASN_JYchemistryColors <- changeColors(ASN_JYchemistryColors,
                                      c("A","F","I","L","M","P","V","W"), 
                                      "black")

ASN_JYchemistryColors <- changeColors(ASN_JYchemistryColors, 
                                      c("C","G","S","T","Y"), 
                                      "forest green")

ASN_JYchemistryColors <- changeColors(ASN_JYchemistryColors, 
                                      c("D","E"), 
                                      "red3")

ASN_JYchemistryColors <- changeColors(ASN_JYchemistryColors, 
                                      c("H","K","R"), 
                                      "blue3")

ASN_JYchemistryColors <- changeColors(ASN_JYchemistryColors, 
                                      c("N","Q"), 
                                      "magenta4")
```

``` r
DiffLogo::seqLogo(pwm=shortH2AalnSplit_freqs_justASN[["H2A.B"]], 
                  alphabet=ASN_JYchemistryColors, 
                  drawLines=20) 
```

    ## [1] "pwm must be of class matrix or data.frame. Trying to convert"

![](logoPlots_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

# Now show difference logos

``` r
diffLogoFromPwm(
    pwm1 = shortH2AalnSplit_freqs_justASN[["H2A"]],
    pwm2 = shortH2AalnSplit_freqs_justASN[["H2A.B"]],
    # ymin=0.1, ymax=-0.1,
    alphabet = ASN_JYchemistryColors)
```

    ## [1] "pwm must be of class matrix or data.frame. Trying to convert"
    ## [1] "pwm must be of class matrix or data.frame. Trying to convert"

![](logoPlots_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

# Finished

show R version used, and package versions

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
    ## [1] grid      stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] ggmsa_1.14.1        DiffLogo_2.32.0     cba_0.2-25         
    ##  [4] proxy_0.4-27        ggseqlogo_0.2       Biostrings_2.76.0  
    ##  [7] GenomeInfoDb_1.44.1 XVector_0.48.0      IRanges_2.42.0     
    ## [10] S4Vectors_0.46.0    BiocGenerics_0.54.0 generics_0.1.4     
    ## [13] lubridate_1.9.4     forcats_1.0.0       stringr_1.5.1      
    ## [16] dplyr_1.1.4         purrr_1.1.0         readr_2.1.5        
    ## [19] tidyr_1.3.1         tibble_3.3.0        ggplot2_3.5.2      
    ## [22] tidyverse_2.0.0    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6            xfun_0.52               lattice_0.22-7         
    ##  [4] tzdb_0.5.0              vctrs_0.6.5             tools_4.5.1            
    ##  [7] yulab.utils_0.2.0       parallel_4.5.1          pkgconfig_2.0.3        
    ## [10] ggplotify_0.1.2         RColorBrewer_1.1-3      lifecycle_1.0.4        
    ## [13] GenomeInfoDbData_1.2.14 compiler_4.5.1          farver_2.1.2           
    ## [16] treeio_1.32.0           ggforce_0.5.0           ggtree_3.16.3          
    ## [19] ggfun_0.2.0             htmltools_0.5.8.1       lazyeval_0.2.2         
    ## [22] yaml_2.3.10             pillar_1.11.0           crayon_1.5.3           
    ## [25] MASS_7.3-65             seqmagick_0.1.7         nlme_3.1-168           
    ## [28] tidyselect_1.2.1        aplot_0.2.8             digest_0.6.37          
    ## [31] stringi_1.8.7           labeling_0.4.3          polyclip_1.10-7        
    ## [34] fastmap_1.2.0           cli_3.6.5               magrittr_2.0.3         
    ## [37] patchwork_1.3.1         ape_5.8-1               withr_3.0.2            
    ## [40] scales_1.4.0            UCSC.utils_1.4.0        timechange_0.3.0       
    ## [43] rmarkdown_2.29          httr_1.4.7              hms_1.1.3              
    ## [46] evaluate_1.0.4          knitr_1.50              gridGraphics_0.5-1     
    ## [49] rlang_1.1.6             Rcpp_1.1.0              tidytree_0.4.6         
    ## [52] glue_1.8.0              R4RNA_1.36.0            tweenr_2.0.3           
    ## [55] rstudioapi_0.17.1       jsonlite_2.0.0          R6_2.6.1               
    ## [58] fs_1.6.6
