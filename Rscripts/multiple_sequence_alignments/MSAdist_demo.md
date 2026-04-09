MSAdist_demo
================
Janet Young

2026-04-06

# Play with MSA2dist package

See
[documentation](https://bioconductor.org/packages/release/bioc/vignettes/MSA2dist/inst/doc/MSA2dist.html).

Load libraries:

``` r
library(tidyverse)
library(MSA2dist)
library(Biostrings)
```

Load example data:

``` r
data(hiv, package="MSA2dist")
data(AAMatrix, package="MSA2dist")
data(woodmouse, package="ape")
```

The `hiv` dataset is a DNAStringSet object, 13 DNA seqs, each 273bp in
length (91aa)

``` r
hiv
```

    ## DNAStringSet object of length 13:
    ##      width seq                                              names               
    ##  [1]   273 GTAGTAATTAGATCTGAAAACTT...AACAATAGTTTTTAATCAATCC U68496
    ##  [2]   273 ATAGTAATTAGATCTGAAAACTT...AACAATAATCTTTAATCAATCC U68497
    ##  [3]   273 GTAGTAATTAGATCTGAAAACTT...AACAATAGTCTTTAATCAATCC U68498
    ##  [4]   273 ATAGTAATTAGATCTGAAAACTT...AACAATAGTCTTTAATCAATCC U68499
    ##  [5]   273 ATAGTAATTAGATCTGAAAACTT...AACAATAGTCTTTAATCAATCC U68500
    ##  ...   ... ...
    ##  [9]   273 ATAGTAATTAGATCTGAAAACTT...AACAATAGTTTTTAATCAATCC U68504
    ## [10]   273 ATAGTAATTAGATCTGAAAACTT...AACAATAGTCTTTAATCAATCC U68505
    ## [11]   273 ATAGTAATTAGATCTGAAAACTT...AACAATAGTTTTTAATCAATCC U68506
    ## [12]   273 GTAGTAATTCGATCTGAAAACTT...AACAATAGTCTTTAATCAATCC U68507
    ## [13]   273 ATAGTAATCAGATCTGAAAACTT...AACAATAATCTTTAATAATTCC U68508

`MSA2dist::cds2aa` translates DNA seqs, has some options regarding which
reading frame. Creates an AAStringSet object

``` r
hiv_aa <- cds2aa(hiv) 
hiv_aa
```

    ## AAStringSet object of length 13:
    ##      width seq                                              names               
    ##  [1]    91 VVIRSENFSNNAKTIIVQLNKSV...TLKQVAEKLREQFIKTIVFNQS U68496
    ##  [2]    91 IVIRSENFSNNAKTIIVQLNKSV...TVKQVAAKLREQFNKTIIFNQS U68497
    ##  [3]    91 VVIRSENFTNNAKTIIVQLNKSV...TLKQVAEKLREQFNTTIVFNQS U68498
    ##  [4]    91 IVIRSENFTNNAKTIIVQLNKSV...TLKQVAEKLREQFNKTIVFNQS U68499
    ##  [5]    91 IVIRSENFTNNAKTIIVHLNESV...TLKQVAEKLREQFNKTIVFNQS U68500
    ##  ...   ... ...
    ##  [9]    91 IVIRSENFTDNAKTIIVQLNKSV...TLRQVAQKLKEQFNRTIVFNQS U68504
    ## [10]    91 IVIRSENFTDNAKTIIVQLNKSV...TLRQVAQKLKEQFNRTIVFNQS U68505
    ## [11]    91 IVIRSENFTDNAKTIIVQLNKSV...TLRQVAKKLKEQFNRTIVFNQS U68506
    ## [12]    91 VVIRSENFTDNAKTIIVQLNKSV...ALRQVAKKLKEQFNRTIVFNQS U68507
    ## [13]    91 IVIRSENFSDNAKTIIVQLNNTV...TLRQVAEKLKEQFNKTIIFNNS U68508

``` r
simple_percent_identity_pairwise <- function(twoSeqs_stringSet) {
    if(length(twoSeqs_stringSet) != 2) {
        stop("\n\nERROR - input does has the wrong number of sequences - function designed to work on just two sequences\n\n")
    }
    if(width(twoSeqs_stringSet)[1] !=  width(twoSeqs_stringSet)[2]) {
        stop("\n\nERROR - seqs don't have the same length as each other\n\n")
    }
    output <- tibble(seq1=names(twoSeqs_stringSet)[1],
                     seq2=names(twoSeqs_stringSet)[2] ) |> 
        ## strip off description, if there is one
        mutate(seq1 = sapply(strsplit(seq1, " "), "[[", 1)) |> 
        mutate(seq2 = sapply(strsplit(seq2, " "), "[[", 1))
    
    twoSeqs <- as.matrix(twoSeqs_stringSet) |> 
        t() |> 
        as_tibble()
    colnames(twoSeqs) <- c("seq1", "seq2")
    
    twoSeqs_degapFullGap <- twoSeqs |> 
        filter(!(seq1=="-" & seq2=="-"))
    output$aln_len <- nrow(twoSeqs_degapFullGap)
    
    output$num_identical <- sum(twoSeqs_degapFullGap$seq1 == twoSeqs_degapFullGap$seq2)
    
    twoSeqs_degapAnyGap <- twoSeqs_degapFullGap |> 
        filter(seq1 != "-" & seq2 != "-")
    output$aln_len_nogaps <- nrow(twoSeqs_degapAnyGap)
    
    output <- output |> 
        mutate(pid_incl_gap = 100*num_identical / aln_len) |> 
        mutate(pid_excl_gap = 100*num_identical / aln_len_nogaps)
    return(output)
}

simple_percent_identity_multiple <- function(one_aln) {
    if(length(unique(width(one_aln)))>1) {
        stop("\n\nERROR - these are not aligned seqs - they have different lengths: ", 
             paste(width(one_aln), collapse=","),
             "\n\n")
    }
    lapply(1:(length(one_aln)-1), function(i) {
        lapply( (i+1):length(one_aln), function(j) {
            simple_percent_identity_pairwise(one_aln[c(i,j)])
        }) |> 
            bind_rows()
    }) |> 
        bind_rows()
}
```

``` r
hiv_aa_myDists <- simple_percent_identity_multiple(hiv_aa[1:4])
as.matrix(hiv_aa_myDists)
```

    ##      seq1     seq2     aln_len num_identical aln_len_nogaps pid_incl_gap
    ## [1,] "U68496" "U68497" "91"    "85"          "91"           "93.40659"  
    ## [2,] "U68496" "U68498" "91"    "75"          "91"           "82.41758"  
    ## [3,] "U68496" "U68499" "91"    "75"          "91"           "82.41758"  
    ## [4,] "U68497" "U68498" "91"    "71"          "91"           "78.02198"  
    ## [5,] "U68497" "U68499" "91"    "73"          "91"           "80.21978"  
    ## [6,] "U68498" "U68499" "91"    "79"          "91"           "86.81319"  
    ##      pid_excl_gap
    ## [1,] "93.40659"  
    ## [2,] "82.41758"  
    ## [3,] "82.41758"  
    ## [4,] "78.02198"  
    ## [5,] "80.21978"  
    ## [6,] "86.81319"

``` r
my_percent_id_to_matrix <- function(my_percent_id_tbl,
                                    col_to_take="pid_excl_gap") {
    my_seqorder <- unique(c(my_percent_id_tbl$seq1,
                            my_percent_id_tbl$seq2))
    ## first make it symmetrical
    my_percent_id_tbl_flip <- my_percent_id_tbl |> 
        dplyr::rename(seq1_new=seq2, 
                      seq2_new=seq1) |> 
        dplyr::rename(seq1=seq1_new, 
                      seq2=seq2_new)
    my_percent_id_tbl_both <- bind_rows(my_percent_id_tbl,
                                        my_percent_id_tbl_flip) 
    # return(my_percent_id_tbl_both)
    my_percent_id_tbl_both <- my_percent_id_tbl_both |> 
        select(seq1, seq2, value=all_of(col_to_take))
    # return(my_percent_id_tbl_both)
    my_matrix <- my_percent_id_tbl_both |>
        pivot_wider(id_cols="seq1",
                    values_from="value",
                    names_from="seq2", 
                    values_fill=0)
    
    my_matrix_v2 <- as.matrix(my_matrix[,my_seqorder])
    rownames(my_matrix_v2) <- my_matrix$seq1
    my_matrix_v2 <- my_matrix_v2[my_seqorder,]
    return(my_matrix_v2)
}
# temp <- hiv_aa_myDists |> my_percent_id_to_matrix() 
temp <- hiv_aa_myDists |> 
    mutate(num_diff=aln_len-num_identical) |> 
    my_percent_id_to_matrix(col_to_take="num_diff") 
temp 
```

    ##        U68496 U68497 U68498 U68499
    ## U68496      0      6     16     16
    ## U68497      6      0     20     18
    ## U68498     16     20      0     12
    ## U68499     16     18     12      0

# Use MSA2dist distance functions

aastring2dist gets a distances matrix. We need to specify scoring
matrix - the vignette uses Grantham scores.

Returns a list of three items \[1\] “distSTRING” “sitesUsed”
“regionUsed”

Something is wrong on my computer! doesn’t display right at all. Should
be a numeric matrix.

``` r
hiv_aa_dist <- aastring2dist(hiv_aa, 
                             score=granthamMatrix())
```

    ## Computing: [========================================] 100% (done)

``` r
# lapply(hiv_aa_dist, class)
# hiv_aa_dist[["distSTRING"]] |> head()
# head(hiv_aa_dist$distSTRING)

# distSTRING = matrix of distances
# sitesUsed = matrix of number of sites used
# regionUsed = IRanges object specifying which part of the alignment was used (default is the whole thing)
```

Show the distances between the first 4 seqs (wrap in as.matrix for now,
because of a weird rendering bug probably with Rstudio)

``` r
as.matrix(hiv_aa_dist[["distSTRING"]][1:4, 1:4])
```

    ##           U68496   U68497    U68498    U68499
    ## U68496  0.000000  4.43956 11.516484  9.879121
    ## U68497  4.439560  0.00000 12.681319 10.406593
    ## U68498 11.516484 12.68132  0.000000  8.769231
    ## U68499  9.879121 10.40659  8.769231  0.000000

Printing this matrix looks really wierd in the Rstudio/Rmd rendering but
it looks OK under most other circumstances, e.g. if I wrap in
`as.matrix()`, when I knit, and if I print in the console. I think
something is going wrong in my Rstudio. Probably need a new version or
to restart my computer or something.

To use a different score matrix, here as an example the AAMatrix from
the R package alakazam is used:

``` r
head(AAMatrix)
```

    ##   A B C D E F G H I J K L M N P Q R S T V W X Y Z * - .
    ## A 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 0 0
    ## B 1 0 1 1 1 1 1 1 1 1 1 1 1 0 1 1 0 1 1 1 1 0 1 1 1 0 0
    ## C 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 0 0
    ## D 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 0 0
    ## E 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 1 0 0
    ## F 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 0 0

These numbers are the parwise scores calculated from the scoring matrix,
divided by the aligned length:

``` r
aa.dist.AAMatrix <- aastring2dist(hiv_aa, score=AAMatrix)
```

    ## Computing: [========================================] 100% (done)

``` r
as.matrix(aa.dist.AAMatrix[["distSTRING"]][1:4, 1:4])
```

    ##            U68496     U68497    U68498    U68499
    ## U68496 0.00000000 0.06593407 0.1758242 0.1758242
    ## U68497 0.06593407 0.00000000 0.2197802 0.1978022
    ## U68498 0.17582418 0.21978022 0.0000000 0.1318681
    ## U68499 0.17582418 0.19780220 0.1318681 0.0000000

Multiply by the aligned length to get the number of mismatches (this is
the same as my function gives):

``` r
as.matrix(aa.dist.AAMatrix[["distSTRING"]][1:4, 1:4] * 91)
```

    ##        U68496 U68497 U68498 U68499
    ## U68496      0      6     16     16
    ## U68497      6      0     20     18
    ## U68498     16     20      0     12
    ## U68499     16     18     12      0

# DNA distance

``` r
dna.dist <- dnastring2dist(hiv, model="K80")
```

    ## Computing: [========================================] 100% (done)

``` r
as.matrix(dna.dist[["distSTRING"]][1:4, 1:4])
```

    ##            U68496     U68497     U68498     U68499
    ## U68496 0.00000000 0.03381189 0.07731910 0.08135801
    ## U68497 0.03381189 0.00000000 0.09396372 0.08960527
    ## U68498 0.07731910 0.09396372 0.00000000 0.05333488
    ## U68499 0.08135801 0.08960527 0.05333488 0.00000000

# Ka/Ks

There are tons of models available for Ka/Ks calculations

``` r
ka_ks_li <- dnastring2kaks(hiv, model="Li")
```

    ## Joining with `by = join_by(seq1, seq2)`
    ## Joining with `by = join_by(seq1, seq2)`
    ## Joining with `by = join_by(seq1, seq2)`

``` r
head( as.matrix( ka_ks_li ) )
```

    ##      Comp1 Comp2 seq1     seq2     ka           ks              vka           
    ## [1,] " 1"  " 2"  "U68496" "U68497" "0.03026357" " 0.0317031856" "0.0003051202"
    ## [2,] " 1"  " 3"  "U68496" "U68498" "0.09777332" " 0.0176141589" "0.0009314173"
    ## [3,] " 1"  " 4"  "U68496" "U68499" "0.10295875" " 0.0176731064" "0.0009595527"
    ## [4,] " 1"  " 5"  "U68496" "U68500" "0.13461355" " 0.0463968986" "0.0013731885"
    ## [5,] " 1"  " 6"  "U68496" "U68501" "0.12607831" " 0.0284429395" "0.0013277282"
    ## [6,] " 1"  " 7"  "U68496" "U68502" "0.17441037" " 0.1092653245" "0.0017871934"
    ##      vks            Ka           Ks              Ka/Ks         
    ## [1,] "0.0004007730" "0.03026357" " 0.0317031856" "   0.9545907"
    ## [2,] "0.0005091970" "0.09777332" " 0.0176141589" "   5.5508365"
    ## [3,] "0.0005787387" "0.10295875" " 0.0176731064" "   5.8257302"
    ## [4,] "0.0020497481" "0.13461355" " 0.0463968986" "   2.9013480"
    ## [5,] "0.0006310195" "0.12607831" " 0.0284429395" "   4.4326750"
    ## [6,] "0.0040766687" "0.17441037" " 0.1092653245" "   1.5962097"

``` r
ka_ks_YN <- dnastring2kaks(hiv, model="YN")
head( as.matrix( ka_ks_YN ) )
```

    ##               Comp1 Comp2 seq1     seq2     Method Ka          Ks         
    ## U68496_U68497 "1"   "2"   "U68496" "U68497" "YN"   "0.0259337" "0.0848079"
    ## U68496_U68498 "1"   "3"   "U68496" "U68498" "YN"   "0.0833439" "0.0402529"
    ## U68496_U68499 "1"   "4"   "U68496" "U68499" "YN"   "0.0888328" "0.0375284"
    ## U68496_U68500 "1"   "5"   "U68496" "U68500" "YN"   "0.11672"   "0.074978" 
    ## U68496_U68501 "1"   "6"   "U68496" "U68501" "YN"   "0.108488"  "0.0673923"
    ## U68496_U68502 "1"   "7"   "U68496" "U68502" "YN"   "0.153774"  "0.135014" 
    ##               Ka/Ks      P-Value(Fisher) Length S-Sites   N-Sites  
    ## U68496_U68497 "0.305794" "0.110988"      "273"  "37.4117" "235.588"
    ## U68496_U68498 "2.07051"  "0.441303"      "273"  "36.6562" "236.344"
    ## U68496_U68499 "2.36709"  "0.471701"      "273"  "38.0351" "234.965"
    ## U68496_U68500 "1.55673"  "0.547288"      "273"  "37.4853" "235.515"
    ## U68496_U68501 "1.60979"  "0.540854"      "273"  "39.4297" "233.57" 
    ## U68496_U68502 "1.13894"  "0.811025"      "273"  "44.966"  "228.034"
    ##               Fold-Sites(0:2:4) Substitutions S-Substitutions N-Substitutions
    ## U68496_U68497 "NA"              "9"           "3"             "6"            
    ## U68496_U68498 "NA"              "20"          "1.42507"       "18.5749"      
    ## U68496_U68499 "NA"              "21"          "1.38331"       "19.6167"      
    ## U68496_U68500 "NA"              "28"          "2.62481"       "25.3752"      
    ## U68496_U68501 "NA"              "26"          "2.49647"       "23.5035"      
    ## U68496_U68502 "NA"              "37"          "5.51119"       "31.4888"      
    ##               Fold-S-Substitutions(0:2:4) Fold-N-Substitutions(0:2:4)
    ## U68496_U68497 "NA"                        "NA"                       
    ## U68496_U68498 "NA"                        "NA"                       
    ## U68496_U68499 "NA"                        "NA"                       
    ## U68496_U68500 "NA"                        "NA"                       
    ## U68496_U68501 "NA"                        "NA"                       
    ## U68496_U68502 "NA"                        "NA"                       
    ##               Divergence-Time
    ## U68496_U68497 "0.0340018"    
    ## U68496_U68498 "0.077558"     
    ## U68496_U68499 "0.081685"     
    ## U68496_U68500 "0.110989"     
    ## U68496_U68501 "0.102552"     
    ## U68496_U68502 "0.150684"     
    ##               Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)
    ## U68496_U68497 "2.90975:2.90975:1:1:1:1"                           
    ## U68496_U68498 "2.92242:2.92242:1:1:1:1"                           
    ## U68496_U68499 "3.24306:3.24306:1:1:1:1"                           
    ## U68496_U68500 "2.4231:2.4231:1:1:1:1"                             
    ## U68496_U68501 "2.49234:2.49234:1:1:1:1"                           
    ## U68496_U68502 "3.27339:3.27339:1:1:1:1"                           
    ##               GC(1:2:3)                              ML-Score AICc
    ## U68496_U68497 "0.300366(0.362637:0.368132:0.17033)"  "NA"     "NA"
    ## U68496_U68498 "0.305861(0.373626:0.373626:0.17033)"  "NA"     "NA"
    ## U68496_U68499 "0.304029(0.373626:0.362637:0.175824)" "NA"     "NA"
    ## U68496_U68500 "0.300366(0.362637:0.362637:0.175824)" "NA"     "NA"
    ## U68496_U68501 "0.302198(0.362637:0.362637:0.181319)" "NA"     "NA"
    ## U68496_U68502 "0.320513(0.368132:0.384615:0.208791)" "NA"     "NA"
    ##               Akaike-Weight Model
    ## U68496_U68497 "NA"          "NA" 
    ## U68496_U68498 "NA"          "NA" 
    ## U68496_U68499 "NA"          "NA" 
    ## U68496_U68500 "NA"          "NA" 
    ## U68496_U68501 "NA"          "NA" 
    ## U68496_U68502 "NA"          "NA"

## average behavior each codon

dnastring2codonmat turns a DNA alignment into a matrix of codons. One
row per codon slice of the alignment, one column per sequence in the
alignment

codonmat2xy does calculations for each codon position. “This function
calculates average behavior of each codon for all pairwise comparisons
for indels, syn, and nonsyn mutations according to Nei and Gojobori
(1986).”

Output:

- n is number of pairwise comparisons

``` r
hiv.xy <-  codonmat2xy(dnastring2codonmat(hiv))
```

    ## Joining with `by = join_by(Codon)`
    ## Joining with `by = join_by(Codon)`
    ## Joining with `by = join_by(Codon)`

``` r
head(as.matrix(hiv.xy))
```

    ##      Codon  n SynSum NonSynSum IndelSum   SynMean NonSynMean IndelMean
    ## [1,]     1 78      0        40        0 0.0000000  0.5128205         0
    ## [2,]     2 78      0         0        0 0.0000000  0.0000000         0
    ## [3,]     3 78     12         0        0 0.1538462  0.0000000         0
    ## [4,]     4 78     12         0        0 0.1538462  0.0000000         0
    ## [5,]     5 78      0         0        0 0.0000000  0.0000000         0
    ## [6,]     6 78      0         0        0 0.0000000  0.0000000         0
    ##      CumSumSynMean CumSumNonSynMean CumSumIndelMean
    ## [1,]     0.0000000        0.5128205               0
    ## [2,]     0.0000000        0.5128205               0
    ## [3,]     0.1538462        0.5128205               0
    ## [4,]     0.3076923        0.5128205               0
    ## [5,]     0.3076923        0.5128205               0
    ## [6,]     0.3076923        0.5128205               0

## Grantham scores

I found [this
info](https://gist.github.com/danielecook/501f03650bca6a3db31ff3af2d413d2a),
not sure how correct it is:

[Grantham
scores](https://www.science.org/doi/10.1126/science.185.4154.862), which
categorize codon replacements into classes of increasing chemical
dissimilarity, were designated conservative (0-50), moderately
conservative (51-100), moderately radical (101-150), or radical (≥151)
according to the classification proposed by Li et al. (28). These
numbers are estimated by looking at the frequency of all pair-wise amino
acid substitutions between species.

# Finished

``` r
sessionInfo()
```

    ## R version 4.5.3 (2026-03-11)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Tahoe 26.4
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
    ##  [1] Biostrings_2.78.0   Seqinfo_1.0.0       XVector_0.50.0     
    ##  [4] IRanges_2.44.0      S4Vectors_0.48.0    BiocGenerics_0.56.0
    ##  [7] generics_0.1.4      MSA2dist_1.14.0     lubridate_1.9.5    
    ## [10] forcats_1.0.1       stringr_1.6.0       dplyr_1.2.1        
    ## [13] purrr_1.2.1         readr_2.2.0         tidyr_1.3.2        
    ## [16] tibble_3.3.1        ggplot2_4.0.2       tidyverse_2.0.0    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] lattice_0.22-9       stringi_1.8.7        hms_1.1.4           
    ##  [4] digest_0.6.39        magrittr_2.0.5       evaluate_1.0.5      
    ##  [7] grid_4.5.3           timechange_0.4.0     RColorBrewer_1.1-3  
    ## [10] iterators_1.0.14     seqinr_4.2-36        fastmap_1.2.0       
    ## [13] foreach_1.5.2        doParallel_1.0.17    ape_5.8-1           
    ## [16] scales_1.4.0         ade4_1.7-24          codetools_0.2-20    
    ## [19] cli_3.6.5            rlang_1.1.7          crayon_1.5.3        
    ## [22] withr_3.0.2          yaml_2.3.12          otel_0.2.0          
    ## [25] parallel_4.5.3       tools_4.5.3          tzdb_0.5.0          
    ## [28] vctrs_0.7.2          R6_2.6.1             lifecycle_1.0.5     
    ## [31] pwalign_1.6.0        MASS_7.3-65          pkgconfig_2.0.3     
    ## [34] pillar_1.11.1        gtable_0.3.6         glue_1.8.0          
    ## [37] Rcpp_1.1.1           GenomicRanges_1.62.1 xfun_0.57           
    ## [40] tidyselect_1.2.1     rstudioapi_0.18.0    knitr_1.51          
    ## [43] farver_2.1.2         nlme_3.1-169         htmltools_0.5.9     
    ## [46] rmarkdown_2.31       compiler_4.5.3       S7_0.2.1
