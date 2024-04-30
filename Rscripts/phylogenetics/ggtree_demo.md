ggtree_demo
================
Janet Young

2024-04-29

# Useful links, general advice

[ggtree manual](http://yulab-smu.top/treedata-book/index.html)

Chapter 1 looks at reading/writing trees. The treeio and ape packages
both have tree read/write functions).

## Read a tree and plot it

read in a newick tree file:

``` r
nwk_file <- system.file("extdata/sample.nwk", package="treeio") 
nwk_tree <- read.tree( nwk_file )
class(nwk_tree)
```

    ## [1] "phylo"

``` r
nwk_tree
```

    ## 
    ## Phylogenetic tree with 13 tips and 12 internal nodes.
    ## 
    ## Tip labels:
    ##   A, B, C, D, E, F, ...
    ## 
    ## Rooted; includes branch lengths.

plot it using ape’s basic plot function

``` r
plot(nwk_tree)
```

![](ggtree_demo_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

plot it using ggtree

``` r
nwk_tree %>% 
    ggtree() +
    geom_tiplab()
```

![](ggtree_demo_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Associate a tree with addtional data

Make a fake data tibble with info on each of the taxa in nkw_tree

``` r
num_taxa <- length(nwk_tree$tip.label)
tip_dat <- tibble( taxon=nwk_tree$tip.label,
                   fake_height= rnorm(n=num_taxa, mean=100, sd=30 ),
                   fake_phenotype= sample( c("teeth","nose","tail"), size=num_taxa, replace=TRUE ) )
tip_dat
```

    ## # A tibble: 13 × 3
    ##    taxon fake_height fake_phenotype
    ##    <chr>       <dbl> <chr>         
    ##  1 A           108.  tail          
    ##  2 B            64.7 nose          
    ##  3 C           116.  teeth         
    ##  4 D           115.  nose          
    ##  5 E           111.  tail          
    ##  6 F            98.6 teeth         
    ##  7 G            45.6 teeth         
    ##  8 H            80.0 teeth         
    ##  9 I            63.1 tail          
    ## 10 J           103.  nose          
    ## 11 K           158.  nose          
    ## 12 L            74.5 tail          
    ## 13 M            83.6 tail

## combine tree and tbl

Associate information tbl with the tree - this creates a `treedata`
object that contains the tree (access it using
`nwk_tree_with_info@phylo`)

Advice:  
- this joining can be flaky. - BEFORE you combine tree with info tbl: -
manipulate tree to get as close to the final product as you can
(e.g. reroot, subset, rotate clades) - manipulate tbl get as close to
the final product as you can (e.g. join a bunch of info togeter) - CHECK
taxon labels in tree and the tbl column you plan to join on match up
OK - CHECK the taxon ID in the join column of tbl is UNIQUE - repeated
values will cause trouble

``` r
nwk_tree_with_info <- left_join(nwk_tree, tip_dat, by=c("label"="taxon"))
class(nwk_tree_with_info)
```

    ## [1] "treedata"
    ## attr(,"package")
    ## [1] "tidytree"

``` r
nwk_tree_with_info
```

    ## 'treedata' S4 object'.
    ## 
    ## ...@ phylo:
    ## 
    ## Phylogenetic tree with 13 tips and 12 internal nodes.
    ## 
    ## Tip labels:
    ##   A, B, C, D, E, F, ...
    ## 
    ## Rooted; includes branch lengths.
    ## 
    ## with the following features available:
    ##   '', 'fake_height', 'fake_phenotype'.
    ## 
    ## # The associated data tibble abstraction: 25 × 5
    ## # The 'node', 'label' and 'isTip' are from the phylo tree.
    ##     node label isTip fake_height fake_phenotype
    ##    <int> <chr> <lgl>       <dbl> <chr>         
    ##  1     1 A     TRUE        108.  tail          
    ##  2     2 B     TRUE         64.7 nose          
    ##  3     3 C     TRUE        116.  teeth         
    ##  4     4 D     TRUE        115.  nose          
    ##  5     5 E     TRUE        111.  tail          
    ##  6     6 F     TRUE         98.6 teeth         
    ##  7     7 G     TRUE         45.6 teeth         
    ##  8     8 H     TRUE         80.0 teeth         
    ##  9     9 I     TRUE         63.1 tail          
    ## 10    10 J     TRUE        103.  nose          
    ## # ℹ 15 more rows

``` r
nwk_tree_with_info %>% 
    ggtree(aes(color=fake_phenotype)) +  # color here colors the branches
    geom_tiplab(aes(color=fake_phenotype)) + # color here colors the tip labels
    geom_tippoint(aes(subset= fake_phenotype=="tail"), color="black") # add dots to taxa, using subset
```

![](ggtree_demo_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
