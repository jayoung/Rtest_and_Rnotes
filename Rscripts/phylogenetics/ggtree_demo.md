ggtree_demo
================
Janet Young

2024-05-10

# Useful links, general advice

The [ggtree manual](http://yulab-smu.top/treedata-book/index.html) is
great (very detailed). - Chapter 1 looks at reading/writing trees. The
treeio and ape packages both have tree read/write functions).

# Read a tree and plot it

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

# Associate a tree with additional data

Make a fake data tibble with info on each of the taxa in nwk_tree

``` r
num_taxa <- length(nwk_tree$tip.label)
tip_dat <- tibble( taxon=nwk_tree$tip.label,
                   fake_height= rnorm(n=num_taxa, mean=100, sd=30 ),
                   fake_phenotype= sample( c("teeth","nose","tail"), 
                                           size=num_taxa, replace=TRUE ),
                   fake_phenotype2=sample( c("fur","scales","spikes"), 
                                           size=num_taxa, replace=TRUE ) )
tip_dat
```

    ## # A tibble: 13 × 4
    ##    taxon fake_height fake_phenotype fake_phenotype2
    ##    <chr>       <dbl> <chr>          <chr>          
    ##  1 A            77.3 tail           scales         
    ##  2 B            49.1 tail           spikes         
    ##  3 C            81.8 tail           fur            
    ##  4 D           143.  teeth          fur            
    ##  5 E           106.  tail           fur            
    ##  6 F            95.2 tail           scales         
    ##  7 G           104.  tail           spikes         
    ##  8 H            90.2 teeth          scales         
    ##  9 I            86.8 nose           spikes         
    ## 10 J           115.  teeth          spikes         
    ## 11 K           128.  nose           fur            
    ## 12 L            56.6 teeth          spikes         
    ## 13 M            77.8 teeth          fur

## combine tree and tbl

We associate info tbl with the tree using `left_join` - this creates a
single `treedata` object that contains info AND tree (we access the tree
like this `nwk_tree_with_info@phylo`)

Advice:  
- this joining process can be flaky - there’s not as much error-checking
built-in as I’d like.  
- BEFORE you combine tree with info tbl:  
- you should manipulate tree to get as close to the final product as you
can (e.g. do any rerooting, subsetting, clade rotation, etc, BEFORE you
left_join)  
- manipulate tbl get as close to the final product as you can (e.g. join
a bunch of info together)  
- CHECK taxon labels in tree and the tbl column you plan to join on
match up OK  
- CHECK the taxon ID in the join column of tbl is UNIQUE - repeated
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
    ##   '', 'fake_height', 'fake_phenotype', 'fake_phenotype2'.
    ## 
    ## # The associated data tibble abstraction: 25 × 6
    ## # The 'node', 'label' and 'isTip' are from the phylo tree.
    ##     node label isTip fake_height fake_phenotype fake_phenotype2
    ##    <int> <chr> <lgl>       <dbl> <chr>          <chr>          
    ##  1     1 A     TRUE         77.3 tail           scales         
    ##  2     2 B     TRUE         49.1 tail           spikes         
    ##  3     3 C     TRUE         81.8 tail           fur            
    ##  4     4 D     TRUE        143.  teeth          fur            
    ##  5     5 E     TRUE        106.  tail           fur            
    ##  6     6 F     TRUE         95.2 tail           scales         
    ##  7     7 G     TRUE        104.  tail           spikes         
    ##  8     8 H     TRUE         90.2 teeth          scales         
    ##  9     9 I     TRUE         86.8 nose           spikes         
    ## 10    10 J     TRUE        115.  teeth          spikes         
    ## # ℹ 15 more rows

## plot tree with annotations

now we can use ggtree to use column data in the joined treedata object
for labels, colors, etc, etc

``` r
nwk_tree_with_info %>% 
    ggtree(aes(color=fake_phenotype)) +  # color here colors the branches
    geom_tiplab(aes(color=fake_phenotype)) + # color here colors the tip labels
    geom_tippoint(aes(subset= fake_phenotype=="tail"), color="black") # add dots to taxa, using subset
```

![](ggtree_demo_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

# add heatmap to the right side of a tree

gheatmap is the function for this

## explore gheatmap demo code from ggtree book

gheatmap demo code from the [ggtree online
book](http://yulab-smu.top/treedata-book/chapter7.html#gheatmap):

first read in their example tree, `beast_tree`:

``` r
beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
beast_tree <- read.beast(beast_file)
```

`beast_tree` is a treedata/tidytree object with 76 tips

next read in their example genotype data:

``` r
genotype_file <- system.file("examples/Genotype.txt", package="ggtree")
genotype <- read.table(genotype_file, sep="\t", stringsAsFactor=F)
colnames(genotype) <- sub("\\.$", "", colnames(genotype))
```

`genotype` is a data.frame with 76 rows and 8 columns

`?gheatmap` help page says that the heatmap data must be a matrix or
data.frame (not a tibble). I checked - the rownames of `genotype` are
identical to the tip labels of the tree, although they’re not in the
same order

``` r
# identical(rownames(genotype), beast_tree@phylo$tip.label)
# FALSE
# identical(sort(rownames(genotype)), sort(beast_tree@phylo$tip.label))
# TRUE
```

Now plot:

``` r
## save a ggtree plot object
p <- ggtree(beast_tree, mrsd="2013-01-01") + 
    geom_treescale(x=2008, y=1, offset=2) + 
    geom_tiplab(size=2)

## add the heatmap to the righthand side. 
gheatmap(p, genotype, offset=15, width=1.5, font.size=3,
         colnames_offset_y= -1) + 
    # scale_y_continuous makes sure we can see the colnames
    scale_y_continuous(expand=c(0.1, 0.1))
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

![](ggtree_demo_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Same thing but also:  
- add x axis scale bar (and heatmap colnames are less ugly now)  
- use a custom color scheme for the heatmap via `scale_fill_manual`.
Default behavior for NAs is that there is no rectangle plotted. If we
use `scale_fill_manual` the default behavior for NAs is now different -
they are dark gray

``` r
p <- ggtree(beast_tree, mrsd="2013-01-01") + 
    geom_tiplab(size=2, align=TRUE, linesize=.5) + 
    theme_tree2()
gheatmap(p, genotype, offset=8, width=0.6, 
         colnames=FALSE, legend_title="genotype") +
    scale_x_ggtree() + 
    scale_y_continuous(expand=c(0, 0.3)) +
    scale_fill_manual(breaks=c("HuH3N2", "pdm", "trig"), 
                      values=c("steelblue", "firebrick", "darkgreen"), 
                      name="genotype")
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.
    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](ggtree_demo_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

## gheatmap example with my fake tree and fake data

``` r
## imagine we're plotting amino acid changes, like Maria's trying to do
# here, NA represents "same as taxon A" (same as human in Maria's case)
# and we'll plot each amino acid as a separate color
tip_dat <- data.frame(row.names=nwk_tree_with_info@phylo$tip.label)
tip_dat$pos1 <- NA
tip_dat$pos2 <- NA
tip_dat$pos3 <- NA
tip_dat[6:8,"pos1"] <- "W"
tip_dat[4:13,"pos2"] <- "T"
tip_dat[9:12,"pos3"] <- "C"

## save a ggtree plot object
p <- nwk_tree_with_info %>% 
    ggtree(aes()) + 
    geom_tiplab(aes())

## add the heatmap to the righthand side. the heatmap data must be matrix or data.frame
# this is the original code
gheatmap(p, tip_dat, 
         offset=3,
         width=0.25, font.size=3) 
```

![](ggtree_demo_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

Or maybe we turn that into WT/ nonWT

``` r
tip_dat2 <- tip_dat %>% 
    mutate(across(everything(),
                  function(x) {
                      case_when( is.na(x) ~ "WT",
                                 TRUE ~ "nonWT")
                  }))
gheatmap(p, tip_dat2, 
         offset=3,
         width=0.25, font.size=3) +
    scale_fill_manual(breaks=c("WT", "nonWT"), 
                      values=c("lightgray", "firebrick"), 
                      name="")
```

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](ggtree_demo_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

It doesn’t work if tip_dat3 is a tibble, although it doesn’t give an
error, it just doesn’t give us a correct heatmap. Also, setting rownames
on a tibble is deprecated.

``` r
tip_dat3 <- tip_dat2 %>% 
    as_tibble(rownames=NA)
## we can check the rownames are present:
# rownames(tip_dat3)

gheatmap(p, tip_dat3, 
         offset=3,
         width=0.25, font.size=3) +
    scale_fill_manual(breaks=c("WT", "nonWT"), 
                      values=c("lightgray", "firebrick"), 
                      name="")
```

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](ggtree_demo_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->
