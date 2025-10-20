ggtree_demo
================
Janet Young

2025-10-20

(this doesn’t work on the fhR 4.2.0 I could run via the Hutch Rstudio
server, with the package versions I have there right now. But fhR 4.4.0
is fine)

# Useful links, general advice

The [ggtree manual](http://yulab-smu.top/treedata-book/index.html) is
great (very detailed). - Chapter 1 looks at reading/writing trees. The
treeio and ape packages both have tree read/write functions).

There is a [ggtree google
group](https://groups.google.com/g/bioc-ggtree) for questions.

There is also a [ggtreeExtra
package](https://www.bioconductor.org/packages/release/bioc/vignettes/ggtreeExtra/inst/doc/ggtreeExtra.html)
that provides a function, `geom_fruit`, to align graphs to a tree,
perhaps more easily than with ggtree.

The [ape package](https://cran.r-project.org/web/packages/ape/) is also
useful for phylogenetics

# Read a tree and plot it

Read in a newick tree file:

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
    ## Rooted; includes branch length(s).

Plot it using ape’s basic plot function

``` r
plot(nwk_tree)
```

![](ggtree_demo_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

Plot it using ggtree

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
    ##  1 A            59.8 tail           spikes         
    ##  2 B           119.  nose           spikes         
    ##  3 C           108.  teeth          scales         
    ##  4 D           108.  tail           spikes         
    ##  5 E            85.5 nose           spikes         
    ##  6 F            50.7 tail           spikes         
    ##  7 G           101.  teeth          fur            
    ##  8 H            89.7 teeth          spikes         
    ##  9 I           122.  tail           fur            
    ## 10 J           119.  nose           spikes         
    ## 11 K            81.5 tail           fur            
    ## 12 L            87.2 tail           fur            
    ## 13 M           146.  nose           fur

## Combine tree and tbl

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
    ## Rooted; includes branch length(s).
    ## 
    ## with the following features available:
    ##   '', 'fake_height', 'fake_phenotype', 'fake_phenotype2'.
    ## 
    ## # The associated data tibble abstraction: 25 × 6
    ## # The 'node', 'label' and 'isTip' are from the phylo tree.
    ##     node label isTip fake_height fake_phenotype fake_phenotype2
    ##    <int> <chr> <lgl>       <dbl> <chr>          <chr>          
    ##  1     1 A     TRUE         59.8 tail           spikes         
    ##  2     2 B     TRUE        119.  nose           spikes         
    ##  3     3 C     TRUE        108.  teeth          scales         
    ##  4     4 D     TRUE        108.  tail           spikes         
    ##  5     5 E     TRUE         85.5 nose           spikes         
    ##  6     6 F     TRUE         50.7 tail           spikes         
    ##  7     7 G     TRUE        101.  teeth          fur            
    ##  8     8 H     TRUE         89.7 teeth          spikes         
    ##  9     9 I     TRUE        122.  tail           fur            
    ## 10    10 J     TRUE        119.  nose           spikes         
    ## # ℹ 15 more rows

## Plot tree with annotations

now we can use ggtree to use column data in the joined treedata object
for labels, colors, etc, etc

``` r
x <- "tail"
nwk_tree_with_info %>% 
    ggtree(aes(color=fake_phenotype)) +  # color here colors the branches
    geom_tiplab(aes(color=fake_phenotype))  + # color here colors the tip labels
    ## there's some weird bug where I got errors if I tried to do fake_phenotype=="tail" but I get around it like this
    geom_tippoint(aes(subset= fake_phenotype==.env$x),
                  color="black") # add dots to taxa, using subset
```

![](ggtree_demo_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

# Add heatmap to the right side of a tree

gheatmap is the function for this

## Explore gheatmap demo code from ggtree book

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

Show that first plot in circle layout (ignore branch lengths)

``` r
## save a ggtree plot object
p <- ggtree(beast_tree, mrsd="2013-01-01", 
            layout="circular",
            branch.length="none") 

## add the heatmap 
gheatmap(p, genotype, 
         offset=1, 
         width=1, font.size=2,
         colnames_offset_y= -1,
         colnames_angle=270) +
    # scale_y_continuous makes sure we can see the colnames
    scale_y_continuous(expand=c(0.01, 0.01))
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

![](ggtree_demo_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

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

![](ggtree_demo_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

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

![](ggtree_demo_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

Note that it doesn’t work if tip_dat3 is a tibble (see plot below, not
useful), although it doesn’t give an error, it just doesn’t give us a
correct heatmap. Also, setting rownames on a tibble is deprecated.

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

![](ggtree_demo_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

# Tree plus data to the right

geom_point to the right

``` r
tr <- rtree(10)
dd <- tibble(id=tr$tip.label) %>% 
    mutate(value=as.integer(str_remove(id, "t")))
```

``` r
p <- ggtree(tr) +
    geom_tiplab() +
    geom_facet(panel = "Data", 
               data = dd,  ## the first column of dd, no matter what it's called, is where ggtree looks for the taxon labels
               geom = geom_point, mapping = aes(x = value)) + 
    theme_light() +   # need theme_bw or similar to get the x axis scales to show
    xlim_tree(10) + # alter tree x scale
    xlim_expand(c(0, 13), "Data") # alter points x scale
p
```

![](ggtree_demo_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

geom_col to the right

``` r
p <- ggtree(tr) +
    geom_tiplab() +
    geom_facet(panel = "Data", data = dd, 
               geom = geom_col, mapping = aes(x = value),
               orientation = "y") + 
    theme_light() +   # need theme_bw or similar to get the x axis scales to show
    xlim_tree(10) + # alter tree x scale
    xlim_expand(c(0, 3), "Data") # alter col x scale
p
```

![](ggtree_demo_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

It CAN handle missing data:

``` r
dd_slice <- slice_head(dd, n=5)

ggtree(tr) +
    geom_tiplab() +
    geom_facet(panel = "Data", data = dd_slice, 
               geom = geom_col, mapping = aes(x = value),
               orientation = "y") + 
    theme_light() +   # need theme_bw or similar to get the x axis scales to show
    xlim_tree(10) + # alter tree x scale
    xlim_expand(c(0, 3), "Data") # alter col x scale
```

![](ggtree_demo_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

Simple geom_fruit demo:

``` r
ggtree(tr) +
    geom_tiplab() +
    geom_fruit(
        data=dd, 
        geom=geom_col, 
        mapping=aes(x=value, y=id)
    ) 
```

![](ggtree_demo_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

Now control how it looks more

``` r
ggtree(tr) +
    geom_tiplab() +
    geom_fruit(
        data=dd, 
        geom=geom_col, 
        mapping=aes(x=value, y=id),
        pwidth=0.8, ## pwidth affects width of the geom_col plot, relative to tree width (default is 0.2)
        offset=0.1, ## add some space between tree and plot
        axis.params = list(
            axis = "x",
            text.size = 4,
            vjust = 1, hjust = 0.5,
            limits = c(0, 10),
            line.size=0.5, line.color = "black",
            title="value", title.size=6
        )
    ) 
```

![](ggtree_demo_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

# Turning a tibble back into a tree

Sometimes we take a tree, use `as_tibble()` (it becomes a `tbl_tree`
object), do a bunch of manipulation (maybe it becomes a classic `tibble`
object).

It’s not totally obvious how to turn it back into a tree. If we use
`as.treedata()` on a `tbl_tree` object , it works fine. But if we use it
on a pure tibble, the branch lengths are lost unless we do
`as.treedata(branch.length, label)` instead.

Demo code for that:

``` r
x <- "(((Strix_aluco:8.2,Asio_otus:4.2):3.1,Athene_noctua:7.3):6.3,Tyto_alba:13.5);"
tree.owls <- read.tree(text= x)
tree.owls <- as.treedata(tree.owls)
```

Do I lose branch length on coercion? no, I don’t. `as_tibble()` actually
returns a `tbl_tree` object. It’s only for pure tibbles that we need to
do something more to ensure we retain branch lengths:
`as.treedata(branch.length, label)`.

``` r
# Do I lose branch length on coercion? no, I don't
tree.owls %>% 
    as_tibble() %>% 
    as.treedata() %>% 
    ggtree() +
    geom_tiplab() +
    geom_treescale() +
    hexpand(0.3) +
    labs(title="tree.owls coerced")
```

![](ggtree_demo_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
(tree.owls %>% 
        as_tibble() %>% 
        as.treedata())@phylo$edge.length
```

    ## [1]  8.2  4.2  7.3 13.5  6.3  3.1

This double coercion is a cheating way to get a pure tibble rather than
a `tbl_tree` object. Now we lose branch lengths:

``` r
tree.owls %>% 
    as_tibble() %>% 
    as_tibble() %>% 
    as.treedata() %>% 
    ggtree() +
    geom_tiplab() +
    geom_treescale() +
    hexpand(0.3) +
    labs(title="tree.owls coerced")
```

![](ggtree_demo_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

But we can avoid losing branch lengths by supplying a couple of
arguments to `as.treedata()`:

``` r
tree.owls %>% 
    as_tibble() %>% 
    as_tibble() %>% 
    as.treedata(branch.length, label) %>% 
    ggtree() +
    geom_tiplab() +
    geom_treescale() +
    hexpand(0.3) +
    labs(title="tree.owls coerced")
```

![](ggtree_demo_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] tidytree_0.4.6     treeio_1.32.0      ggtreeExtra_1.18.0 ggtree_3.16.3     
    ##  [5] ape_5.8-1          lubridate_1.9.4    forcats_1.0.0      stringr_1.5.2     
    ##  [9] dplyr_1.1.4        purrr_1.1.0        readr_2.1.5        tidyr_1.3.1       
    ## [13] tibble_3.3.0       ggplot2_3.5.2      tidyverse_2.0.0   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] yulab.utils_0.2.1  utf8_1.2.6         rappdirs_0.3.3     generics_0.1.4    
    ##  [5] ggplotify_0.1.3    stringi_1.8.7      lattice_0.22-7     hms_1.1.3         
    ##  [9] digest_0.6.37      magrittr_2.0.4     evaluate_1.0.5     grid_4.5.1        
    ## [13] timechange_0.3.0   RColorBrewer_1.1-3 fastmap_1.2.0      jsonlite_2.0.0    
    ## [17] ggnewscale_0.5.2   aplot_0.2.9        scales_1.4.0       lazyeval_0.2.2    
    ## [21] cli_3.6.5          rlang_1.1.6        withr_3.0.2        yaml_2.3.10       
    ## [25] tools_4.5.1        parallel_4.5.1     tzdb_0.5.0         gridGraphics_0.5-1
    ## [29] vctrs_0.6.5        R6_2.6.1           lifecycle_1.0.4    fs_1.6.6          
    ## [33] ggfun_0.2.0        pkgconfig_2.0.3    pillar_1.11.1      gtable_0.3.6      
    ## [37] glue_1.8.0         Rcpp_1.1.0         xfun_0.53          tidyselect_1.2.1  
    ## [41] rstudioapi_0.17.1  knitr_1.50         farver_2.1.2       patchwork_1.3.2   
    ## [45] htmltools_0.5.8.1  nlme_3.1-168       labeling_0.4.3     rmarkdown_2.29    
    ## [49] compiler_4.5.1
