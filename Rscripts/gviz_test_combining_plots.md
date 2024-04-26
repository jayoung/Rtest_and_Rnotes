gviz_test_combining_plots
================
Janet Young

2024-04-15

Code is from
[here](https://www.bioconductor.org/packages/release/bioc/vignettes/Gviz/inst/doc/Gviz.html#7_Composite_plots_for_multiple_chromosomes)

``` r
#
chr <- as.character(unique(seqnames(cpgIslands))) #chr7
gen <- genome(cpgIslands) # hg19
atrack <- AnnotationTrack(cpgIslands, name = "CpG")
plotTracks(list(atrack))
```

![](gviz_test_combining_plots_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->
Test using viewports to combine plots - it works

``` r
chroms <- c("chr1", "chr2", "chr3", "chr4")
itracks <- lapply(chroms, function(chr) {
    IdeogramTrack(genome = "mm9", chromosome = chr)
} )

## AnnotationTrack - gr1 - 4 regions of width 100bp-1kb, with strand but not score, on separate chromosomes
gr1 <- GRanges(seqnames = chroms, 
               ranges = IRanges(start = 1,  width = c(100, 400, 200,1000)),
               strand = c("+", "+", "-", "+"))
maTrack <- AnnotationTrack(range=gr1, genome = "mm9", 
                           chromosome = "chr1", name = "foo")

## DataTrack - gr2 - 170 regions on all 4 chromosomes, all 9bp wide, with values column but no score, all within the 4 regions in ggr1
gr2 <- GRanges(seqnames = rep(chroms, c(10, 40, 20, 100)),
               ranges = IRanges(start = c(seq(1, 100, len = 10),
                                          seq(1, 400, len = 40), 
                                          seq(1, 200, len = 20),
                                          seq(1, 1000, len = 100)), 
                                width = 9), values = runif(170))
mdTrack <- DataTrack(
    range = gr2,
    data = "values", chromosome = "chr1", genome = "mm9", name = "bar")

## axis track - just a little scale bar
mgTrack <- GenomeAxisTrack(scale = 50, labelPos = "below", exponent = 3)

## tile 4 plots, one per item in chroms
ncols <- 2
nrows <- length(chroms) %/% ncols
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrows, ncols)))
for(i in seq_along(chroms)) {
    pushViewport(viewport(layout.pos.col = ((i - 1) %% ncols) + 1,
                          layout.pos.row = (((i) - 1) %/% ncols) + 1))
    plotTracks(list(itracks[[i]], maTrack, mdTrack, mgTrack), 
               chromosome = chroms[i], add = TRUE, title=chroms[i])
    popViewport(1)
}
```

![](gviz_test_combining_plots_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->