# Rtest_and_Rnotes
my R playground and notes

rhino location: `~/FH_fast_storage/git_more_repos/Rtest_and_Rnotes`


## Current learning
[Tidyverse style guide](https://style.tidyverse.org/syntax.html).  I got up to the end of the 'syntax' section.  Perhaps see also the [Advanced R style guide](http://adv-r.had.co.nz/Style.html).

## Resources

Nick Tierney's (mostly) [rstats blog](https://www.njtierney.com)

R For The Rest Of Us [resources](https://rfortherestofus.com/resources)

### haven't looked yet: future

[The Elements of Data Analytic Style](https://leanpub.com/datastyle) by Jeff Leek (a Leanpub book)


## Things I've learned


### Rstudio tricks

"Reindenting your code only shifts things around horizontally. If you want more powerful code reformatting, try using “Code > Reformat Code” (or use ⌘⇧A on macOS or ctrl + shift + A on Windows). It’s a more aggressive form of reformatting that will add extra line breaks and other things to make the code more readable."

### Debugging

Three options:    
- `browser()` (place inside a function, temporarily)    
- `debug(myFunction)` plus `undebug(myFunction)`    
- `debugonce()`    
See (`explore_debugging_functions.R`)[Rscripts/explore_debugging_functions.R] for details.

### Misc 

The 'embracing' operator (`{{ }}`), and unquoting using !! and !!! - see [`testCode.R`](Rscripts/testCode.R) for details.

& versus && (and | versus ||):  use the short form for bitwise operation on vectors. Use the long form when we want a single TRUE/FALSE answer.   `any` and `all` functions run OR and AND on all elements of a vector

```
x <- c(TRUE,TRUE,FALSE,FALSE)
y <- c(TRUE,FALSE,TRUE,FALSE)
x & y
# x && y # this is no good!
any(x)  ## TRUE
all(x)  ## FALSE
```

The `switch` function - a multiway `if` statement, I think?
```
centre <- function(x, type) {
  switch(type,
         mean = mean(x),
         median = median(x),
         trimmed = mean(x, trim = .1))
}
x <- rcauchy(10)
centre(x, "mean")
centre(x, "median")
centre(x, "trimmed")
```

In `switch`, if there's an empty argument, it 'falls-through' to the next thing (e.g. here, `myFunc("a")` returns the same thing as `myFunc("b")`).
Not also that we can add `call. = FALSE` to a `stop` statement to modify the error message that'll be produced
```
myFunc <- function(x) {
  switch(x, 
    a = ,
    b = 1, 
    c = 2,
    stop("Unknown `x`", call. = FALSE)
  )
}
```


## Useful packages

### pretty tables in Rmarkdown (etc)

`kable/kableExtra`

`flextable` - see [Rscripts/flextable_demo.md](Rscripts/flextable_demo.md) 

(and I think some others)


### fastq files

shortRead package

[ngsReports package](https://www.bioconductor.org/packages/devel/bioc/vignettes/ngsReports/inst/doc/ngsReportsIntroduction.html#loading-fastqc-data-into-r) - ~/public_databases/NCBI/SRA/data/mammalian_expression_profiles/human/human_spermMicrobiomeTotalRNAseq/fastqc/parseFastqc.R


### multiple sequence alignments (MSAs)

a [list of tools](https://cmdcolin.github.io/awesome-genome-visualization/?latest=true&tag=MSA
) for viewing MSAs

viewing MSA alongside a tree: try ggtree() with option msaplot(). also "ggtreeExtra() has a different way to do it which is probably more flexible"


### phylogenetics

ape

[ggtree](https://yulab-smu.top/treedata-book/)  (also tidytree and treeio) (can parse PAML and Hyphy output as well as make some very nice plots).   [ggtree publication](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12628)

[ggtreeExtra](https://bioconductor.org/packages/devel/bioc/vignettes/ggtreeExtra/inst/doc/ggtreeExtra.html).  ggtree can use geom_facet to align associated graphs to the tree but it only works with rectangular, roundrect, ellipse and slanted layouts. ggtreeExtra allows graphs on a tree in rectangular, circular, fan and radial layouts

in ~/domesticated_capsid/Rreports/RTL3_frameshift_plots_v2_aln28.Rmd I got a tree of >5000 mammal species from Upham publication, and extracted the species I want

### plotting genes etc:

[gggenes](https://cran.r-project.org/web/packages/gggenes/vignettes/introduction-to-gggenes.html) - demo shows it only plotting gene arrows, not other data.  Also suggests: `gggenomes` for visualising comparative genomics, `plasmapR` for quickly drawing plasmid maps from GenBank files

[ggbio](https://www.bioconductor.org/packages/release/bioc/html/ggbio.html) - not sure whether it is still maintained  

[ggcoverage](https://github.com/showteeth/ggcoverage)

[Gviz](https://www.bioconductor.org/packages/release/bioc/vignettes/Gviz/inst/doc/Gviz.html) (that's what I used for the tetrahymena project, and Michelle's project, and SATAY data)

[GenVisR](https://bioconductor.org/packages/release/bioc/html/GenVisR.html)

for genomes with karyotypes: [karyoploteR](https://bioconductor.org/packages/release/bioc/vignettes/karyoploteR/inst/doc/karyoploteR.html)

rtracklayer can make plots by interacting with a UCSC browser

igvR can interact with IGV

[tidyGenomeBrowser](https://github.com/MalteThodberg/tidyGenomeBrowser)

[GenomicPlot](https://bioconductor.org/packages/release/bioc/vignettes/GenomicPlot/inst/doc/GenomicPlot_vignettes.html) is more for making metaplots combining data over multiple features

Explored a few options in April 2024 for the SATAY data - see ~/FH_fast_storage/forOtherPeople/forGrantKing/SATAY/janet_Rscripts/ files browser_style_plots_failed_attempts.Rmd and browser_style_plots.Rmd


### miscellaneous

wordclouds - `wordcloud` and `wordcloud2` packages.  see ``/Volumes/malik_h/user/jayoung/presentations/MalikLab/otherSlides_mine/KennedyHighSchoolVisit_2021_Dec7/Hutch_wordCloud.R`

violin plots - will probably use `ggplot - geom_violin()`. Some other options are vioplot::vioplot(), DescTools::PlotViolin(),  easyGgplot2::ggplot2.violinplot(),  UsingR::violinplot().  I noted a while back that I like PlotViolin better than vioplot, but when I run it on large datasets it is very slow if I allow it to use its default bandwidth selection algorithm.   If I specify the bw="nrd0" option, it is MUCH quicker. See also [here](http://www.sthda.com/english/wiki/ggplot2-violin-plot-easy-function-for-data-visualization-using-ggplot2-and-r-software)

