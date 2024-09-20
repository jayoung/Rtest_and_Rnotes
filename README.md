# Rtest_and_Rnotes

my R playground and notes

A lot of this is quite old, and most of it is only for myself, not meant for other people, although there are a few files in here I use for teaching.

rhino location: `~/FH_fast_storage/git_more_repos/Rtest_and_Rnotes`


## Current learning
[Tidyverse style guide](https://style.tidyverse.org/syntax.html).  I got up to the end of the 'syntax' section.  Perhaps see also the [Advanced R style guide](http://adv-r.had.co.nz/Style.html).

## Resources

[R primers](https://r-primers.andrewheiss.com)

Nick Tierney's (mostly) [rstats blog](https://www.njtierney.com)

R For The Rest Of Us [resources](https://rfortherestofus.com/resources)

Advice on [making figures](https://github.com/MichaelClerx/making-figures/tree/main)

Some notes on [good coding practices](https://r4ds.hadley.nz/workflow-scripts.html) - using Rmarkdown, clean environments, reproducibility, Rprojects

Rmarkdown:

* Rstudio's [intro to Rmarkdown](https://rmarkdown.rstudio.com/lesson-1.html)    
* [intro2r chapter 8](https://intro2r.com/rmarkdown_r.html)    
* detailed [Rmarkdown guide](https://bookdown.org/yihui/rmarkdown/)    

miscellaneous coll-looking R tips from [Luke Pemberton](https://lpembleton.rbind.io/ramblings/R/), like embedding smaller plots as insets on top of bigger ones, including colors in titles, nice axis formatting, etc, etc

### haven't looked yet: future

[The Elements of Data Analytic Style](https://leanpub.com/datastyle) by Jeff Leek (a Leanpub book)


## Things I've learned


### R does rounding weirdly!

See [here](https://psiaims.github.io/CAMIS/R/rounding.html)

"The `round()` function in Base R will round to the nearest whole number and ‘rounding to the even number’ when equidistant, meaning that exactly 12.5 rounds to the integer 12. Note that the janitor package in R contains a function `round_half_up()` that rounds away from zero. in this case it rounds to the nearest whole number and ‘away from zero’ or ‘rounding up’ when equidistant, meaning that exactly 12.5 rounds to the integer 13.""

### Rstudio tricks

"Reindenting your code only shifts things around horizontally. If you want more powerful code reformatting, try using “Code > Reformat Code” (or use ⌘⇧A on macOS or ctrl + shift + A on Windows). It’s a more aggressive form of reformatting that will add extra line breaks and other things to make the code more readable."

#### Snippets  

Example: type `fun` and press the `tab` key, and R provides the skeleton of a new function

To see all snippets:  Tools - Edit Code Snippets


### Debugging

Three options:    
- `browser()` (place inside a function, temporarily)    
- `debug(myFunction)` plus `undebug(myFunction)`    
- `debugonce()`    
See (`explore_debugging_functions.R`)[Rscripts/explore_debugging_functions.R] for details.

### Miscellaneous 

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

On a Mac, R packages go here - `/Library/Frameworks/R.framework/Versions` - in subdirectories by version. After installing new R, can delete old packages to save disk space


### combining plots

patchwork package is great

### importing images and combining with R plots

here's how you'd combine an imported image (import using magick package) with a ggplot:

```{r}
library(patchwork)
library(ggplot2)
library(magick)

plt1 <- image_read("https://bellard.org/bpg/2.png") %>%
  image_ggplot()
plt2 <- iris %>% 
  ggplot(aes(x=Sepal.Length, y=Sepal.Width)) +
  geom_point()

plt1 | plt2
```

### pretty tables in Rmarkdown (etc)

`kable/kableExtra`

`flextable` - see [Rscripts/flextable_demo.md](Rscripts/flextable_demo.md) 

[`reactable`](https://glin.github.io/reactable/index.html)

`emphatic` - see [Rscripts/emphatic_demo.md](Rscripts/emphatic_demo.md) 

[`tinytable`](https://vincentarelbundock.github.io/tinytable/)

`gt` can make [multicolumn tables](https://viz.aweatherman.com/viz/538-caption/plot.html), i.e. can wrap a very long table.  That same [tutorial](https://viz.aweatherman.com/viz/538-caption/plot.html) shows how to make a multicolumn table, how to include little logos withiin each cell, and how to make a nice-looking two-part footnote. 

(and I think some others)

### Volcano plots

[EnhancedVolcano](https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html)

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

[Gviz](https://www.bioconductor.org/packages/release/bioc/vignettes/Gviz/inst/doc/Gviz.html) is what I used for the tetrahymena project, and Michelle's project, and SATAY data. Seems to be maintained and very functional.

[GenVisR](https://bioconductor.org/packages/release/bioc/html/GenVisR.html)

for genomes with karyotypes: [karyoploteR](https://bioconductor.org/packages/release/bioc/vignettes/karyoploteR/inst/doc/karyoploteR.html)

rtracklayer can make plots by interacting with a UCSC browser

igvR can interact with IGV

[tidyGenomeBrowser](https://github.com/MalteThodberg/tidyGenomeBrowser)

[GenomicPlot](https://bioconductor.org/packages/release/bioc/vignettes/GenomicPlot/inst/doc/GenomicPlot_vignettes.html) is more for making metaplots combining data over multiple features

Explored a few options in April 2024 for the SATAY data - see ~/FH_fast_storage/forOtherPeople/forGrantKing/SATAY/janet_Rscripts/ files browser_style_plots_failed_attempts.Rmd and browser_style_plots.Rmd


#### wordclouds

`wordcloud` and `wordcloud2` packages.  see ``/Volumes/malik_h/user/jayoung/presentations/MalikLab/otherSlides_mine/KennedyHighSchoolVisit_2021_Dec7/Hutch_wordCloud.R`

#### violin plots

use `ggplot - geom_violin()`. 

Some other options are vioplot::vioplot(), DescTools::PlotViolin(),  easyGgplot2::ggplot2.violinplot(),  UsingR::violinplot().  

Before the days of ggplot, I noted that I like PlotViolin better than vioplot, but when I run it on large datasets it is very slow if I allow it to use its default bandwidth selection algorithm.   If I specify the bw="nrd0" option, it is MUCH quicker. 

See also [here](http://www.sthda.com/english/wiki/ggplot2-violin-plot-easy-function-for-data-visualization-using-ggplot2-and-r-software)

### other packages

`flowchart` and `ggflowchart` packages

`ggarrow` and `arrowheadr` packages for nicer looking arrows

#### ggplot themes

https://rfortherestofus.com/2019/08/themes-to-improve-your-ggplot-figures/
