# Rtest_and_Rnotes

my R playground and notes

A lot of this is quite old, and most of it is only for myself, not meant for other people, although there are a few files in here I use for teaching.

rhino location: `~/FH_fast_storage/git_more_repos/Rtest_and_Rnotes`


## Very partial table of contents of this repo

There are many files in this repo, but these are some that will probably be useful frequently:

- [R tips and tricks](https://github.com/jayoung/Rtest_and_Rnotes/blob/main/notes/R_tips_and_tricks.md). In a bit of a random order right now.
- [ggplot tips and tricks](https://github.com/jayoung/Rtest_and_Rnotes/blob/main/Rscripts/ggplot/ggplot_tips_and_tricks.md)
- [ggtree_demo.md](https://github.com/jayoung/Rtest_and_Rnotes/blob/main/Rscripts/phylogenetics/ggtree_demo.md)
- [tutorials_JY_forOthers](https://github.com/jayoung/Rtest_and_Rnotes/tree/main/tutorials_JY_forOthers) folder has example code I used to teach others
- [render_from_command_line_NOTES.md](https://github.com/jayoung/Rtest_and_Rnotes/blob/main/notes/render_from_command_line_NOTES.md)
- [package_installation_NOTES.md](https://github.com/jayoung/Rtest_and_Rnotes/blob/main/notes/package_installation_NOTES.md)

Some notes I'm storing in other repos:

- [AI_notes.md](https://github.com/jayoung/MalikLab_bioinformaticsResources/blob/main/janets_NOTES_forMyself/programming_and_statistics/AI_notes.md)

## Quick things I often need to remember or copy/paste:

To Knit/render an Rmd doc from the linux command line, I have a couple of scripts in this repo (in `useful_functions`) that I can run from gizmo/rhino: `render_Rmd.perl` (run scripts in parallel) and `render_Rmd_series.perl` (run scripts in series) 

Insert an image into a .md or .Rmd document:
```
![Caption for the picture.](/path/to/image.png){#id .class width=50% height=50%}
```

Rotate x axis labels on a ggplot:
```
+ theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
```


## My coding intentions

Use `janitor::tabyl` more for cross-tables (see [demo](https://github.com/jayoung/Rtest_and_Rnotes/blob/main/Rscripts/janitor_cross_tables_demo.md))

Use native pipe (switch default in Rstudio)

Use renv for new projects

## Current learning to do list

see another list [here](https://github.com/jayoung/thoughts/blob/main/learning_to_do.md#r-learning-and-statistics)

[Tidyverse style guide](https://style.tidyverse.org/syntax.html).  I got up to the end of the 'syntax' section.  Perhaps see also the [Advanced R style guide](http://adv-r.had.co.nz/Style.html).

R stuff:
for file lists, check out `dir()` rather than `list.files()`:
```
files <- dir(here("data", "participants"), pattern="*.csv")
```

To read multiple files all at once, check out `purrr::map_df`:
```
data <- files |>
    map_df(~read_csv(file=here("data", "particants", .x)))
```


there's a package called 'readxl' - it's part of tidyverse, but not core tidyverse. Has a function called `readxl::read_excel()` with sheet option.

`furrr` R package is like `purrr`, but parallelized 

joins using plyranges R package - different types https://www.bioconductor.org/packages/devel/bioc/vignettes/plyranges/inst/doc/an-introduction.html#9_Learning_more

writing R packages and documentation - https://style.tidyverse.org/documentation.html
- R packages that help develop and maintain and test packages are  {devtools}, {testthat}, {usethis} 
- software design principals - rather than my own single bloated package, think about a universe of smaller packages. Functions and packages both benefit from being as modular as possible. https://milesmcbain.xyz/posts/data-analysis-reuse/
- R package writing beginners tips https://www.youtube.com/watch?v=F1GJSn9SqTk

R learning - data cleaning: [a primer](https://rfortherestofus.com/2019/12/how-to-clean-messy-data-in-r/), [some tips](https://rfortherestofus.com/2020/02/data-cleaning-tips-in-r/) and [a 15 min video](https://rfortherestofus.com/2020/01/recoding-numeric-values-to-character-values-automatically-in-r/)

check out these two packages for making interactive graphs: ggiraph and plotly. video tutorial https://rfortherestofus.com/2021/03/interactive-graphs-in-r/

## Resources

[R primers](https://r-primers.andrewheiss.com)

Nick Tierney's (mostly) [rstats blog](https://www.njtierney.com)

R For The Rest Of Us [resources](https://rfortherestofus.com/resources)

Advice on [making figures](https://github.com/MichaelClerx/making-figures/tree/main)

The R graph gallery(https://r-graph-gallery.com) and a [list of packages that extend ggplot](https://r-graph-gallery.com/best-dataviz-packages.html)

[r-packages.io](https://r-packages.io/) - searchable documentation of all CRAN packages

[Tutorial](https://www.jumpingrivers.com/blog/intro-to-theme-ggplot2-r/) on using `theme()` to change ggplot appearance/aesthetics.

Some notes on [good coding practices](https://r4ds.hadley.nz/workflow-scripts.html) - using Rmarkdown, clean environments, reproducibility, Rprojects

Rmarkdown:

* Rstudio's [intro to Rmarkdown](https://rmarkdown.rstudio.com/lesson-1.html)    
* [intro2r chapter 8](https://intro2r.com/rmarkdown_r.html)    
* detailed [Rmarkdown guide](https://bookdown.org/yihui/rmarkdown/)    

miscellaneous coll-looking R tips from [Luke Pemberton](https://lpembleton.rbind.io/ramblings/R/), like embedding smaller plots as insets on top of bigger ones, including colors in titles, nice axis formatting, etc, etc

Blog post on working with [regular expressions](https://blog.djnavarro.net/posts/2024-12-16_regex-backreferences/)

Bioconductor support:

- [https://support.bioconductor.org/](https://support.bioconductor.org/)
- Zulip workspace (a bit like Slack)

[CRAN task views](http://cran.r-project.org/web/views/) may help identify the right package(s) to use for a certain task. See also [Bioconductor workflow vignettes](https://www.bioconductor.org/packages/release/BiocViews.html#___Workflow)

Small daily lessons:
- [25-day advent calendaR](https://github.com/kiirsti/ggplot_adventcalendaR) providing an introduction to ggplot2
- The 25 days of chRistmas: [an R advent calendaR](https://kiirstio.wixsite.com/kowen/post/the-25-days-of-christmas-an-r-advent-calendar)

### haven't looked yet: future

[The Elements of Data Analytic Style](https://leanpub.com/datastyle) by Jeff Leek (a Leanpub book)

Data Science resources [list](https://nrennie.rbind.io/data-science-resources/)




## List of useful packages

### patchwork to combine plots

patchwork package is great. See [demo code](https://github.com/jayoung/Rtest_and_Rnotes/blob/main/Rscripts/patchwork_package_demo_and_tips.md).

### more control over axes and legends

[legendry](https://teunbrand.github.io/legendry/articles/keys.html) looks useful

### importing images and combining with R plots

here's how you'd combine an imported image (import using magick package) with a ggplot:

```{r}
library(patchwork)
library(ggplot2)
library(magick)

plt1 <- image_read("https://bellard.org/bpg/2.png") |>
  image_ggplot()
plt2 <- iris |> 
  ggplot(aes(x=Sepal.Length, y=Sepal.Width)) +
  geom_point()

plt1 | plt2
```

### pretty tables in Rmarkdown (etc)

`kable/kableExtra`

`flextable` - see [Rscripts/table_display_flextable_demo.md](Rscripts/table_display_flextable_demo.md) 

[`reactable`](https://glin.github.io/reactable/index.html)

`emphatic` - see [Rscripts/emphatic_demo.md](Rscripts/emphatic_demo.md) 

[`tinytable`](https://vincentarelbundock.github.io/tinytable/)

`gt` - see [Rscripts/table_display_gt_demo.md](Rscripts/table_display_gt_demo.md). 

- It can make [multicolumn tables](https://viz.aweatherman.com/viz/538-caption/plot.html), i.e. can wrap a very long table.  
- That same [tutorial](https://viz.aweatherman.com/viz/538-caption/plot.html) shows how to make a multicolumn table, how to include little logos within each cell, and how to make a nice-looking two-part footnote.    
- And [here](https://rfortherestofus.com/2023/10/ggplots-in-gt-tables) is a tutorial to show how to use `gt` to put little graphs in some cells of a table

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

Maybe try [msa package](https://bioconductor.org/packages/devel/bioc/vignettes/msa/inst/doc/msa.pdf)

### phylogenetics

ape

[ggtree](https://yulab-smu.top/treedata-book/)  (also tidytree and treeio) (can parse PAML and Hyphy output as well as make some very nice plots).   [ggtree publication](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12628)

[ggtreeExtra](https://bioconductor.org/packages/devel/bioc/vignettes/ggtreeExtra/inst/doc/ggtreeExtra.html).  ggtree can use geom_facet to align associated graphs to the tree but it only works with rectangular, roundrect, ellipse and slanted layouts. ggtreeExtra allows graphs on a tree in rectangular, circular, fan and radial layouts

in ~/domesticated_capsid/Rreports/RTL3_frameshift_plots_v2_aln28.Rmd I got a tree of >5000 mammal species from Upham publication, and extracted the species I want

Perhaps the RERconverge package "for associating evolutionary rates with convergent traits"

### orthology

[Orthology.eg.db](https://bioconductor.org/packages/release/data/annotation/html/Orthology.eg.db.html). Uses NCBI orthology data. Tested it a bit April 2025. Not useful for human-mosquito (Aedes aegypti) but probably OK for human-mouse. See `~/FH_fast_storage/paml_screen/Priya_Shah_YFV/Rreports/orthology_test.Rmd`

Inparanoid downloads

### plotting genes etc:

[gggenes](https://cran.r-project.org/web/packages/gggenes/vignettes/introduction-to-gggenes.html) - demo shows it only plotting gene arrows, not other data.  Also suggests: `gggenomes` for visualising comparative genomics, `plasmapR` for quickly drawing plasmid maps from GenBank files

[ggbio](https://www.bioconductor.org/packages/release/bioc/html/ggbio.html) - not sure whether it is still maintained  

[ggcoverage](https://github.com/showteeth/ggcoverage)

[Gviz](https://www.bioconductor.org/packages/release/bioc/vignettes/Gviz/inst/doc/Gviz.html) is what I used for the tetrahymena project, and Michelle's project, and SATAY data. Seems to be maintained and very functional.

[GenVisR](https://bioconductor.org/packages/release/bioc/html/GenVisR.html)

for genomes with karyotypes: [karyoploteR](https://bioconductor.org/packages/release/bioc/vignettes/karyoploteR/inst/doc/karyoploteR.html). Perhaps also try `GenomeGraphs` or `HilbertVis`

rtracklayer can make plots by interacting with a UCSC browser

igvR can interact with IGV

[tidyGenomeBrowser](https://github.com/MalteThodberg/tidyGenomeBrowser)

[GenomicPlot](https://bioconductor.org/packages/release/bioc/vignettes/GenomicPlot/inst/doc/GenomicPlot_vignettes.html) is more for making metaplots combining data over multiple features

Explored a few options in April 2024 for the SATAY data - see  [browser_style_plots_failed_attempts.Rmd](https://github.com/gak2882/SATAY/blob/main/janet_Rscripts/browser_style_plots_failed_attempts.md) and [browser_style_plots.Rmd](https://github.com/gak2882/SATAY/blob/main/janet_Rscripts/browser_style_plots_demo.md)

[genoPlotR](https://genoplotr.r-forge.r-project.org/screenshots.php)

### wordclouds

`wordcloud` and `wordcloud2` packages.  see `/fh/fast/malik_h/user/jayoung/presentations/MalikLab/otherSlides_mine/KennedyHighSchoolVisit_2021_Dec7/Hutch_wordCloud.R`

### violin plots

use `ggplot - geom_violin()`. 

Some other options are vioplot::vioplot(), DescTools::PlotViolin(),  easyGgplot2::ggplot2.violinplot(),  UsingR::violinplot().  
Before the days of ggplot, I noted that I like PlotViolin better than vioplot, but when I run it on large datasets it is very slow if I allow it to use its default bandwidth selection algorithm.   If I specify the bw="nrd0" option, it is MUCH quicker. 

See also [here](http://www.sthda.com/english/wiki/ggplot2-violin-plot-easy-function-for-data-visualization-using-ggplot2-and-r-software)

### other packages

`flowchart` and `ggflowchart` packages

`ggarrow` and `arrowheadr` packages for nicer looking arrows

[`gcplyr` package](https://mikeblazanin.github.io/gcplyr/) for microbial growth curves. Can help read platereader data (with metadata) in and get it in a tidy format. Can model various parameters of growth curves, "like growth rate/doubling time, maximum density (carrying capacity), lag time, area under the curve, diauxic shifts, extinction, and more without fitting an equation for growth to your data."

### ggplot themes

https://rfortherestofus.com/2019/08/themes-to-improve-your-ggplot-figures/


### dumbbell plots

[A nice alternative](https://albert-rapp.de/posts/ggplot2-tips/15_alternative_paired_bars/15_alternative_paired_bars) to paired bar charts. That tutorial shows other alternatives too: arrow plots, slope plots. 

