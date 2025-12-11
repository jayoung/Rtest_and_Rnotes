# Rtest_and_Rnotes

my R playground and notes

A lot of this is quite old, and most of it is only for myself, not meant for other people, although there are a few files in here I use for teaching.

rhino location: `~/FH_fast_storage/git_more_repos/Rtest_and_Rnotes`


## Links, and frequently used notes

Many files in this repo, but these are some that will probably be useful frequently:

- [ggplot tips and tricks](https://github.com/jayoung/Rtest_and_Rnotes/blob/main/Rscripts/ggplot/ggplot_tips_and_tricks.md)
- [ggtree_demo.md](https://github.com/jayoung/Rtest_and_Rnotes/blob/main/Rscripts/phylogenetics/ggtree_demo.md)

Some notes I'm storing in other repos:

- [AI_notes.md](https://github.com/jayoung/MalikLab_bioinformaticsResources/blob/main/janets_NOTES_forMyself/programming_and_statistics/AI_notes.md)

## Intentions

Use `janitor::tabyl` more for cross-tables (see [demo](https://github.com/jayoung/Rtest_and_Rnotes/blob/main/Rscripts/janitor_cross_tables_demo.md))

Use native pipe (switch default in Rstudio)

Use renv for new projects

## Current learning to do list

see another list [here](https://github.com/jayoung/thoughts/blob/main/learning_to_do.md#r-learning-and-statistics)

[Tidyverse style guide](https://style.tidyverse.org/syntax.html).  I got up to the end of the 'syntax' section.  Perhaps see also the [Advanced R style guide](http://adv-r.had.co.nz/Style.html).



R stuff:
for file lists, check out `dir()`:
```
files <- dir(here("data", "participants"), pattern="*.csv")
```
for reading multiple files, check out `purrr::map_df`:
```
data <- files |>
    map_df(~read_csv(file=here("data", "particants", .x)))
```

downloading files within R:
```
download.file(url, destfile="data-raw/name-of-file.xlsx")
```


there's a package called 'readxl' part of tidyverse, but not core tidyverse. readxl::read_excel(). Has a sheet option

`furrr` R package is like `purrr`, but parallelized 


figure out gheatmap in ggtree


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

### haven't looked yet: future

[The Elements of Data Analytic Style](https://leanpub.com/datastyle) by Jeff Leek (a Leanpub book)

Data Science resources [list](https://nrennie.rbind.io/data-science-resources/)

## Positron

An IDE for R or python.  Perhaps a future replacement for Rstudio but as of 2025 it's not able to interact with data stored on the Hutch servers, although I could use it on my Mac, I think.

Its' data viewer looks nice (`View()`) - includes little graphs


I am trying it on the work laptop (Dec 2025). Some notes from the walkthrough it gives on migrating from Rstudio:

- "Positron doesn't have an exact equivalent to RStudio Projects, but the concept in Positron that is most analogous to an RStudio Project is a workspace. You can read more in the VS Code documentation about what exactly a workspace is, but in general think of a workspace as about the same thing as a folder, which is about the same thing as an RStudio Project."  More info on [how to think about Rprojects in Positron](https://positron.posit.co/migrate-rstudio-rproj.html)
- "Air" is a formatting tool for R code


Databot - command-shift-P gives us the command pane, and there's an option to open databot. Select which LLM to use (choose Claude Sonnet), then we type requests in plain english. It generates code but doesn't run it without our permission.  You would probably want to copy the useful code chunks into your own qmd/Rmd document.  I have installed the extension but I haven't set it up to connect to a particular LLM yet (because I don't have accounts with the LLMs). It can make plots and tables and plain-english summaries of what the tables show.

After it runs it offers options, like to put its code/findings into markdown files

Databot demo videos from Ted Laderas (Hutch DASL): Using Databot on the NHANES dataset: 
- [part 1](https://www.youtube.com/watch?v=qs2GozYUUOk)
- [part 2](https://www.youtube.com/watch?v=lT2J71Jg_ug)
- [part 3](https://www.youtube.com/watch?v=j0KdDMIgcLY)



## Things I've learned

### Pipes:

I was using `magrittr`'s pipe in my code:   `%>%`.   

I recently switch to the native R pipe: `|>`. 

Can insert a pipe in Rstudio using Command-shift-m on the mac. There is a setting in R studio to tell it which of those two styles of pipe you want to use.

### R does rounding weirdly!

See [here](https://psiaims.github.io/CAMIS/R/rounding.html)

"The `round()` function in Base R will round to the nearest whole number and ‘rounding to the even number’ when equidistant, meaning that exactly 12.5 rounds to the integer 12. Note that the janitor package in R contains a function `round_half_up()` that rounds away from zero. in this case it rounds to the nearest whole number and ‘away from zero’ or ‘rounding up’ when equidistant, meaning that exactly 12.5 rounds to the integer 13.""

### Rstudio tricks

"Reindenting your code only shifts things around horizontally. If you want more powerful code reformatting, try using “Code > Reformat Code” (or use ⌘⇧A on macOS or ctrl + shift + A on Windows). It’s a more aggressive form of reformatting that will add extra line breaks and other things to make the code more readable."

In settings, there's an option to turn on rainbow parentheses to help see pairings.

[Snippets](https://rfortherestofus.com/2024/04/snippets-rstudio):

- e.g. the `fun` snippet:  within a code chunk, type `fun` and press tab - the skeleton of a function appears
- press tab again and you move within the snippet to the next piece you might fill in. Shift-tab also does something (not sure exactly what)
- `Tools menu - Edit code snippets` shows what snippets are available
- Markdown snippets are also useful (e.g. place an image). Here we need to do `shift-tab` to activate.  E.g. `r-shift-tab` inserts an R code chunk, if we do it from within a markdown area (i.e. outside an existing R code chunk)

#### Rstudio keyboard shortcuts

highlight a function, press function-F1, and it brings up the help page

shift-command-M types a pipe (there's a setting for whether you want that to be `%>%` or `|>`)

tab adds an indent to one or more selected lines of code, shift+tab removes one indent

(see tips [here](https://rfortherestofus.com/2023/11/rstudio-hotkeys))

#### Rstudio Snippets  

Example: type `fun` and press the `tab` key, and R provides the skeleton of a new function

To see all snippets:  Tools - Edit Code Snippets


### reprex

To create reproducible code + output we can use `reprex`

First, we write the code the demonstrates the problem

Load the library:
```
library(reprex)
```

Copy the code you want to the clipboard, and enter `reprex()` in the R Console. In RStudio, you’ll see a preview of your rendered reprex, but it is also now on the clipboard ready to paste.


### Debugging

Three options:    
- `browser()` (place inside a function, temporarily)    
- `debug(myFunction)` plus `undebug(myFunction)` (the Rstudio console window has a 'stop' button to exit debugging)   
- `debugonce()`    
See (`explore_debugging_functions.R`)[Rscripts/explore_debugging_functions.R] for details.

### Insert an image into an Rmd document

(or any md document, I think)

A simple image:
```
![Caption for the picture.](/path/to/image.png)
```

Modify size (specify pixels, probably)
```
![Caption for the picture.](/path/to/image.png){#id .class width=30 height=20px}
```

Modify size by simple scaling
```
![Caption for the picture.](/path/to/image.png){#id .class width=50% height=50%}
```

### Miscellaneous 



To extract some list elements using their indices in a pipe, use `magrittr::extract()`. This is the equivalent of using `[`.   Example:

```
my_list |> magrittr::extract(3:5) 
```

The related function `magrittr::extract2()` is the equivalent of `[[`, so we would use it to extract a single element from a list. Example:

```
my_list |> magrittr::extract2(1) 
```

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

`unz()` function lets you read files within a zip file without even unpacking it. See also `gzfile()` and `bzfile()`

## Useful packages

On a Mac, R packages go here - `/Library/Frameworks/R.framework/Versions` - in subdirectories by version. After installing new R, can delete old packages to save disk space


### combining plots

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

for genomes with karyotypes: [karyoploteR](https://bioconductor.org/packages/release/bioc/vignettes/karyoploteR/inst/doc/karyoploteR.html)

rtracklayer can make plots by interacting with a UCSC browser

igvR can interact with IGV

[tidyGenomeBrowser](https://github.com/MalteThodberg/tidyGenomeBrowser)

[GenomicPlot](https://bioconductor.org/packages/release/bioc/vignettes/GenomicPlot/inst/doc/GenomicPlot_vignettes.html) is more for making metaplots combining data over multiple features

Explored a few options in April 2024 for the SATAY data - see ~/FH_fast_storage/forOtherPeople/forGrantKing/SATAY/janet_Rscripts/ files browser_style_plots_failed_attempts.Rmd and browser_style_plots.Rmd


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


## Package installation note

There are some packages that won't install within an rhino/gizmo Rstudio-server session, because we need to load additional modules first, e.g. flextable and gdtools require the cairo module.

Here's what I'm trying, from a gizmo command line session
```
module purge
module load cairo/1.18.0-GCCcore-13.3.0
module load fhR/4.4.0-foss-2023b

R

install.package("gdtools")
install.package("flextable")
q(save="no")

module purge
```
Seems like that worked, and I can load those libraries within Rstudio-server without needing to load the cairo module.


## Knit/render an Rmd doc from the linux command line

Easiest way: on gizmo/rhino:

- `useful_functions/render_Rmd.perl` is a generic perl script that will take the name of one or more Rmd scripts as input, and make and run an sbatch script to render the report for each (in parallel, i.e. scripts are not linked).

- `useful_functions/render_Rmd_series.perl` is a generic perl script that will take the name of several Rmd scripts as input, and make and run an sbatch script to render the report for each, one after the other, stopping if any of them fail (scripts are in a linked series).

(there are also links to those scripts in `~/bin`)

Demo script and output, in `useful_functions/knit_using_shellScript`.

By default, if I render that way, plot images are saved as svg files, and those can be too big to sync to github.  Rendering using Rstudio saves png files, which are smaller. Dan Tenenbaum told me how to fix that, by adding this code in my Rmd doc:
```{r setup, include=FALSE}
ragg_png = function(..., res = 192) {
  ragg::agg_png(..., res = res, units = "in")
}
knitr::opts_chunk$set(dev = "ragg_png", fig.ext = "png")
```

I put those code lines in a (hidden) file called `~/.ragg_png_functions_from_dan.R` (also linked in an unhidden way, at `~/FH_fast_storage/git_more_repos/Rtest_and_Rnotes/useful_functions/knit_using_shellScript/ragg_png_functions_from_dan.R`).

I will source that from `.Rprofile` (global on rhino, and local if there is one), and then I shouldn't need to include it in my code.  Not putting in my Mac global for now - don't think it's relevant.

```
source("/home/jayoung/.ragg_png_functions_from_dan.R")
```

See also advice on ["Two Hidden Ways to Set Global Chunk Options for knitr"](https://yihui.org/en/2023/10/opts-chunk/).

### other notes about render from command line

Alternatives:

On my mac:
```
cd /Volumes/malik_h/user/jayoung/git_more_repos/Rtest_and_Rnotes

## I end up with TWO output files: .html AND .md - not sure why
## outputs go in the same dir as the input script

Rscript -e 'rmarkdown::render("Rscripts/miscellaneous_testCode.Rmd", output_format="github_document", clean=TRUE)'
```

On gizmo/rhino, from the command line without sbatch:
```
cd ~/FH_fast_storage/git_more_repos/Rtest_and_Rnotes

module load fhR/4.4.1-foss-2023b-R-4.4.1
module load Pandoc/2.13

## again we get TWO output files
Rscript -e 'rmarkdown::render("Rscripts/miscellaneous_testCode.Rmd", output_format="github_document", clean=TRUE)'

## note: running R this way, I do not have /home/jayoung/R/x86_64-pc-linux-gnu-library/4.4 in .libPaths(), whereas when I run it from Hutch RStudio-server I do.
## I don't know where that comes from - I suspect something to do with Rstudio server setup 
Rscript -e '.libPaths()'

module purge
```


On gizmo/rhino, using sbatch and a shell script, and run it like this: `cd ~/FH_fast_storage/git_more_repos/Rtest_and_Rnotes ; sbatch example_knit_batch_script.sh` - the script looks like this:
```
#!/bin/bash
source /app/lmod/lmod/init/profile
module load fhR/4.4.1-foss-2023b-R-4.4.1
module load Pandoc/2.13
Rscript -e 'rmarkdown::render("Rscripts/miscellaneous_testCode.Rmd", output_format="github_document", clean=TRUE)'
module purge
```



On gizmo/rhino, it can also be done with `sbatch --wrap` but it is annoying due to all the quote and escapes:
```
cd ~/FH_fast_storage/git_more_repos/Rtest_and_Rnotes
sbatch --cpus-per-task=1 --wrap="/bin/bash -c \"source /app/lmod/lmod/init/profile ; module load fhR/4.4.1-foss-2023b-R-4.4.1 ; module load Pandoc/2.13 ; Rscript -e 'rmarkdown::render(\\\"Rscripts/miscellaneous_testCode.Rmd\\\", output_format=\\\"github_document\\\", clean=TRUE)' ; module purge\""
```

For Tamanash's DMS data, I wrote an R wrapper script called [`DMS_remakeReport_standaloneScript.R`](https://github.com/jayoung/DMSanalysis/blob/main/Rscripts/DMS_remakeReport_standaloneScript.R) that generated reports for multiple data sets by (1) loading the Rdata file (2) rendering the Rmd script (which runs with that loaded data in the environment). I further wrote a perl script called [`remakeReport_wrapper.pl`](https://github.com/jayoung/DMSanalysis/blob/main/bin/remakeReport_wrapper.pl) that calls `DMS_remakeReport_standaloneScript.R` on various Rdata files and makes a shell script to run each in a separate sbatch job.
