# Knit/render an Rmd doc from the linux command line

Easiest way: on gizmo/rhino:

- `useful_functions/render_Rmd.pl` is a generic perl script that will take the name of one or more Rmd scripts as input, and make and run an sbatch script to render the report for each (in parallel, i.e. scripts are not linked).

- `useful_functions/render_Rmd_series.pl` is a generic perl script that will take the name of several Rmd scripts as input, and make and run an sbatch script to render the report for each, one after the other, stopping if any of them fail (scripts are in a linked series).

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
