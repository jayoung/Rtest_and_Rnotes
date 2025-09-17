#!/bin/bash
source /app/lmod/lmod/init/profile
module load fhR/4.4.1-foss-2023b-R-4.4.1
module load Pandoc/2.13
Rscript -e 'rmarkdown::render("Rscripts/miscellaneous_testCode.Rmd", output_format="github_document", clean=TRUE)'
module purge

