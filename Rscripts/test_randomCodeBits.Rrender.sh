#!/bin/bash
source /app/lmod/lmod/init/profile
module purge
module load fhR/4.4.1-foss-2023b-R-4.4.1
module load Pandoc/2.13
Rscript -e 'rmarkdown::render("test_randomCodeBits.Rmd", output_format="github_document", clean=TRUE)' > test_randomCodeBits.Rrender.Rout.txt 2> test_randomCodeBits.Rrender.Rerr.txt
rm test_randomCodeBits.html
module purge
