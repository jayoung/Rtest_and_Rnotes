#!/bin/bash
source /app/lmod/lmod/init/profile
module purge
module load fhR/4.4.1-foss-2023b-R-4.4.1
module load Pandoc/2.13
Rscript -e 'rmarkdown::render("script_to_knit_without_raggChunk.Rmd", output_format="github_document", clean=TRUE)' > script_to_knit_without_raggChunk.Rrender.Rout.txt 2> script_to_knit_without_raggChunk.Rrender.Rerr.txt
rm script_to_knit_without_raggChunk.html
module purge
