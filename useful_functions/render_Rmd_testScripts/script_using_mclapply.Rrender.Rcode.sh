#!/bin/bash


Rscript -e 'rmarkdown::render("script_using_mclapply.Rmd", output_format="github_document", clean=TRUE)' > script_using_mclapply.Rrender.Rout.txt 2> script_using_mclapply.Rrender.Rerr.txt

rm script_using_mclapply.html 2>&1 >> script_using_mclapply.Rrender.log.txt
