#!/bin/bash


echo 'Running script script_using_mclapply.Rmd' >> zzz_Rmd_series.Rrender.log.txt
    Rscript -e 'rmarkdown::render("script_using_mclapply.Rmd", output_format="github_document", clean=TRUE)' > script_using_mclapply.Rrender.Rout.txt 2> script_using_mclapply.Rrender.Rerr.txt

STATUS=$?
if [ $STATUS -eq 0 ]
then
  echo '  SUCCEEDED' >> zzz_Rmd_series.Rrender.log.txt
else
  echo '  FAILED' >> zzz_Rmd_series.Rrender.log.txt
  exit 1
fi
rm script_using_mclapply.html 2>&1 >> zzz_Rmd_series.Rrender.log.txt

echo 'All scripts ran OK' >> zzz_Rmd_series.Rrender.log.txt

