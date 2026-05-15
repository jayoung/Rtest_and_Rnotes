#!/bin/bash


echo 'Running script script1.Rmd' >> zzz_Rmd_series.Rrender.log.txt
    Rscript -e 'rmarkdown::render("script1.Rmd", output_format="github_document", clean=TRUE)' > script1.Rrender.Rout.txt 2> script1.Rrender.Rerr.txt

STATUS=$?
if [ $STATUS -eq 0 ]
then
  echo '  SUCCEEDED' >> zzz_Rmd_series.Rrender.log.txt
else
  echo '  FAILED' >> zzz_Rmd_series.Rrender.log.txt
  exit 1
fi
rm script1.html 2>&1 >> zzz_Rmd_series.Rrender.log.txt


echo 'Running script script2.Rmd' >> zzz_Rmd_series.Rrender.log.txt
    Rscript -e 'rmarkdown::render("script2.Rmd", output_format="github_document", clean=TRUE)' > script2.Rrender.Rout.txt 2> script2.Rrender.Rerr.txt

STATUS=$?
if [ $STATUS -eq 0 ]
then
  echo '  SUCCEEDED' >> zzz_Rmd_series.Rrender.log.txt
else
  echo '  FAILED' >> zzz_Rmd_series.Rrender.log.txt
  exit 1
fi
rm script2.html 2>&1 >> zzz_Rmd_series.Rrender.log.txt

echo 'All scripts ran OK' >> zzz_Rmd_series.Rrender.log.txt

