#!/bin/bash
source /app/lmod/lmod/init/profile
module purge
module load fhR/4.4.1-foss-2023b-R-4.4.1
module load Pandoc/2.13


echo 'Running script script_to_knit_without_raggChunk.Rmd' >> zzz_Rmd_series.Rrender.log.txt
Rscript -e 'rmarkdown::render("script_to_knit_without_raggChunk.Rmd", output_format="github_document", clean=TRUE)' > script_to_knit_without_raggChunk.Rrender.Rout.txt 2> script_to_knit_without_raggChunk.Rrender.Rerr.txt

STATUS=$?
if [ $STATUS -eq 0 ]
then
  echo '  SUCCEEDED' >> zzz_Rmd_series.Rrender.log.txt
else
  echo '  FAILED' >> zzz_Rmd_series.Rrender.log.txt
  exit 1
fi
rm script_to_knit_without_raggChunk.html 2>&1 >> zzz_Rmd_series.Rrender.log.txt

module purge

echo 'All scripts ran OK' >> zzz_Rmd_series.Rrender.log.txt

