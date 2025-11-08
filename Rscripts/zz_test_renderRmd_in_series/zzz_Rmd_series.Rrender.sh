#!/bin/bash
source /app/lmod/lmod/init/profile
module purge
module load fhR/4.4.1-foss-2023b-R-4.4.1
module load Pandoc/2.13


echo 'Running script test_randomCodeBits_part1.Rmd' >> zzz_Rmd_series.Rrender.log.txt
Rscript -e 'rmarkdown::render("test_randomCodeBits_part1.Rmd", output_format="github_document", clean=TRUE)' > test_randomCodeBits_part1.Rrender.Rout.txt 2> test_randomCodeBits_part1.Rrender.Rerr.txt

STATUS=$?
if [ $STATUS -eq 0 ]
then
  echo '  SUCCEEDED' >> zzz_Rmd_series.Rrender.log.txt
else
  echo '  FAILED' >> zzz_Rmd_series.Rrender.log.txt
  exit 1
fi
rm test_randomCodeBits_part1.html 2>&1 >> zzz_Rmd_series.Rrender.log.txt


echo 'Running script test_randomCodeBits_part2.Rmd' >> zzz_Rmd_series.Rrender.log.txt
Rscript -e 'rmarkdown::render("test_randomCodeBits_part2.Rmd", output_format="github_document", clean=TRUE)' > test_randomCodeBits_part2.Rrender.Rout.txt 2> test_randomCodeBits_part2.Rrender.Rerr.txt

STATUS=$?
if [ $STATUS -eq 0 ]
then
  echo '  SUCCEEDED' >> zzz_Rmd_series.Rrender.log.txt
else
  echo '  FAILED' >> zzz_Rmd_series.Rrender.log.txt
  exit 1
fi
rm test_randomCodeBits_part2.html 2>&1 >> zzz_Rmd_series.Rrender.log.txt

module purge

echo 'All scripts ran OK' >> zzz_Rmd_series.Rrender.log.txt

