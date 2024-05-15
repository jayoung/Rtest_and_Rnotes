## I think I was testing this out when I was working on Tamanash's first round of DMS data. I had a nice Rmd script that would generate a report on data passed into it, and I wanted to run the report in an automated way.  I think I deleted the paired testReport.Rmd file.

myDat <- cars

## myDat is in this script, but stays in the environment when I call the markdown this way
rmarkdown::render("Rscripts/testReport.Rmd", output_dir=".")
