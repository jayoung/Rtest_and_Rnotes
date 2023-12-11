myDat <- cars

## myDat is in this script, but stays in the environment when I call the markdown this way
rmarkdown::render("Rscripts/testReport.Rmd", output_dir=".")
