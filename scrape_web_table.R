library(XML)
library(RCurl)
library(rlist)
theurl <- getURL("https://en.wikipedia.org/wiki/Brazil_national_football_team",.opts = list(ssl.verifypeer = FALSE) )
tables <- readHTMLTable(theurl)
tables <- list.clean(tables, fun = is.null, recursive = FALSE)
n.rows <- unlist(lapply(tables, function(t) dim(t)[1]))

x <- tables[[which.max(n.rows)]]


vep_summary_file <- "/Volumes/malik_h/user/jayoung/bat_reproduction/data/fastq/artibeus_jamaicensis/genotyping/refseq_reference/genotyping_4_2023_Jan27/zz_testVEP/testRegion_47snps.snpsANDindels.filt.VEP.vcf_summary.html" 
vep_summary <- readHTMLTable(vep_summary_file)


# the tables I want are a different style/method (drawTable - google.visualization.arrayToDataTable) from the two tables I get using readHTMLTable (class="stats_table")

