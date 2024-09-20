# my_tbl a tibble containing the data. 
# groupA and groupB - colnames of columns containing the two groups of data we want to compare.
# pseudoCountForLog - for ratios, you might want to add a small number to numerator and denominator to avoid problems of dividing by zero.  Makes ratios from small counts more conservative.
# paired -  are observations in the two groups paired?
doTtests <- function(my_tbl, groupA, groupB, 
                     pseudoCountForLog=0,
                     paired=FALSE) {
  # check all sample names supplied are in the data
  if( sum(!groupA %in% colnames(my_tbl))>0 ) {
    stop("\n\nERROR - you specified samples in groupA that aren't in the data\n\n")
  }
  if( sum(!groupB %in% colnames(my_tbl))>0 ) {
    stop("\n\nERROR - you specified samples in groupB that aren't in the data\n\n")
  }
  my_tbl_A <- my_tbl %>% 
    select(all_of(groupA))
  my_tbl_B <- my_tbl %>% 
    select(all_of(groupB))
  # do t.test for each row
  results <- lapply(1:dim(my_tbl)[1], function(i) {
    res <- t.test(x=unlist(my_tbl_A[i,]), 
                  y=unlist(my_tbl_B[i,]),
                  paired=paired)
    return(res)
  })
  # get p-values and rowMeans
  outputToAdd <- tibble(
    meanA=rowMeans(my_tbl_A),
    meanB=rowMeans(my_tbl_B),
    log2fc_AvsB=log2( (meanB+pseudoCountForLog) / 
                        (meanA+pseudoCountForLog) ),
    tTestPval= sapply(results, function(x) {x$p.value})
  )
  output <- bind_cols(my_tbl, outputToAdd)
  # do multiple testing correction
  output$tTestPvalAdj <- p.adjust(output$tTestPval, 
                                  method = "BH")
  return(output)
}

