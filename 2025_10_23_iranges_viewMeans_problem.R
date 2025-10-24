### https://github.com/Bioconductor/IRanges/issues/63

## load IRanges
library(IRanges)

### can't create a Views object yet - works if I load GenomicRanges
my_views <- Views(c(1:3,NA,NA,6:10), 
                  start=1:6, width=5)



## load GenomicRanges
library(GenomicRanges)

my_views <- Views(c(1:3,NA,NA,6:10), 
                  start=1:6, width=5)

my_views_withoutNA <- Views(c(1:10), 
                            start=1:6, width=5)

## this works
viewSums(my_views, na.rm=TRUE)

## viewMeans works if na.rm=FALSE
viewMeans(my_views, na.rm=FALSE)

## but fails if na.rm=TRUE, even if the object doesn't have any NAs
viewMeans(my_views, na.rm=TRUE)
viewMeans(my_views_withoutNA, na.rm=TRUE)

## viewApply gives warnings, and all NA results if we don't coerce to integer
viewApply(my_views, function(x) { mean(x, na.rm=FALSE) })
viewApply(my_views, function(x) { mean(x, na.rm=TRUE) })
viewApply(my_views_withoutNA, function(x) { mean(x, na.rm=FALSE) })
viewApply(my_views_withoutNA, function(x) { mean(x, na.rm=TRUE) })

## viewApply works if we coerce to integer first
viewApply(my_views, function(x) { mean(as.integer(x), na.rm=FALSE) })
viewApply(my_views, function(x) { mean(as.integer(x), na.rm=TRUE) })

sessionInfo()

