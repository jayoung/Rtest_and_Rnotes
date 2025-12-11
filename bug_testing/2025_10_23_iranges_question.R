# https://github.com/Bioconductor/IRanges/issues/62
# 
# hi there,
# 
# I've recently been loving switching back and forth between tibbles and GRanges objects with ease. But I just realized there's no equivalent method to coerce tibble (or data.frame) to IRanges. Not that it's very difficult to do with a couple more lines of code, but it'd be nice to have the shortcut.
# 
# Here's a reprex.
# 
# thanks!
# 
# Janet


library(GenomicRanges)
library(dplyr)
x <- tibble(start=1:3, end=4:6)
y <- x |> mutate(seqnames="chr1")
GRanges(y)


IRanges(x)

IRanges(as.data.frame(x))

sessionInfo()
