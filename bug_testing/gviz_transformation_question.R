# April 15, 2024
# https://github.com/ivanek/Gviz/issues/90

library(Gviz)

## generate test data
# middle position has two overlapping regions, but all positions sum to the same count. Looks fine on a linear scale but not on a log scale
test_sites <- GRanges(seqnames="chr1",
                      ranges=IRanges(start=c(10,20,20,30), width=1),
                      counts=c(100,50,50,100)) 
test_region <- GRanges(seqnames="chr1",
                       ranges=IRanges(start=1, end=40)) 

## make axis track
myTrack1 <- GenomeAxisTrack(test_region)

## make raw counts track
myTrack2a <- DataTrack(
    test_sites, 
    name="counts",
    region=test_region,
    type = c("h"),
    window=width(test_region), 
    aggregation="sum",
    ylim=c(0,150)
)

## use transformation option to get it on a log10 scale. Now the middle position looks like it has a higher score.
myTrack2b <- DataTrack(
    test_sites, 
    name="log10 counts v1 ",
    region=test_region,
    type = c("h"),
    window=width(test_region), 
    aggregation="sum",
    transformation=function(x) {log10(x)},
    ylim=c(0,3.5)
)

## SOLUTION - do the log10 scale transformation during window aggregation instead
myTrack2c <- DataTrack(
    test_sites, 
    name="log10 counts v2",
    region=test_region,
    type = c("h"),
    window=width(test_region), 
    aggregation=function(x) {log10(sum(x))},
    ylim=c(0,3.5)
)

myTrackList <- list(myTrack2c,myTrack2b,myTrack2a,myTrack1)

plotTracks(myTrackList)

png(filename="bug_testing/gviz_transformation_question.png")
plotTracks(myTrackList)
dev.off()

sessionInfo()

