##############################################
# Example 1 - find the nearest TSS for the peaks
##############################################
library(ChIPpeakAnno)
data(myPeakList)
data(TSS.human.GRCh37)
annotatedPeak = annotatePeakInBatch (myPeakList[1:6,], AnnotationData = TSS.human.GRCh37)
#The annotated peaks can be saved as an Excel file for biologists to view easily.
write.table(as.data.frame(annotatedPeak), file="annotatedPeakList.xls", sep="\t", row.names=FALSE) 

##############################################
# Example 2 - Plot the distribution of the peaks relative to the TSS 
# Gives a birds-eye view of the peak distribution relative to the genomic features of interest.
##############################################
data(annotatedPeak)
### STAT1 data 
y = annotatedPeak$distancetoFeature[!is.na(annotatedPeak$distancetoFeature) & annotatedPeak$fromOverlappingOrNearest == "NearestStart"] 

hist(y, xlab="Distance To Nearest TSS", main="", breaks=1000, xlim=c(min(y)-100, max(y)+100)) 

temp = as.data.frame(annotatedPeak)
plot(density(y))

y = annotatedPeak$distancetoFeature[!is.na(annotatedPeak$distancetoFeature) & annotatedPeak$fromOverlappingOrNearest == "NearestStart" & abs(annotatedPeak$distancetoFeature) <10000]

pie(table(temp[as.character(temp$fromOverlappingOrNearest) == "Overlapping" | (as.character(temp$fromOverlappingOrNearest) == "NearestStart" & !temp$peak %in%  temp[as.character(temp$fromOverlappingOrNearest) == "Overlapping", ]$peak) ,]$insideFeature))

##############################################
# Example 3 - Obtain annotation on-line using getAnnotation
##############################################
mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") 
#Annotation = getAnnotation(mart,  featureType="TSS")
Annotation = getAnnotation(mart,  featureType="miRNA") ##### RangedData
as.data.frame(Annotation)[1:10,]

##############################################
# Example 4 - Label the peaks from your experiment with a list of peaks in the literature
# including both nearest and overlapping sites
##############################################
myexp =  RangedData(IRanges(start=c(1543200,1557200,1563000,1569800, 167889600,100,1000), end=c(1555199,1560599,1565199,1573799, 167893599,200,1200), names=c("p1","p2","p3","p4","p5","p6", "p7")), strand=as.integer(1),space=c(6,6,6,6,5,4,4))

literature = RangedData(IRanges(start=c(1549800,1554400,1565000,1569400,167888600,120,800), end=c(1550599,1560799,1565399,1571199,167888999,140,1400), names=c("f1","f2","f3","f4","f5","f6","f7")), strand=c(1,1,1,1,1,-1,-1), space=c(6,6,6,6,5,4,4))

annotatedPeak1= annotatePeakInBatch(myexp, AnnotationData = literature, output="both", maxgap=1000, multiple=TRUE)

pie(table(as.data.frame(annotatedPeak1)$insideFeature))

as.data.frame(annotatedPeak1)

##### different parameter setting
annotatedPeak1= annotatePeakInBatch(myexp, AnnotationData = literature, output="overlapping", maxgap=1000, multiple=TRUE)

as.data.frame(annotatedPeak1)

annotatedPeak1= annotatePeakInBatch(myexp, AnnotationData = literature, output="nearestStart", PeakLocForDistance ="middle", FeatureLocForDistance="middle")

as.data.frame(annotatedPeak1)

##############################################
# Example 5 - BED2RangedData and GFF2RangedData
##############################################
test.bed = data.frame(cbind(chrom = c("4", "6"), chromStart=c("100", "1000"),chromEnd=c("200", "1100"), name=c("peak1", "peak2")))

test.rangedData = BED2RangedData(test.bed)

as.data.frame(annotatePeakInBatch(test.rangedData, AnnotationData = literature))

test.GFF = data.frame(cbind(seqname  = c("chr4", "chr4"), source=rep("Macs", 2), feature=rep("peak", 2), start=c("100", "1000"), end=c("200", "1100"), score=c(60, 26), strand=c(1, 1), frame=c(".", 2), group=c("peak1", "peak2")))

test.rangedData = GFF2RangedData(test.GFF)

as.data.frame(annotatePeakInBatch(test.rangedData, AnnotationData = literature))

##############################################
# Example 6 - Determine the significance of the overlapping and 
# visualize the overlap as a Venn diagram among different datasets
##############################################
data(Peaks.Ste12.Replicate1)
data(Peaks.Ste12.Replicate2)
data(Peaks.Ste12.Replicate3)


makeVennDiagram(RangedDataList(Peaks.Ste12.Replicate1, Peaks.Ste12.Replicate2, Peaks.Ste12.Replicate3), NameOfPeaks = c("Replicate1","Replicate2","Replicate3"), maxgap = 0, totalTest = 1580)

#makeVennDiagram(RangedDataList(myexp, literature), NameOfPeaks = c("myexp","literature"), maxgap = 0, totalTest = 1580)

#MergedPeaks = findOverlappingPeaks(myexp, literature, maxgap = 0, multiple = F, NameOfPeaks1 = "R1", NameOfPeaks2 = "R2")$MergedPeaks


#Combine the Overlapping Peaks Across Replicates
MergedPeaks = findOverlappingPeaks(findOverlappingPeaks(Peaks.Ste12.Replicate1, Peaks.Ste12.Replicate2, maxgap = 0, multiple = F, NameOfPeaks1 = "R1", NameOfPeaks2 = "R2")$MergedPeaks, Peaks.Ste12.Replicate3, maxgap = 0, multiple = F, NameOfPeaks1 = "R1R2",  NameOfPeaks2 = "R3")$MergedPeak 
as.data.frame(MergedPeaks)

##############################################
# Example 7 - Obtain the sequences around the binding sites for PCR amplification or motif discovery
##############################################
peaks = RangedData(IRanges(start = c(100, 500), end = c(300, 
 600), names = c("peak1", "peak2")), space = c("NC_008253", 
 "NC_010468"))
 
library(BSgenome.Ecoli.NCBI.20080805)

peaksWithSequences = getAllPeakSequence(peaks, upstream = 100, downstream = 100, genome = Ecoli)

write2FASTA(peaksWithSequences, file="test.fa", width=50)
#available.genomes() 


##############################################
# Example 8 - Obtain enriched GO terms near the peaks 
##############################################
data(annotatedPeak)
library(org.Hs.eg.db) 

enrichedGO <- getEnrichedGO (annotatedPeak[1:6,], orgAnn="org.Hs.eg.db", maxP=0.1, multiAdj =TRUE, minGOterm=1,  multiAdjMethod="BH")

##############################################
# Example 9 - Find the peaks with bi-directional promoters with summary statistics (>= version1.9.6)
##############################################

temp = peaksNearBDP(annotatedPeak[1:10,], AnnotationData=TSS.human.GRCh37, MaxDistance=5000, PeakLocForDistance="middle", FeatureLocForDistance = "TSS")
c(temp$percentPeaksWithBDP, temp$n.peaksWithBDP, temp$n.peaks)

temp$peaksWithBDP

##############################################
# Example 10 - Summarize the occurrence of motifs in peaks (>= version1.9.6)
##############################################

filePath = system.file("extdata", "examplePattern.fa", package="ChIPpeakAnno")

summarizePatternInPeaks(patternFilePath=filePath,skip=0L, format="fasta", BSgenomeName=Ecoli, peaks=peaks, outfile="testExample10.xls", append=FALSE)

file.remove("testExample10.xls")


