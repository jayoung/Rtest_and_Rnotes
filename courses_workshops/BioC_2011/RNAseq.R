library(Biobase)
library(DEXSeq)
library(pasilla)

## DEXSeq exons
## DESeq genes

##### cuffdiff is a tool to find differential isoforms but it finds many positives on replicate samples

openVignette("DEXSeq")
openVignette("pasilla")


data("pasillaExons", package="pasilla")

# pasilla package also has its own vignette that describes data cleaning and construction of exon counts etc

pData(pasillaExons)
head(fData(pasillaExons))  #contains info on the counting bins. in this case exons

pasillaExons <- estimateSizeFactors(pasillaExons)   #takes geometric means of all samples to make a fake reference, and then gets all the ratios per exon and then takes median or something of those ratios as the sizeFactors

showMethods("estimateSizeFactors", inc=TRUE)

sizeFactors(pasillaExons)

## data stays as it is but these size factors 


pasillaExons <- estimateDispersions(pasillaExons)

head(modelFrameForGene(pasillaExons, "FBgn0010909"))

testGeneForDEU(pasillaExons, "FBgn0010909")
 ### DEU differential exon usage
 
pasillaExons <- testForDEU(pasillaExons)


res1 <- DEUresultTable(pasillaExons)

table(res1$padjust < 0.1)

### this is a nice plot to look at for many applications
hist(res1$pvalue,breaks=500)

# padjust is BH corrected (?)
hist(res1$padjust,breaks=50, xlim=c(0,0.5), ylim=c(0,20) )

design(pasillaExons)

formuladispersion <- count ~ sample + (exon + type) * condition
pasillaExons <- estimateDispersions(pasillaExons, formula = formuladispersion)


formula0 <- count ~ sample + type * exon + condition  ####### null model

formula1 <- count ~ sample + type * exon + condition * I(exon == exonID)  ##### now there's an effect by exon

pasillaExons <- testForDEU(pasillaExons, formula0 = formula0, formula1 = formula1)

res2 <- DEUresultTable(pasillaExons)

table(res2$padjust < 0.1)

plotDEXSeq(pasillaExons, "FBgn0010909", cex.axis = 1.2, cex = 1.3, lwd = 2, legend = TRUE)

## html report generation is very nice.
