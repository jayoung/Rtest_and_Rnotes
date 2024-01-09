Several BioConductor packages are available for displaying genomic information graphically. The following commands illustrate some of the chromosome plotting utilities from the geneplotter package.


library(annotate)
library(geneplotter)
library(hgu95av2)
newChrom <- buildChromLocation("hgu95av2")
newChrom
cPlot(newChrom) # This displays all genes on the chromosomes of an organisms. Genes encoded by the antisense strands are represented by lines below the chromosomes.

data(sample.ExpressionSet)
myeset <- sample.ExpressionSet
cColor(featureNames(sample.ExpressionSet), "red", newChrom) # This highlights in the above plot a set of genes of interest in red color (e.g. expressed genes of an experiment).

cPlot(newChrom,c("1","2"), fg="yellow", scale="relative")
cColor(featureNames(myeset), "red", newChrom) # Plots just a specific set of chromosomes.

#### OK - not very pretty.  plots tick marks on chr axes.






#######


library(GenomeGraphs)
library(biomaRt)
data("exampleData", package = "GenomeGraphs")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
minbase <- 180292097
maxbase <- 180492096
genesplus <- makeGeneRegion(start = minbase,
                            end = maxbase, strand = "+", chromosome = "3",
                            biomart = mart)
genesmin <- makeGeneRegion(start = minbase,
                           end = maxbase, strand = "-", chromosome = "3",
                           biomart = mart)
seg <- makeSegmentation(segStart[[1]], segEnd[[1]],
                        segments[[1]], dp = DisplayPars(color = "black",
                                    lwd = 2, lty = "solid"))
cop <- makeGenericArray(intensity = cn,
                        probeStart = probestart, segmentation = seg,
                        dp = DisplayPars(size = 3, color = "seagreen",
                          type = "dot"))
cop <- makeGenericArray(intensity = cn,
                        probeStart = probestart, 
                        dp = DisplayPars(size = 3, color = "seagreen",
                          type = "dot"))
ideog <- makeIdeogram(chromosome = 3)
expres <- makeGenericArray(intensity = intensity,
                           probeStart = exonProbePos,
                           dp = DisplayPars(color = "darkred",
                             type = "point"))
genomeAxis <- makeGenomeAxis(add53 = TRUE,
                             add35 = TRUE)
gdPlot(list(a = ideog, b = expres, c = cop,
            d = genesplus, e = genomeAxis, f = genesmin),
       minBase = minbase, maxBase = maxbase,
       labelCex = 2)
       
       
       
       
ideog <- makeIdeogram(chromosome = 1)
ideog2 <- makeIdeogram(chromosome = 2)
ideog3 <- makeIdeogram(chromosome = 3)
ideog4 <- makeIdeogram(chromosome = 4)
gdPlot(list("1"= ideog, "2" = ideog2, "3" =ideog3, "4"=ideog4 ), minBase = minbase, maxBase = maxbase)