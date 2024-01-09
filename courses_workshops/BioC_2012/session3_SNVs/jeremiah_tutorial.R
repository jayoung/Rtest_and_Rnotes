################################################################################
#Variant calling in R, some background and an Introduction to SNVsOmuC and gmapR
################################################################################

#Introduction

##This doc outlines the basic usages of the SNVsOmuC and gmapR packages as well
##as some basic usage of  other BioC packages for the manipulation of short-read alignments
##and variant data.
##I cover the general workflow for loading data, calling single sample variants and
##tumor-specific somatic mutations or other sample-specific variant types (eg RNA editing).

###manipulating bam files
source("~/.Rprofile")
download.file("http://bioconductor.org/help/course-materials/2012/BioC2012/jeremiah2.R", destfile="jeremiah.R")
library(SNVsOmuC)
library(gmapR)
library(Rsamtools)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)

bamfile = file.path(system.file(package = "SNVsOmuC", mustWork=TRUE),
    "extdata/merged/SRC111111.TESTPREPENDSTR.merged.uniques.bam")

############################################################
#A quick overview of some tools for the summary of BAM files
############################################################

##Reading bam files using Rsamtools
nrows=100
bf <- BamFile(file = bamfile, yieldSize=nrows)
bam <- scanBam(bf)
bam

## checking the quality type in the bam file
nrows = 100
bf <- BamFile(file = bamfile, yieldSize=nrows)
?ScanBamParam
param <- ScanBamParam(what="qual")
bam <- scanBam(bf, param=param)
bam
quals <- as.character(bam[[1]]$qual)
length <- length(quals)
max_vec <- array()
min_vec <- array()
for(i in seq_along(length)){
    max_vec[i]<-max(as.integer(charToRaw(quals[i])))
      min_vec[i]<-min(as.integer(charToRaw(quals[i])))
  }
min <- min(min_vec)
max <- max(max_vec)
max
min

qual <- detectQuality(bamfile)
qual

##getting the total coverage
?readBamGappedAlignments
chr_ga <-  readBamGappedAlignments(bamfile)
pos <- chr_ga[strand(chr_ga)=="+"]
neg <- chr_ga[strand(chr_ga)=="-"]
pos_grl <- grglist(pos, drop.D.ranges=TRUE)
neg_grl <- grglist(neg, drop.D.ranges=TRUE)
pos_cov <- coverage(pos_grl)
neg_cov <- coverage(neg_grl)
cov <- pos_cov+neg_cov
plot(x=1, pch = NA, ylim = c(-120, 200), xlim = c(0,1200), ylab = "stranded covergage", xlab = "Position")
points(as.numeric(pos_cov$chr12[6646350:6647650]), type = 'l', col = "red")
points(-as.numeric(neg_cov$chr12[6646350:6647650]), type = 'l', col = "blue")
points(as.numeric(cov$chr12[6646350:6647650]), type = 'l', col = "grey")
abline(h=0)

##creating pileup data

?applyPileups
?PileupParam

which <- GRanges(seqnames = "chr12", IRanges(start =1, width = 6647650))
source("/home/jeremiah/pileup_to_gr.R")
bu_gr <- pileupAsGRanges(bamfile, regions = which)

bases <- t(as.matrix(values(bu_gr)[1:4]))
bases[,800:810]
barplot(bases[,1:1272], col=c("red", "blue", "green", "black"), border=NA)
points(as.numeric(cov$chr12[6646350:6647650]), type = 'l', col = "grey")

###############################################
#Intro to the new packages gmapR and SNVsOmuC
###############################################

##GSNAP alignments in R
dnaStringSet <-Hsapiens$chr4
makeTestFasta <- function(seqStarts, pairNum) {
    seqEnds <- seqStarts + 74
      nucleotides <- sapply(seq_along(seqStarts),
                                                    function(i) as.character(subseq(dnaStringSet, seqStarts[i], seqEnds[i]))
                                                    )
      fastaHeaders <- paste0(">", "fakeRead", seq_along(nucleotides), "#0/", pairNum)
      fastaLines <- paste(fastaHeaders, nucleotides, sep="\n")
      tmpFasta <- paste0(tempfile(), ".", pairNum, ".fa")
      writeLines(fastaLines, con=tmpFasta)
      tmpFasta
  }
seqStarts <- seq(1, 100000, by = 10000) + length(dnaStringSet) * 2 / 3
tmpFasta1 <- makeTestFasta(seqStarts, '1')
tmpFasta2 <- makeTestFasta(seqStarts + 125, '2')


### this part actually runs GSNAP - sets up parameters, etc.  All of GSNAP is installed when you install the package (no external installation needed)

library(GmapGenome.Hsapiens.UCSC.hg19)

genome <- GmapGenome.Hsapiens.UCSC.hg19

?GsnapParam

gsnapParam <- GsnapParam(genome = genome,
                         unique_only = FALSE,
                         max_mismatches = NULL,
                         suboptimal_levels = 0L, mode = "standard",
                         npaths = 10L,
                         novelsplicing = FALSE, splicing = NULL,
                         nthreads = 1L,
                         batch = "2")

gsnapOutput <- gsnap(input_a=tmpFasta1,
                                          input_b=tmpFasta2,
                                          params=gsnapParam)


path <- path(gsnapOutput)
##bam <- paste(path, dir(path)[8], sep = '/')
chr_ga <-  readBamGappedAlignments(bam)
grl <- grglist(chr_ga, drop.D.ranges=TRUE)
cov <- coverage(grl)

##Calling single-sample vairants

##The gmapR package has a function which converts the tally format from gmap
##(similar to samtools pileup) into a variant GRanges object. This object contains
##all of the information needed to call variants. The convience function below combines
##the gmapR tally function with the variant calling functionality here to simultaneously
##tally the BAM file and call the variants.

##genome <- library(GmapGenome.Hsapiens.UCSC.hg19)

?BamTallyParam

### take note of what qualtiy-encoding is used for high_qual_cutoff

param <- BamTallyParam(cycle_breaks = c(0L, 15L, 60L, 75L),
                                              high_quality_cutoff = as.integer(56),
                                              minimum_mapq = as.integer(13),
                                              concordant_only = FALSE, unique_only = FALSE,
                                              primary_only = FALSE,
                                              min_depth = 0L, variant_strand = 1L,
                                              ignore_query_Ns = TRUE,
                                              indels = FALSE)
tallies_gr <-bam_tally(bamfile,
                                              genome,
                                              param)
tallies_gr

## ncycles is how many different positions in reads was it seen in 

##Calling variants on a subset of the genome
### e.g. to allow parallelizing over many machines by splitting up genome
which <- GRanges(seqnames = "chr12", IRanges(start =6647109, width = 6647453))
param <- BamTallyParam(which=which, cycle_breaks = c(0L, 15L, 60L, 75L),
                       high_quality_cutoff = as.integer(56),
                       minimum_mapq = as.integer(13),
                       concordant_only = FALSE, unique_only = FALSE,
                       primary_only = FALSE,
                       min_depth = 0L, variant_strand = 1L,
                       ignore_query_Ns = TRUE,
                       indels = FALSE)
tallies_gr <-bam_tally(bamfile,
                       genome,
                       param)
tallies_gr

####################################################################
#Provide a method for filtering the pileup data and calling variants
####################################################################

var <- variantFilter(tallies_gr, useQual=TRUE)
var$filtered_granges

##Convenience method for variant calling

##A convenience method is provided to call variants across the entire genome
##given a BAM file as input. This section assumes that you have installed gmapR
##and built a genome iit file.

##variantsGR <- callVariantsP(bam=bamfile, genome=genome(genome),
##                             genome_dir=path(directory(genome)),mapq=13, force_qual=TRUE)
##variantsGR$filtered_granges

##A real example...
### for now filter getting just chrs 20, 21, 22 and variants with at least 15 reads
variantsGR<- get(load("/home/jeremiah/SRR305173/results/SRR305173.filtered_variants_granges.RData"))
chr2X <- variantsGR[seqnames(variantsGR) %in% c("20", "21", "22")]
chr2X <- chr2X[values(chr2X)$count.total>15]

### plot to see allele frequency among reads (not the way I'd look at it)
plot(log10(values(chr2X)$count), log10(values(chr2X)$count.total), pch = 19, cex= 0.01)

### to filter for things that are at 50% or higher (note: not just a simple 0.5 threshold)
pval <- (pbinom(q=values(chr2X)$count, size=values(chr2X)$count.total,prob= .5))
good <- chr2X[pval>0.001]
plot(log10(values(good)$count), log10(values(good)$count.total), pch = 19, cex= 0.01)

hist(values(good)$count/values(good)$count.total,
          breaks=50, main="SRR305173 SNPs on Chrom 20-22")

##Two sample comparison methods
tumor <- variantsGR
tumor <- tumor[seqnames(tumor) %in% c("20", "21", "22")]
normal <- get(load("/home/jeremiah/SRR305174/results/SRR305174.filtered_variants_granges.RData"))
normal <- normal[seqnames(normal) %in% c("20", "21", "22")]
raw <-get(load("/home/jeremiah/raw_chr_subset_gr.RData"))
cov<- get(load("/home/jeremiah/SRR305174/SRR305174.coverage.RData"))


##Calling sample-specific mutations

##next we want to compare the tumor and normal samples to find the tumor-specific
##mutations...

##however, before doing so we need to know that the two samples are actually from
##the same patient

##Checking that we have data from the same person

concord <- variantConcordance(tumor, normal)
concord

#reported 97% concordance: fairly usual

##now that we know we have tumor and normal data from the same sample we need one
##more bit of data before we get the tumor-specific variants

##Actually calling the sample-specific variants

TS<- tumorNormalCompare(tumor_gr=tumor,
                        normal_gr=normal, normal_raw=raw,
                        normal_cov=cov)


##Writing a vcf file

##Next we want to write out our tumor-specific varints out to a vcf file

cgpGr2vcf(TS, file = "~/temp_ts.vcf", sample_id="Test_sam", project = "SNVsOmuC_Vignette")

### currently just a single sample at a time - need to do once per sample of interest

##annotaing variants and loading VCF files
library(VariantAnnotation)
fl <- "/home/jeremiah/temp_ts.vcf.gz"
vcf <- readVcf(fl, "hg19")
vcf
rowData(vcf)

#######################################################################
## We can also use VariantAnnotation to identify potentially functional
## using the GRanges objects generated by gmapR and SNVsOmuC directly
#######################################################################

seqlevels(TS) <- paste("chr", seqlevels(TS), sep = "")
genome(TS) <- "hg19"
seqlevels(TS)<-seqlevels(TS)[1:24]    # to avoid warnings
strand(TS)<-"*"   ## otherwise we would only get + strand genes
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

coding <- predictCoding(TS, varAllele=DNAStringSet(values(TS)$alt), txdb, seqSource=Hsapiens)
coding


##Making some plots

##Finally we want make some plots to check out our data.

##The first plot shows us the mutation transition/transversion rate matrix

plotTitv(variantsGR, main = "Single-sample SNVs")  ### TiTv transition/transversion ratios
plotTitv(TS, main = "tumor specific mutations")


##and finally we want to plot our variants on the genome

seqlevels(tumor) <- paste("chr", seqlevels(tumor), sep = "")
genome(tumor) <- "hg19"
seqlevels(tumor)<-seqlevels(tumor)[1:24] 

plotTumor(tumor,TS)

