### R code from vignette source 'tutorial.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width=50)
library(RangesRNAseqTutorial)


###################################################
### code chunk number 2: rnaseq-readGappedAlignments
###################################################
bams <- dir(system.file("bam", package="leeBamViews"), 
            full = TRUE, pattern = "bam$")
reads_ga <- readGappedAlignments(bams[[1]])
head(reads_ga, 1)

head(width(reads_ga),1)

###################################################
### code chunk number 3: rnaseq-grglist
###################################################
reads_grl <- grglist(reads_ga)

table(elementLengths(reads_grl)) ### shows that all are length 1, i.e. all unspliced

###################################################
### code chunk number 4: expr-cds
###################################################
library(TxDb.Scerevisiae.UCSC.sacCer2.sgdGene)
sc2_txdb <- TxDb.Scerevisiae.UCSC.sacCer2.sgdGene
sc2_cds <- cds(sc2_txdb)   ### GRanges
sc2_cds_gene <- cdsBy(sc2_txdb, by = "gene") ### GRangesList

sc2_codingexons_tx <- cdsBy(sc2_txdb, by="tx")
### should still work even for human genes - should still give individual exons

table(sapply(sc2_cds_gene,length)) # some genes have >1 CDS

###################################################
### code chunk number 5: expr-cds-tx
###################################################
sc2_cds_tx <- cdsBy(sc2_txdb, by = "tx")

class(sc2_cds_tx)  ### GRangesList

###################################################
### code chunk number 6: gene-ids
###################################################
sc2_tx <- transcripts(sc2_txdb, 
                      columns = c("tx_id", "gene_id"))
tx_id_match <- match(names(sc2_cds_tx), 
                     values(sc2_tx)$tx_id)
gene_id <- values(sc2_tx)$gene_id[tx_id_match]


###################################################
### code chunk number 7: chipseq-gene-id
###################################################

### check they were all unique - in this case they were
all(elementLengths(gene_id) <= 1)
values(sc2_tx)$gene_id <- drop(gene_id)


###################################################
### code chunk number 8: reconcile-seqlevels
###################################################
sc2_cds_gene <- keepSeqlevels(sc2_cds_gene, "chrXIII")
sc2_cds_gene <- renameSeqlevels(sc2_cds_gene,
                                c(chrXIII = "Scchr13"))


###################################################
### code chunk number 9: countOverlaps
###################################################
sc2_counts <- countOverlaps(sc2_cds_gene, reads_grl, 
                            ignore.strand = TRUE)
### for now this is ANY overlap. not always sufficient for the question being asked

###################################################
### code chunk number 10: countOverlaps-all
###################################################
sc2_counts <- seqapply(bams, function(bam) {
  countOverlaps(sc2_cds_gene, 
                grglist(readGappedAlignments(bam)),
                ignore.strand = TRUE)
})


###################################################
### code chunk number 11: SummarizedExperiment
###################################################
assay <- list(counts = do.call(cbind, as.list(sc2_counts)))
rowData <- sc2_cds_gene
samples <- gsub("_.*", "", basename(bams))
colData <- DataFrame(genotype = gsub(".$", "", samples), 
                     lane = gsub(".*(.)$", "\\1", samples), 
                     row.names = samples)
se <- SummarizedExperiment(assay, rowData, colData)

se   ## has a show method.

?SummarizedExperiment
### useful to store lots of things about expt incl metadata
### Michael really recommends using this format


###################################################
### code chunk number 12: rnaseq-rpkm
###################################################

### some accessors
sc2_count_mat <- assay(se)
k <- sum(width(rowData(se)))
m <- countBam(BamFileList(bams))[,"records"] / 1e6
rpkm <- sc2_count_mat / k / 
  rep(m, each = nrow(sc2_count_mat))

### summarizeOverlaps is a shortcut (discards multiply mapping reads) counts overlpas between a set of BAM files and a set of genomic features (has some default setting should look at)


###################################################
### code chunk number 13: feeding-deseq
###################################################
library(DESeq)
rowDataDf <- as.data.frame(unlist(range(rowData(se))))
cond <- as.data.frame(colData(se))
featureData <- AnnotatedDataFrame(rowDataDf)
deseq <- newCountDataSet(assay(se), cond = cond,
                         featureData = featureData)


###################################################
### code chunk number 14: feeding-edger
###################################################
library(edgeR)
edger <- DGEList(assay(se), group = colData$genotype,
                 genes = rowDataDf)


###################################################
### code chunk number 15: rnaseq-junctions
###################################################
hits <- findOverlaps(reads_grl, sc2_cds_gene, 
                     ignore.strand = TRUE)
values(reads_grl)$hits <- as(hits, "List")
head(table(elementLengths(values(reads_grl)$hits)))


###################################################
### code chunk number 16: aldoa-roi
###################################################
aldoa_eg <- org.Hs.egSYMBOL2EG$ALDOA
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
aldoa_exons <- exons(txdb, vals = list(gene_id = aldoa_eg),
                     columns = c("gene_id", "tx_id"))


###################################################
### code chunk number 17: aldoa-exons-tx
###################################################


aldoa_range <- range(aldoa_exons)
## shows span of those exons


aldoa_exons_tx <- multisplit(aldoa_exons, 
                             values(aldoa_exons)$tx_id)
# for each transcript, show all exons


###################################################
### code chunk number 18: aldoa-cds-tx
###################################################
aldoa_cds <- cds(txdb, vals = list(gene_id = aldoa_eg),
                 columns = c("gene_id", "tx_id"))
aldoa_cds_tx <- multisplit(aldoa_cds, 
                           values(aldoa_cds)$tx_id)
aldoa_exons_tx <- aldoa_exons_tx[names(aldoa_cds_tx)]


###################################################
### code chunk number 19: ccds-only
###################################################


aldoa_exons_tx <- aldoa_exons_tx[c(1, 3, 4, 6)]
aldoa_cds_tx <- aldoa_cds_tx[c(1, 3, 4, 6)]

## CCDS is a higher confidence subset of UCSC known genes


###################################################
### code chunk number 20: aldoa-bams
###################################################
extdatadir <- system.file("extdata", 
                          package = "RangesRNAseqTutorial")
files <- tools::list_files_with_exts(extdatadir, "bam")
names(files) <- tools::file_path_sans_ext(basename(files))
bamFiles <- BamFileList(files)


###################################################
### code chunk number 21: normal-bam
###################################################
bam <- bamFiles$normal


###################################################
### code chunk number 22: read-bam-paired
###################################################
param <- ScanBamParam(tag = "XS", what = "isize", 
                      which = aldoa_range)
ga <- readGappedAlignmentPairs(path(bam), param = param)


temp1 <- left(ga)
temp2 <- GenomicRanges::left(ga)

### but there is a conflict with left from another package (locfit) so need to specify the package

###################################################
### code chunk number 23: isize
###################################################
fraglen <- abs(values(first(ga))$isize)

### not all aligners populate the isize field
### it calculates genomic extent of paired end, removes size of introns in cigar field but does NOT remove size of introns that are in the inter-read gap.  so there will be a subset of reads that are quite overestimated in size (apparently it's easy to identify those as outliers)

## alternatively: an exercise to calculate it, not worrying about the introns

fragsNoIntrons <-  end(GenomicRanges::right(ga)) + 1 - start(GenomicRanges::left(ga))

head(cigar(GenomicRanges::right(ga)))



###################################################
### code chunk number 24: fraglen-solution
###################################################
fraglen_simple <- end(right(ga)) - start(left(ga)) + 1L


###################################################
### code chunk number 25: fraglen-plot
###################################################
calcFragLenCDF <-  function(x) {
  tab <- tabulate(x, max(x))
  DataFrame(fraglen = seq_len(max(x)), cumfrac = cumsum(tab) / sum(tab))
}
plotFragLenCDF <- function(x, xlim, ...) {
  if (missing(xlim))
    xlim <- c(x$fraglen[which(x$cumfrac > 0)[1]], max(x$fraglen))
  p <- qplot(fraglen, cumfrac, data = as.data.frame(x), geom = "line", 
             ylab = "Fraction of reads with fragment length <= X",
             xlab = "Fragment length", xlim = xlim, ...)
  ## if (!is.null(x$sample))
  ##   p <- p + geom_line(aes(group = sample))
  p
}
## fraglen_cdf <- calcFragLenCDF(fraglen)
## plotFragLenCDF(fraglen_cdf, xlim = c(0, 2000))
plotFragLenDensity <- function(fraglen, xlim, ...) {
  p <- qplot(fraglen, geom = "density", 
             xlab = "Fragment length", xlim = xlim, ...)
  p
}
plotFragLenDensity(fraglen, xlim = c(0, 2000))


###################################################
### code chunk number 26: findSpliceOverlaps-normal
###################################################
hits <- findSpliceOverlaps(bamFiles$normal, aldoa_exons_tx, 
                           cds = aldoa_cds_tx, 
                           pairedEnd = TRUE)

hits ### to show

###################################################
### code chunk number 27: findSpliceOverlaps-hidden
###################################################
hits_normal <- hits
hits_tumor <- findSpliceOverlaps(bamFiles$tumor, aldoa_exons_tx, 
                                 cds = aldoa_cds_tx,
                                 pairedEnd = TRUE)


###################################################
### code chunk number 28: count-splice-hits
###################################################
countSpliceHits <- function(hits) {
  coding_hits <- hits[values(hits)$coding]
  compatible_hits <- 
    coding_hits[values(coding_hits)$compatible]
  unique_hits <- coding_hits[values(coding_hits)$unique]
  cbind(coding = countSubjectHits(coding_hits),
        compatible = countSubjectHits(compatible_hits),
        unique = countSubjectHits(unique_hits))
}
counts <- countSpliceHits(hits)

#### chunk 29 was in response to a challenge on slide 53

###################################################
### code chunk number 29: summarizeSpliceOverlaps
###################################################
summarizeSpliceOverlaps <- function(bams, tx, cds) {
  counts <- lapply(bams, function(bam) {
    hits <- findSpliceOverlaps(bam, tx, cds = cds)
    counts <- countSpliceHits(hits)
    split(counts, colnames(counts)[col(counts)])
  })
  assays <- do.call(mapply, 
                    c(cbind, counts, SIMPLIFY = FALSE))
  colData <- DataFrame(tumorStatus = c("tumor", "normal"))
  rownames(colData) <- colData$tumorStatus
  SummarizedExperiment(assays, tx, colData)
}
se <- summarizeSpliceOverlaps(bamFiles, aldoa_exons_tx, 
                              aldoa_cds_tx)


###################################################
### code chunk number 30: order-se
###################################################
uc <- assay(se, "unique")
uc_ord <- order(rowSums(uc), decreasing = TRUE)
uc_top <- uc[head(uc_ord, 2),]
fisher.test(uc_top)$estimate


###################################################
### code chunk number 31: propagateXS
###################################################
propagateXS <- function(ga) {
  first_xs <- values(first(ga))$XS
  last_xs <- values(last(ga))$XS
  xs <- first_xs
  xs[is.na(xs)] <- last_xs[is.na(xs)]
  values(ga)$XS <- xs
  ga
}
#### the per read information should get propagated up to the gappedAlignemetPairs, if you want to use that information

###################################################
### code chunk number 32: loadAlignments
###################################################
loadAlignments <- function(bam) {
  ga <- readGappedAlignmentPairs(path(bam), 
                                 param = param)
  propagateXS(ga)
}
ga_normal <- loadAlignments(bamFiles$normal)
ga_tumor <- loadAlignments(bamFiles$tumor)


###################################################
### code chunk number 33: resolveStrand
###################################################
resolveStrandFromXS <- function(reads, 
                                xs = values(reads)$XS) 
{
  strand <- ifelse(!is.na(xs), xs, "*")
  strand(reads) <- relist(Rle(strand, elementLengths(reads)), 
                          reads)
  reads
}


###################################################
### code chunk number 34: loadReads
###################################################
getReadRanges <- function(ga) {
  resolveStrandFromXS(grglist(ga))
}
reads_normal <- getReadRanges(ga_normal)
reads_tumor <- getReadRanges(ga_tumor)


###################################################
### code chunk number 35: get-unique-reads
###################################################
getUniques <- function(hits) {
  unique(queryHits(hits[values(hits)$unique & 
                        subjectHits(hits) %in% c(1, 4)]))
}
unique_normal <- getUniques(hits_normal)
unique_tumor <- getUniques(hits_tumor)
unique_reads_normal <- reads_normal[unique_normal]
unique_reads_tumor <- reads_tumor[unique_tumor]
unique_reads <- 
  mstack(normal = unlist(unique_reads_normal),
         tumor = unlist(unique_reads_tumor))

### mstack - useful for plotting


###################################################
### code chunk number 36: ggbio-isoform-coverage
###################################################
read_track <- autoplot(unique_reads, stat = "coverage",
                       facets = name ~ .)
read_track ## makes the plot
                       
tx_track <- autoplot(aldoa_exons_tx, geom = "alignment", 
                     ylab = "")
tx_track
                     
cov_tracks <- tracks(read_track, tx_track, 
                     heights = c(3, 1))
cov_tracks ### plots the two together

###autoplot is from ggbio

###################################################
### code chunk number 37: ggbio-zoom
###################################################
xlim(cov_tracks) <- c(30075000, 30080000)
cov_tracks
### manul zoom

###################################################
### code chunk number 38: ggbio-zoom-coverage
###################################################

## or automatic zoom to region of high coverage
cov_chr16 <- coverage(unique_reads)$chr16
roi <- range(ranges(slice(cov_chr16, 1000)))
roi <- roi + 500
xlim(cov_tracks) <- roi
cov_tracks


###################################################
### code chunk number 39: splices
###################################################
splices <- resolveStrandFromXS(introns(ga_normal), 
                               values(ga_normal)$XS)


###################################################
### code chunk number 40: grkey
###################################################
gr2key <- function(x) {
  paste(seqnames(x), start(x), end(x), strand(x), 
        sep = ":")
}
### hack to make identifiers for each intron

###################################################
### code chunk number 41: key2gr
###################################################
key2gr <- function(x, ...) {
  key_mat <- matrix(unlist(strsplit(x, ":", fixed=TRUE)), 
                    nrow = 4)
  GRanges(key_mat[1,],
          IRanges(as.integer(key_mat[2,]), 
                  as.integer(key_mat[3,])),
          key_mat[4,], ...)
}


###################################################
### code chunk number 42: txkey
###################################################

# psetdiff gets gaps on a GRangesList, i.e. the introns

introns <- psetdiff(range(aldoa_exons_tx), aldoa_exons_tx)
introns_flat <- unlist(introns, use.names = FALSE)
tx_keys <- gr2key(introns_flat)


###################################################
### code chunk number 43: splice-summary
###################################################
splices_flat <- unlist(splices, use.names = FALSE)
splice_table <- table(gr2key(splices_flat))
splice_summary <- 
  key2gr(names(splice_table), 
         score = as.integer(splice_table),
         novel = !names(splice_table) %in% tx_keys,
         seqlengths = seqlengths(splices))


###################################################
### code chunk number 44: summarizeSplicing
###################################################
summarizeSplices <- function(ga) {
  splices <- resolveStrandFromXS(introns(ga), values(ga)$XS)
splices_flat <- unlist(splices, use.names = FALSE)
splice_table <- table(gr2key(splices_flat))
splice_summary <- 
  key2gr(names(splice_table), 
         score = as.integer(splice_table),
         novel = !names(splice_table) %in% tx_keys,
         seqlengths = seqlengths(splices))
  splice_summary
}


###################################################
### code chunk number 45: combine-splices
###################################################
unique_splices_normal <- 
  summarizeSplices(ga_normal[unique_normal])
unique_splices_tumor <- 
  summarizeSplices(ga_tumor[unique_tumor])
unique_splices <- mstack(normal = unique_splices_normal,
                         tumor = unique_splices_tumor)


###################################################
### code chunk number 46: ggbio-splices
###################################################
read_track <- autoplot(unique_splices, geom = "arch",
                       aes(height = score), 
                       color = "deepskyblue3",
                       ylab = "coverage",
                       facets = name ~ .) + 
  stat_coverage(unique_reads, facets = name ~ .)
tracks(read_track, tx_track, heights = c(3, 1),
       xlim = roi)


###################################################
### code chunk number 47: ggbio-splices-size
###################################################
read_track <- autoplot(unique_splices, geom = "arch",
                       aes(size = score, 
                           height = width / 5), 
                       color = "deepskyblue3",
                       ylab = "coverage",
                       facets = name ~ .) + 
  stat_coverage(unique_reads, facets = name ~ .)
tracks(read_track, tx_track, heights = c(3, 1),
       xlim = roi)


###################################################
### code chunk number 48: novel-splices
###################################################
splices_normal <- summarizeSplices(ga_normal)
splices_tumor <- summarizeSplices(ga_tumor)
all_splices <- mstack(normal = splices_normal,
                      tumor = splices_tumor)
novel_splices <- 
  all_splices[values(all_splices)$novel &
              values(all_splices)$score >= 9]
unique_novel_splices <- c(unique_splices, novel_splices)
values(unique_novel_splices)$annotated <- 
  !values(unique_novel_splices)$novel


###################################################
### code chunk number 49: ggbio-novel-junctions
###################################################
novel_track <- autoplot(unique_novel_splices, 
                        geom = "arch",
                        aes(size = score, 
                            height = width / 5,
                            color = annotated), 
                        ylab = "coverage",
                        facets = name ~ .) +
  stat_coverage(unique_reads, facets = name ~ .)
tracks(novel_track, tx_track, heights = c(3, 1),
       xlim = roi)

