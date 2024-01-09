library("genoset")
openVignette("genoset")

data(genoset)

genoset.ds = GenoSet(locData = locData.rd, foo = fake.lrr, pData = fake.pData,
annotation = "SNP6", universe = "hg19")


baf.ds = BAFSet(locData = locData.rd, lrr = fake.lrr, baf = fake.baf,
pData = fake.pData, annotation = "SNP6", universe = "hg19")


cn.ds = CNSet(locData = locData.rd, cn = fake.lrr, pData = fake.pData,
annotation = "SNP6", universe = "hg19")

cn.ds

rle.baf.ds = BAFSet(locData = locData.rd, lrr = DataFrame(K = Rle(c(rep(1.5, 300), rep(2.3, 700))), L = Rle(c(rep(3.2, 700), rep(2.1, 300))), M = Rle(rep(1.1, 1000)), row.names = rownames(fake.lrr)), baf = DataFrame(K = Rle(c(rep(0.05, 600), rep(0.5, 400))), L = Rle(c(rep(0, 700), rep(0.5, 300))), M = Rle(rep(1, 1000)), row.names = rownames(fake.baf)), pData = fake.pData, annotation = "SNP6")

?BAFSet

gene.rd = RangedData(ranges = IRanges(start = c(3.5e+07, 1.27e+08), end = c(35500000, 1.29e+08), names = c("HER2", "CMYC")), space = c("chr17", "chr8"))

gene.ds = cn.ds[gene.rd, ]



library(DNAcopy) 

cbs.cna = CNA(cn(cn.ds), chr(cn.ds), pos(cn.ds)) 

cbs.smoothed.CNA = smooth.CNA(cbs.cna) 

cbs.segs = segment(cbs.cna)

###### can take a while, so use multicore. does t-test on windows of probes throughout genome.

assayDataElement(cn.ds, "cn.segs") = runCBS(cn(cn.ds), locData(cn.ds))


library(BSgenome.Hsapiens.UCSC.hg19) 

cn.ds2 = cn.ds 
cn.ds2 = loadGC(cn.ds2, expand = 1e+06, bsgenome = Hsapiens) 
head(locData(cn.ds2)$gc)


genoPlot(cn.ds, 1)  # 1 means sample 1
genoPlot(cn.ds, 1, element = "cn.segs", add = TRUE, col = "red")

genoPlot(cn.ds, 1, chr = "chr12") 
genoPlot(cn.ds, 1, chr = "chr12", element = "cn.segs", add = TRUE, col = "red")


mbaf.data = baf2mbaf(baf(baf.ds), hom.cutoff = 0.9) 

assayDataElement(baf.ds, "mbaf") = mbaf.data

### compress in memory buy converting to Rle
mbaf.data = DataFrame(sapply(colnames(mbaf.data), function(x) { Rle(mbaf.data[, x])  }, USE.NAMES = TRUE, simplify = FALSE)) 

as.numeric(object.size(assayDataElement(baf.ds, "mbaf"))) / as.numeric(object.size(mbaf.data))


assayDataElement(baf.ds, "baf.segs") = runCBS(assayDataElement(baf.ds, "mbaf"), locData(baf.ds))

assayDataElement(baf.ds, "lrr.segs") = runCBS(lrr(baf.ds), locData(baf.ds))


gain.list = lapply(sampleNames(baf.ds), function(sample.name) { as.logical(assayDataElement(baf.ds, "lrr.segs")[, sample.name] >
0.3)}) 

gain.mat = do.call(cbind, gain.list) 

gain.freq = rowMeans(gain.mat, na.rm = TRUE)
