library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(VariantAnnotation)

## read in test vcf file and sort out seqlevels
fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
vcf <- readVcf(fl, "hg19")
seqlevels(vcf) <- "chr22"

### for my fly SNPs, locateVariants was failing to annotate some SNPs within non-coding transcripts. I wanted to troubleshoot this using a human SNP in a non-coding transcript.  Here, it gives me an annotation of "intergenic" instead of no annotation at all, which is what I got for the fly SNPs. 

##### solution may be to define txdb using cached version, to get around a problem described here https://github.com/Bioconductor/VariantAnnotation/issues/87 and here https://github.com/Bioconductor/VariantAnnotation/issues/76

### it seemed flaky though

cache <- new.env(parent = emptyenv())

# CodingVariants
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
z <- cdsBy(txdb)
z <- z[lengths(z) > 0L]
cache[["cdsbytx"]] <- z

# IntronVariants
z <- intronsByTranscript(txdb)
z <- z[lengths(z) > 0L]
cache[["intbytx"]] <- z

# ThreeUTRVariants 
z <- threeUTRsByTranscript(txdb)
z <- z[lengths(z) > 0L]
cache[["threeUTRbytx"]] <- z

# FiveUTRVariants 
z <- fiveUTRsByTranscript(txdb)
z <- z[lengths(z) > 0L]
cache[["fiveUTRbytx"]] <- z

# IntergenicVariants 
z <- transcriptsBy(txdb, "gene")
z <- z[lengths(z) > 0L]
cache[["txbygene"]]  <- z

# SpliceSiteVariants 
# no need to populate cache further ?

# PromoterVariants 
z <- transcripts(txdb)
z <- z[lengths(z) > 0L]
cache[["tx"]] <- splitAsList(z, seq_len(length(z))) 
names(cache[["tx"]]) <- z$tx_id


#### get exons for ALL transcripts (not just CDS exons)
all_exons <- exonsBy(txdb, by="tx", use.names=TRUE)

### are there some

# https://genome.ucsc.edu/cgi-bin/hgc?hgsid=2328384550_Yt4aI4a3Ve4GUanl2xbsVmrMGE1l&db=hg19&c=chr11&l=94884568&r=94896739&o=94892022&t=94892023&g=dbSnp155Common&i=rs143452847
test_nonCodingRNA_SNP <- GRanges(
    seqnames="chr11",
    ranges=IRanges(start=94892023, width=1),
    REF=DNAStringSet(c("T")),
    ALT=DNAStringSet(c("C")) )

## intergenic
test_nonCodingRNA_SNP_loc <- locateVariants(test_nonCodingRNA_SNP, txdb, AllVariants(), cache = cache)

test_nonCodingRNA_SNP_loc
# GRanges object with 1 range and 9 metadata columns:
#     seqnames    ranges strand |   LOCATION  LOCSTART    LOCEND
# <Rle> <IRanges>  <Rle> |   <factor> <integer> <integer>
#     [1]    chr11  94892023      * | intergenic      <NA>      <NA>
#     QUERYID      TXID         CDSID      GENEID
# <integer> <integer> <IntegerList> <character>
#     [1]         1      <NA>                      <NA>
#     PRECEDEID                FOLLOWID
# <CharacterList>         <CharacterList>
#     [1] 143684,143686,84441,... 10888,143689,154810,...


subsetByOverlaps(cache[["tx"]],test_nonCodingRNA_SNP)
# GRangesList object of length 1:
#     $`42651`
# GRanges object with 1 range and 2 metadata columns:
#     seqnames            ranges strand |     tx_id     tx_name
# <Rle>         <IRanges>  <Rle> | <integer> <character>
#     [1]    chr11 94883703-94892312      + |     42651  uc001pfi.1
# -------
#     seqinfo: 93 sequences (1 circular) from hg19 genome

"42651" %in% names(cache[["cdsbytx"]])
# FALSE

subsetByOverlaps(all_exons,test_nonCodingRNA_SNP)

# GRangesList object of length 1:
#     $uc001pfi.1
# GRanges object with 3 ranges and 3 metadata columns:
#     seqnames            ranges strand |   exon_id   exon_name exon_rank
# <Rle>         <IRanges>  <Rle> | <integer> <character> <integer>
#     [1]    chr11 94883703-94884912      + |    150379        <NA>         1
# [2]    chr11 94890792-94890967      + |    150380        <NA>         2
# [3]    chr11 94891773-94892312      + |    150381        <NA>         3
# -------
#     seqinfo: 93 sequences (1 circular) from hg19 genome


subsetByOverlaps(cache[["txbygene"]],test_nonCodingRNA_SNP)
# GRangesList object of length 0:
#     <0 elements>

# chr11:94892023-94892023

## get rowRanges, as in the vignette
rd <- rowRanges(vcf)

### use just the first 10 snps
test_snps <- rd[1:10]

## show that they're in an intron:
txdb_introns <- intronsByTranscript(txdb, use.names=TRUE)
findOverlaps(test_snps, txdb_introns)
# Hits object with 10 hits and 0 metadata columns:
#     queryHits subjectHits
# <integer>   <integer>
#     [1]         1       75253
# [2]         2       75253
# [3]         3       75253
# [4]         4       75253
# [5]         5       75253
# [6]         6       75253
# [7]         7       75253
# [8]         8       75253
# [9]         9       75253
# [10]        10       75253

#### locateVariants returns no annotation for these SNPs, whether we use AllVariants() or IntronVariants()
test_loc <- locateVariants(test_snps, txdb, AllVariants())
test_loc
# 
# GRanges object with 0 ranges and 9 metadata columns:
#     seqnames    ranges strand | LOCATION  LOCSTART    LOCEND   QUERYID      TXID
# <Rle> <IRanges>  <Rle> | <factor> <integer> <integer> <integer> <integer>
#     CDSID      GENEID       PRECEDEID        FOLLOWID
# <IntegerList> <character> <CharacterList> <CharacterList>
#     -------
#     seqinfo: no sequences

test_loc_2 <- locateVariants(test_snps, txdb, IntronVariants())
test_loc_2
# GRanges object with 0 ranges and 9 metadata columns:
#     seqnames    ranges strand | LOCATION  LOCSTART    LOCEND   QUERYID      TXID
# <Rle> <IRanges>  <Rle> | <factor> <integer> <integer> <integer> <integer>
#     CDSID      GENEID       PRECEDEID        FOLLOWID
# <IntegerList> <character> <CharacterList> <CharacterList>
#     -------
#     seqinfo: no sequences


test_loc_3 <- locateVariants(test_snps, txdb, AllVariants(), cache = cache)
test_loc_3
# GRanges object with 10 ranges and 9 metadata columns:
#     seqnames    ranges strand | LOCATION  LOCSTART    LOCEND   QUERYID        TXID
# <Rle> <IRanges>  <Rle> | <factor> <integer> <integer> <integer> <character>
#     rs7410291    chr22  50300078      - |   intron     10763     10763         1       75253
# rs147922003    chr22  50300086      - |   intron     10755     10755         2       75253
# rs114143073    chr22  50300101      - |   intron     10740     10740         3       75253
# rs141778433    chr22  50300113      - |   intron     10728     10728         4       75253
# rs182170314    chr22  50300166      - |   intron     10675     10675         5       75253
# rs115145310    chr22  50300187      - |   intron     10654     10654         6       75253
# rs186769856    chr22  50300268      - |   intron     10573     10573         7       75253
# rs77627744    chr22  50300346      - |   intron     10495     10495         8       75253
# rs193230365    chr22  50300423      - |   intron     10418     10418         9       75253
# rs9627788    chr22  50300438      - |   intron     10403     10403        10       75253
# CDSID      GENEID       PRECEDEID        FOLLOWID
# <IntegerList> <character> <CharacterList> <CharacterList>
#     rs7410291                     79087                                
# rs147922003                     79087                                
# rs114143073                     79087                                
# rs141778433                     79087                                
# rs182170314                     79087                                
# rs115145310                     79087                                
# rs186769856                     79087                                
# rs77627744                     79087                                
# rs193230365                     79087                                
# rs9627788                     79087                                
# -------
#     seqinfo: 1 sequence from an unspecified genome; no seqlengths
