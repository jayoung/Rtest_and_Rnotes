library(BSgenome.Hsapiens.UCSC.hg38)

windows10kb <-  tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg38), 
                           tilewidth=10000, 
                           cut.last.tile.in.chrom=TRUE)

system.time( x <- getSeq(Hsapiens, windows10kb[1:100000]) )
#  user  system elapsed 
# 6.625   1.202   8.157 

system.time( y <- Views(Hsapiens, windows10kb[1:100000]) )
#  user  system elapsed 
# 0.023   0.002   0.025 
system.time( y2 <- as(y, "DNAStringSet" ) )
             
sessionInfo()
# 
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Monterey 12.6
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
# [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] BSgenome.Hsapiens.UCSC.hg38_1.4.4 BSgenome_1.64.0                  
# [3] rtracklayer_1.56.1                Biostrings_2.64.1                
# [5] XVector_0.36.0                    GenomicRanges_1.48.0             
# [7] GenomeInfoDb_1.32.3               IRanges_2.30.1                   
# [9] S4Vectors_0.34.0                  BiocGenerics_0.42.0              
# 
# loaded via a namespace (and not attached):
#     [1] rstudioapi_0.14             zlibbioc_1.42.0             GenomicAlignments_1.32.1   
# [4] BiocParallel_1.30.3         lattice_0.20-45             rjson_0.2.21               
# [7] tools_4.2.1                 grid_4.2.1                  SummarizedExperiment_1.26.1
# [10] parallel_4.2.1              Biobase_2.56.0              matrixStats_0.62.0         
# [13] yaml_2.3.5                  crayon_1.5.2                BiocIO_1.6.0               
# [16] Matrix_1.5-1                GenomeInfoDbData_1.2.8      restfulr_0.0.15            
# [19] bitops_1.0-7                codetools_0.2-18            RCurl_1.98-1.9             
# [22] DelayedArray_0.22.0         compiler_4.2.1              MatrixGenerics_1.8.1       
# [25] Rsamtools_2.12.0            XML_3.99-0.11              


